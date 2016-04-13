#!/usr/bin/python

"""
This script simulates a user-like relaxation process. It takes JSON config and runs a series of queries until
the conditions are met.
"""

import json
import os
import sys
import subprocess
import time
import random
from scidb_http import SciDBConnection

# logging
import logging
logger = logging.getLogger('relax')
logger_stream = logging.StreamHandler()
logger_stream.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
logger.addHandler(logger_stream)
del logger_stream
logger.setLevel(logging.INFO)

# seed the random with current time
random.seed()


class JsonConfig(object):
    """
    Holds relaxation parameters taken from the config file.
    """
    def __init__(self, file_name):
        """Read the JSON config file."""
        with open(file_name, 'r') as f:
            self.js = json.load(f)

    def get_param(self, path):
        """Get param under the specified path.

        The parameter path is specified in the dot-style, e.g., "relax.card". If the path doesn't exist,
        KeyError is raised.

        The path can also be submitted as a sequence of steps, e.g., ['relax', 'card'].
        """
        if isinstance(path, str):
            steps = path.split('.')
        else:
            steps = path
        res = self.js
        for s in steps:
            res = res[s]
        return res

    def set_param(self, path, val):
        """Set param at path to val."""
        param = self.get_param(path[:-1])
        param[path[-1]] = val

    def __getitem__(self, item):
        return self.get_param(item)

    def __str__(self):
        return self.js


class RelaxDegree(object):
    """
    Computes relaxation degree for given result tuples.
    """
    def __init__(self, relax_config):
        """"Init a new RD computer.
        Initialized by the relaxation configuration.
        """
        self.config = relax_config
        self.total_constraints = len(self.config)
        self.dist_w = float(self.config['dist_w'])

    def _violated(self, name, val):
        """Check if the constraint <name> is violated.

        If the constraint is violated returns the normalized distance for the constraint. Otherwise, returns None.
        """
        const = self.config[name]
        l, r = const['valid']
        max_l, max_r = const['relax']
        if val >= l and val <= r:
            return None
        elif val < l:
            return float(l - val) / (l - max_l)
        else:
            # val > r
            return float(val - r) / (max_r - r)

    def rd(self, res):
        """Compute RD of the result tuple.

        Assuming the result is given as a map: function name -> value
        """
        viol_count = 0.0
        max_rd = 0.0
        for func, val in res.items():
            if func in self.config.keys():
                # existing constraint
                viol = self._violated(func, val)
                if viol:
                    max_rd = max(max_rd, viol)
                    viol_count += 1
        return self.dist_w * max_rd + (1 - self.dist_w) * viol_count / self.total_constraints


class Constraint(object):
    """Simulates constraint with values and relaxation."""
    class Parameter(object):
        """Nested class to hold parameters."""
        def __init__(self, name, val, val_limit, step):
            self.name = name
            self.val = val
            self.limit = val_limit
            self.step = step

        def can_step(self):
            """Check if a relaxation step is possible."""
            return self.step > 0 and self.val < self.limit or self.step < 0 and self.val > self.limit

        def do_step(self):
            """Do a relaxation step."""
            self.val += self.step
            if self.step > 0 and self.val > self.limit:
                self.val = self.limit
            elif self.step < 0 and self.val < self.limit:
                self.val = self.limit

        def __str__(self):
            return "%s=%d" % (self.name, self.val)

    def __init__(self, const_config):
        self.id = int(const_config['id'])
        self.params = []
        for p_name, p in const_config['params'].items():
            self.params.append(self.Parameter(p_name, p['start'], p['limit'], p['step']))

    def can_relax(self):
        """Check if relaxing if possible."""
        for p in self.params:
            if p.can_step():
                return True
        return False

    def relax(self):
        """Relax the constraint."""
        for p in self.params:
            p.do_step()

    def dump_to_jc(self, jc):
        """Dump the constraint values to JsonConfig."""
        for p in self.params:
            path = p.name.split('.')
            ref = jc.get_param(path[:-1])
            ref[path[-1]] = p.val

    def __str__(self):
        return 'id=%d,' % self.id + ','.join([str(p) for p in self.params])


class QueryRunner(object):
    """Runs the query and provides results in form of list of tuples."""
    def __init__(self, relax_config):
        # we're going to use only the binary part, so don't have to provide correct values
        self.scidb = SciDBConnection('localhost', 1239, None, None)
        with open(relax_config.get_param('config.query_file'), 'r') as f:
            self.query = f.read().strip()
        self.logger = logging.getLogger('relax.query')
        self.logger.info('Going to query with: %s' % self.query)
        # the first host is the local one
        self.remote_hosts = relax_config.get_param('config.hosts')[1:]
        self.task_dir = relax_config.get_param('config.shared_task_dir')
        self.task_file = relax_config.get_param('config.shared_task_file')

    def _copy_task(self, task_json):
        """Save JSON task to a file and copy to hosts."""
        local_task_path = os.path.join(self.task_dir, self.task_file)
        # dunp the JSON task locally
        with open(local_task_path, 'wb') as f:
            self.logger.info('Dumping task to ' + local_task_path)
            json.dump(task_json, f, sort_keys=True, indent=4)
        # copy task remotely via scp
        scp_args = ['scp', local_task_path]
        for host in self.remote_hosts:
            self.logger.info('Uploading task file to ' + host)
            self.logger.info(subprocess.check_output(scp_args + ['%s:%s' % (host, self.task_dir)],
                                                     stderr=subprocess.STDOUT))

    @staticmethod
    def _str_to_map(res_str):
        res = {}
        for r in res_str.strip().split(', '):
            (name, val) = r.split('=')
            res[name] = int(val)
        return res

    def run_query(self, task_json):
        # First, set the task file
        self._copy_task(task_json)
        # Then, query
        self.logger.info('Querying SciDB...')
        query_array = iter(self.scidb.query_interactive(self.query))
        result_list = []
        try:
            while True:
                next_res = query_array.next()
                result_list.append((time.time(), self._str_to_map(next_res)))
        except StopIteration:
            # no more results
            pass
        self.logger.info('SciDB query finished')
        self.scidb.close()
        return result_list


class Relaxator(object):
    """Keeps information about the relaxation and computes the next step."""
    def __init__(self, config):
        """Create a new relaxator.

        Expecting RelaxConfig as config.
        """
        self.config = config
        self.logger = logging.getLogger('relax.relaxator')
        query_file_name = self.config.get_param('config.task_file')
        self.query = JsonConfig(query_file_name)
        self.method = self.config.get_param('config.method')
        if self.method not in ['all', 'rr', 'random']:
            raise ValueError('Unknown relaxation method ' + self.method)
        self.target_card = self.config.get_param('config.card')
        self.constraints = []
        for c in self.config['constraints']:
            self.constraints.append(Constraint(c))
        self.query_res = []  # tuples: (start_time, end_time, JSON task, all results)
        self.query_runner = QueryRunner(config)
        self.last_relaxed_ind = -1  # index of the last relaxed constraint (for round-robin)
        self.rd_compute = RelaxDegree(self.config['rd'])

    def can_relax(self):
        """Check if we can relax further."""
        for c in self.constraints:
            if c.can_relax():
                return True
        return False

    def relax(self):
        """Does single relax step."""
        if self.can_relax():
            if self.method == 'all':
                # must relax all constraints
                for c in self.constraints:
                    if c.can_relax():
                        c.relax()
            elif self.method == 'rr':
                # round-robin
                self.last_relaxed_ind += 1
                if self.last_relaxed_ind == len(self.constraints):
                    self.last_relaxed_ind = 0
                # we definitely can_relax(), so we have a constraint...
                while not self.constraints[self.last_relaxed_ind].can_relax():
                    self.last_relaxed_ind += 1
                self.constraints[self.last_relaxed_ind].relax()
            else:
                # random
                # to avoid re-generating random when we have a small number of constraints left, copy and filter
                const_can_relax = [c for c in self.constraints if c.can_relax()]
                const = random.choice(const_can_relax)
                const.relax()

    def query_db(self):
        """Query with the current relaxation."""
        for c in self.constraints:
            c.dump_to_jc(self.query)
        start_time = time.time()
        self.logger.info('Querying with constraints: %s' % '|'.join([str(c) for c in self.constraints]))
        query_res = self.query_runner.run_query(str(self.query))
        # query_res consists of tuples: (time, result)
        end_time = time.time()
        rds = [self.rd_compute.rd(r[1]) for r in query_res]
        self.query_res.append({'task': self.query, 'start': start_time, 'end': end_time, 'res': query_res, 'rds': rds})

    def goal_ok(self):
        """Check if we have reached the target."""
        return len(self.query_res) and len(self.query_res[-1]['res']) > self.target_card

    def stats(self):
        """Return stats about relaxation runs."""
        if len(self.query_res) == 0:
            return None
        relax_start = self.query_res[0]['start']
        time_first = -1.0
        total_time = 0.0
        rd_times = {}
        for qr in self.query_res:
            total_time += qr['end'] - qr['start']
            if time_first == -1.0 and len(qr['res']) > 0:
                time_first = qr['res'][0][0] - relax_start
            for i, rd in enumerate(qr['rd']):
                if rd not in rd_times:
                    rd_times[rd] = qr['res'][i] [0]- relax_start
        return {'total_time': total_time, 'ttf': time_first, 'rd': rd_times}

    def dump_res(self, f):
        for qr in self.query_res:
            f.write('Run (time=%.3f, res_num=%d):\n' % (qr['end'] - qr['start'], len(qr['res'])))
            for i, r in enumerate(qr['res']):
                r_str = ','.join(['%s=%s' % (k, str(v)) for k, v in r.items()])
                r_str += 'rd=%.3f' % qr['rds'][i][1]
                f.write(r_str)
                f.write('\n\nTask:\n')
                f.write(str(qr['task']))
                f.write('\n----------------\n')
        f.write('\nStats follows:\n')
        f.write(str(self.stats()))


def _main(config_file):
    config = JsonConfig(config_file)
    relax = Relaxator(config)
    relax.query_db()
    while not relax.goal_ok() and relax.can_relax():
        relax.relax()
        relax.query_db()
    if not relax.goal_ok():
        logger.info('Goal not completed')
    else:
        logger.info('Goal completed')
    relax.dump_res(sys.stdout)
    return 0

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print 'Usage: relax.py <config_file>'
        sys.exit(1)
    sys.exit(_main(sys.argv[1]))
