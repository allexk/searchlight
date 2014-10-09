#!/usr/bin/python

import os
import datetime
import re
import sys

# String patterns for parsing
START_TIME_RE = re.compile('([\d-]+) ([\d:,]+).*Registered query')
SOLVER_TIME_RE = re.compile('([\d-]+) ([\d:,]+).*Solver total time: ([\d.]+)s')
SEARCH_END_RE = re.compile('([\d-]+) ([\d:,]+).*Broadcasting end-of-search...')
VAL_TOTAL_RE = re.compile('Adapter\(validator\).*time: ([\d.]+)s')
VAL_TOTAL_HELPER_RE = re.compile('Adapter\(validator helper_\d+\).*time: ([\d.]+)s')
VAL_TOTAL_CPU_RE = re.compile('Total aggregates CPU time: ([\d.]+)s')
SEARCH_STATS_RE = re.compile('Main search stats: fails=(\d+), true fails=(\d+),'
    ' candidates=(\d+)')
MSG_FORW_RE = 'Forwards sent \(total\)=([\d\s,]+),'
RES_TIME_RE = re.compile('([\d-]+) ([\d:,]+): .+')

class ParseDict(dict):
    """
    Dictionary for parsing SL logs. It has two functions for parsing lines
    belonging to SL and SL results logs, It tries all matching functions
    in order and updates itself with new information.
    """
    def _parse_start_time(self, l):
        """Parses start time."""
        m = re.search(START_TIME_RE, l)
        if m:
            self['start_time'] = self._parse_timepoint(m.group(1), m.group(2))
    
    def _parse_end_time(self, l):
        """Parses end time."""
        # can use SOLVER_TIME_RE for that, since the times match
        m = re.search(SOLVER_TIME_RE, l)
        if m:
            self['end_time'] = self._parse_timepoint(m.group(1), m.group(2))
                
    def _parse_solver_time(self, l):
        """Parses total solver time."""
        m = re.search(SOLVER_TIME_RE, l)
        if m:
            self['solver_time'] = self._parse_duration(m.group(3))
        
    def _parse_search_time(self, l):
        """Parses total search time."""
        m = re.search(SEARCH_END_RE, l)
        if m:
            self['search_time'] = self._parse_timepoint(m.group(1), m.group(2))
            
    def _parse_validator_times(self, l):
        """Parser validator time (inc. helpers)."""
        if not 'val_times' in self:
            self['val_times'] = []
        m = re.search(VAL_TOTAL_RE, l)
        if m:
            t = self._parse_duration(m.group(1))
            # main validator time is always at the beginning
            self['val_times'].insert(0, t)
        m = re.search(VAL_TOTAL_HELPER_RE, l)
        if m:
            t = self._parse_duration(m.group(1))
            self['val_times'].append(t)
            
    def _parse_candidates(self, l):
        """Parses candidate information."""
        m = re.search(SEARCH_STATS_RE, l)
        if m:
            self['cand_produced'] = int(m.group(3))
        
        m = re.search(MSG_FORW_RE, l)
        if m:
            self['cand_forwarded'] = [int(c) for c in m.group(1).split(', ')]

    def _parse_result_time(self, l):
        """Parses time of a single result."""
        m = re.search(RES_TIME_RE, l)
        if m:
            if not 'res_times' in self:
                self['res_times'] = []
            dt = self._parse_timepoint(m.group(1), m.group(2))
            self['res_times'].append(dt)

    def _parse_timepoint(self, date, time):
        """Parses time out of the date and time strings.
        
        Returns datetime object. Assumes first two space-separated fields are
        date and time, e.g., 2014-09-24 00:04:46,293
        """
        return datetime.datetime.strptime(' '.join([date, time]),
                                          '%Y-%m-%d %H:%M:%S,%f')
    
    def _parse_duration(self, dur):
        """Parses duration in seconds."""
        return datetime.timedelta(seconds=float(dur))
    
    def parse_line(self, l):
        """Tries to parse some info from a line."""
        self._parse_start_time(l)
        self._parse_end_time(l)
        self._parse_solver_time(l)
        self._parse_search_time(l)
        self._parse_validator_times(l)
        self._parse_candidates(l)
        self._parse_result_time(l)

def process_run_dir(dir):
    """ Process a run dir from a Searchlight query and returns some stats.
        
        The stats dictionary contains:
            total_time  -- wall-clock time of the query
            search_time -- wall-clock time of search (sans remaining validation)
            solver_times -- individual solver times
            validator_times -- individual validator times
                (list of lists: [val time, helper times...])
            cands_produced -- number of candidates produced
            cands_forwarded -- number of candidates forwarded
            cands_validated -- number of candidates validated locally
    """
    # determine instance dirs
    dirs = [d for d in os.listdir(dir) if os.path.isdir(os.path.join(dir, d))]
    inst_stats = []
    for d in dirs:
        inst = int(d)
        print 'Parsing logs for instance %d' % inst
        stats = ParseDict()
        rd = os.path.join(dir, d)
        with open(os.path.join(rd, 'searchlight.log'), 'r') as sl_log:
            for l in sl_log:
                stats.parse_line(l)
        if inst == 0:
            # results come only from the coordinator
            with open(os.path.join(rd,
                                   'searchlight_results.log'), 'r') as slr_log:
                for l in slr_log:
                    stats.parse_line(l)
        inst_stats.insert(inst, stats)

    # we should merge stats from instances together
    res = {}
    coord_stats = inst_stats[0]
    res['total_time'] = coord_stats['end_time'] - coord_stats['start_time']
    res['search_time'] = coord_stats['search_time'] - coord_stats['start_time']
    
    res['solver_times'] = [st['solver_time'] for st in inst_stats]
    res['validator_times'] = [st['val_times'] for st in inst_stats]
    
    # candidates require more attention...
    res['cands_produced'] = [st['cand_produced'] for st in inst_stats]
    res['cands_forwarded'] = [sum(st['cand_forwarded']) for st in inst_stats]
    res['cands_validated'] = [t[0] - t[1] for t in zip(res['cands_produced'],
                                                       res['cands_forwarded'])]
    # result durations from the coordinator
    if 'res_times' in coord_stats:
        times = coord_stats['res_times']
        delays = [times[i] - times[i - 1] for i in range(1, len(times))]
        res['first_delay'] = times[0] - coord_stats['start_time']
        delays.insert(0, res['first_delay'])

        # don't count zero delays (i.e., count batches of results)
        delays = [d for d in delays if not d == datetime.timedelta()]

        res['min_delay'] = min(delays)
        res['max_delay'] = max(delays)
        res['avg_delay'] = sum(delays, datetime.timedelta()) / len(delays)
    else:
        res['first_delay'] = 0
        res['min_delay'] = 0
        res['max_delay'] = 0
        res['avg_delay'] = 0

    for st in inst_stats:
        forw = st['cand_forwarded']
        for i in range(len(forw)):
            res['cands_validated'][i] += forw[i]
    
    return res

def avg_stat(stats, key, init=0):
    """Averages a stat with the key. Initial sim value is specified as init."""
    l = [s[key] for s in stats]
    return sum(l, init) / len(l)

def avg_list_stat(stats, key, init=0, sub_ind=None):
    """Averages a list statistic with the key.
    
    Basically, element-wise avg(). If the stat is a list, only the first
    element will be averaged.
    """
    res = []
    for i in range(len(stats[0][key])):
        l = [rs[key][i] for rs in stats]
        if isinstance(l[0], list):
            l = [sl[0] for sl in l]
        res.append(sum(l, init) / len(l))
    return res

def parse_experiment(dir):
    """
    Parses logs for each run in the specified directory. Some stats from
    different runs are averaged.
    """
    print 'Processing experiment:', dir
    if not os.path.exists(os.path.join(dir, 'run-1')):
        # no individual runs
        print 'Processing single-run experiment...'
        return process_run_dir(dir)

    # we're here if we have multiple runs
    run_dirs = [d for d in os.listdir(dir)
                if os.path.isdir(os.path.join(dir, d)) and d.startswith('run-')]
    run_dirs.sort()
    runs_stats = []
    for rd in run_dirs:
        print '\nProcessing run:', rd
        rs = process_run_dir(os.path.join(dir, rd))
        print '\nStats for run:', rd
        print '--------------------'
        print_stats(rs)
        runs_stats.append(rs)

    # merge
    res = {}
    res['total_time'] = avg_stat(runs_stats, 'total_time',
                                 datetime.timedelta())
    res['search_time'] = avg_stat(runs_stats, 'search_time',
                                  datetime.timedelta())
    res['solver_times'] = avg_list_stat(runs_stats, 'solver_times',
                                        datetime.timedelta())
    res['validator_times'] = avg_list_stat(runs_stats, 'validator_times',
                                           datetime.timedelta())
    res['cands_produced'] = avg_list_stat(runs_stats, 'cands_produced')
    res['cands_forwarded'] = avg_list_stat(runs_stats, 'cands_forwarded')
    res['cands_validated'] = avg_list_stat(runs_stats, 'cands_validated')
    
    res['solver_times'] = avg_list_stat(runs_stats, 'solver_times',
                                        datetime.timedelta())

    res['first_delay'] = avg_stat(runs_stats, 'first_delay',
                                           datetime.timedelta())
    res['min_delay'] = avg_stat(runs_stats, 'min_delay',
                                           datetime.timedelta())
    res['max_delay'] = avg_stat(runs_stats, 'max_delay',
                                           datetime.timedelta())
    res['avg_delay'] = avg_stat(runs_stats, 'avg_delay',
                                           datetime.timedelta())
    return res

def print_stats(stats):
    """Pretty-prints the SL statistics."""
    print 'Total time:           ', stats['total_time']
    print 'Search time:          ', stats['search_time']
    print 'Solver times:         ', [str(t) for t in  stats['solver_times']]
    print 'Validator times:      ',
    for vt in stats['validator_times']:
        if isinstance(vt, list):
            print [str(t) for t in vt],
        else:
            print vt,
    print
    print 'Candidates produced:  ', stats['cands_produced']
    print 'Candidates forwarded: ', stats['cands_forwarded']
    print 'Candidates validated: ', stats['cands_validated']

    print
    print 'First result delay: ', stats['first_delay']
    print 'Minimum result delay: ', stats['min_delay']
    print 'Maximum result delay: ', stats['max_delay']
    print 'Average result delay: ', stats['avg_delay']

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: stats.py <experiment dir>'
        sys.exit(1)

    if not os.path.isdir(sys.argv[1]):
        print 'Cannot find directory:', sys.argv[1]

    stats = parse_experiment(sys.argv[1])
    print '\nAverage stats:'
    print_stats(stats)
    sys.exit(0)
