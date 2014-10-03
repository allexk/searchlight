#!/usr/bin/python

"""
This is a Fabric script that contains two tasks to speed up running experiments.

The first one is "prepare", which flushes the buffers, copies the task config
and starts the server.

The second is "stop", which stops the server and retrieves remote logs from
the hosts. It also retrieves the task config from the first host.
"""

from fabric.api import env, run, sudo, get, put, settings, cd, task, \
    execute, parallel
import os

# default SciDb dir
DEFAULT_SCIDB_DIR = '/opt/scidb/14.3'

# default data prefix dir
DEFAULT_DATA_PREFIX = '/mnt/sdb/scidb'

# default remote dir to put the task to
DEFAULT_REMOTE_TASK_DIR = '~'

# default Searchlight repo dir
DEFAULT_SL_REPO_DIR = 'src/searchlight'

# default sudo password filr
DEFAULT_SUDO_PASSWD_FILE = 'sudo_passwd'

# default file with queries
DEFAULT_QUERY_FILE = 'query.txt'

# default task file to run in a query 
DEFAULT_TASK_FILE = 'sl.json'

# default scidb user to login and run SciDb
DEFAULT_SCIDB_USER = 'scidb'

# Some utility functions

def get_hosts(cluster):
    """
    Returns a list of hosts belonging to the cluster.
    
    The list is specified in the file '<cluster>.hosts'. One host per line.
    """
    hosts = []
    with open('%s.hosts' % cluster, 'r') as hosts_file:
        hosts = [l.strip() for l in hosts_file if l.strip()]
    return hosts

def get_sorted_dirs(path):
    """
    Returns sub-dirs in alphabetical order in the specified path.
    
    It is assumed to be a part of the task running on some host.
    """
    with cd(path):
        dirs = run('ls', quiet=True).split()
        dirs.sort()
    return dirs

def get_sudo_passwd():
    """
    Get sudo password for remote hosts from the local file, if exists.
    """
    passwd = None
    with open(DEFAULT_SUDO_PASSWD_FILE, 'r') as passwd_file:
        passwd = passwd_file.read().strip()
    return passwd

# Sub-tasks ran from other tasks

def stop_scidb(db):
    """ Stops the specified SciDb cluster"""
    with cd(DEFAULT_SCIDB_DIR):
        with settings(warn_only=True):
            run('python bin/scidb.py stop_all %s' % db)

@parallel
def drop_caches():
    """ Drops OS caches on a host """
    sudo('echo 3 > /proc/sys/vm/drop_caches')

def start_scidb(db):
    """ Starts the specified SciDb cluster """
    with cd(DEFAULT_SCIDB_DIR):
        run('python bin/scidb.py start_all %s' % db)

def copy_logs(cluster, log_dir):
    """
    Retrieves *.log files of the current host to the log_dir.
    
    The function tries to retrieve logs for all servers and
    instances on the current host. It puts the logs into
    ever increasing dirs, like 0,1,2,...,. Thus it's important
    to specify hosts in the same order as in the cluster config
    to avoid misinterpretation of results.
    
    Params:
        cluster -- cluster name (to find the logs dir)
        log_dir -- local directory to put the logs in
    """
    cluster_dir = os.path.join(DEFAULT_DATA_PREFIX, cluster)
    server_dirs = get_sorted_dirs(cluster_dir)
    for sd in server_dirs:
        sd = os.path.join(cluster_dir, sd)
        print "Retrieving logs from server path: %s..." % sd
        inst_dirs = get_sorted_dirs(sd)
        for id in inst_dirs:
            id = os.path.join(sd, id)
            local_id = os.path.join(log_dir, str(env.inst_count))
            env.inst_count += 1
            if not os.path.exists(local_id):
                os.mkdir(local_id)
            print "Retrieving logs from %s to %s..." % (id, local_id)
            with cd(id):
                get('*.log', local_id)
                run('rm -f *.log', quiet=True)

def get_remote_task(log_dir, task_file):
    """Retrieves remote task definition from the host"""
    get(os.path.join(DEFAULT_REMOTE_TASK_DIR, task_file), log_dir)

def run_query(query_file):
    """
    Runs AFL/AQL queries from the specified file.
    
    One query per line. Empty lines are ignored.
    """
    with open(query_file, 'r') as f:
        with cd(DEFAULT_SCIDB_DIR):
            for query in f:
                query = query.strip()
                if query:
                    run('bin/iquery -a -q "%s"' % query)

@parallel
def copy_to_remote(remote_dir, local_file):
    "Places the specified file to a remote host in the specified dir."
    put(local_file, remote_dir)

@parallel
def update_src(repo_dir):
    "Updates a remote repo by pulling from git and building the tree."
    with cd(os.path.join(repo_dir, 'build')):
        run('git pull')
        run('cmake .')
        run('make')
        run('sudo make install')

def set_env(cluster):
    """Sets environment before running a task."""
    env.hosts = env.hosts or get_hosts(cluster)
    env.password = env.password or get_sudo_passwd()

# Main tasks

@task
def copy(cluster, file, remote_dir):
    "Copies a file to a remote dir."
    set_env(cluster)
    execute(copy_to_remote, remote_dir, file)

@task
def update_src(cluster):
    "Updates and rebuilds a repository."
    set_env(cluster)
    execute(DEFAULT_SL_REPO_DIR)

@task
def prepare(cluster, task_file=DEFAULT_TASK_FILE):
    """
    Prepares the specified SciDb cluster for running a task.
    
    It does the following:
        1) OS caches are dropped
        2) task file is copied to all servers
        3) SciDb cluster is started from the first server
        
    Params:
        cluster   -- the name of the cluster
        task_file -- file with the JSON task specification
        user      -- user to login and run SciDb from
    """
    set_env(cluster)
    execute(drop_caches)
    execute(copy_to_remote, DEFAULT_REMOTE_TASK_DIR, task_file)
    execute(start_scidb, cluster, hosts=[env.hosts[0]])

@task
def query(cluster, query_file=DEFAULT_QUERY_FILE):
    """
    Runs AFL/AQL queries specified in the query_file on the remote cluster.
    
    The query file is supposed to have a single query per line. Empty lines are
    ignored.
    """
    set_env(cluster)
    execute(run_query, query_file, hosts=[env.hosts[0]])

@task
def stop(cluster, log_dir, task_file=DEFAULT_TASK_FILE):
    """
    Stops the SciDb cluster and retrieves logs to the specified dir
    
    Params:
        cluster   -- the name of the cluster
        log_dir   -- directory to copy logs to
        task_file -- task file name (for copying from the remote)
        user      -- user name for log in
    """
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    set_env(cluster)
    env.inst_count = 0  # total instance counter (user variable)
    execute(stop_scidb, cluster, hosts=[env.hosts[0]])
    execute(copy_logs, cluster, log_dir)
    execute(get_remote_task, log_dir, task_file, hosts=[env.hosts[0]])

@task
def sl(cluster, log_dir, task_file=DEFAULT_TASK_FILE,
       query_file=DEFAULT_QUERY_FILE):
    """
    Prepares SciDb cluster for running the task, runs the queries, stops it and
    retrieves the logs.
    
    Multiple executions of this task with the same parameters will result
    in running the same query multiple times. This can be used to run the same
    query in an experiment.
    """
    # should enforce 'localhost' since sub-tasks might change env.hosts
    execute(prepare, cluster, task_file, hosts='localhost')
    execute(query, cluster, query_file, hosts='localhost')
    execute(stop, cluster, log_dir, task_file, hosts='localhost')
    
@task
def exp(cluster, log_dir, runs, task_file=DEFAULT_TASK_FILE,
        query_file=DEFAULT_QUERY_FILE):
    """
    Runs an experiment, which is a query executed specified number of times.
    """
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    for run in range(int(runs)):
        run_log_dir = os.path.join(log_dir, 'run-%d' % (run + 1))
        # doesn't make any sense to run on any host except localhost
        print 'Running experiment: run=%d' % (run + 1)
        execute(sl, cluster, run_log_dir, task_file, query_file,
                hosts='localhost')
