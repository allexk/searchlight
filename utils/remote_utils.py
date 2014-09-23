#!/usr/bin/python

"""
This is a Fabric script that contains two tasks to speed up running experiments.

The first one is "prepare", which flushes the buffers, copies the task config
and starts the server.

The second is "stop", which stops the server and retrieves remote logs from
the hosts. It also retrieves the task config from the first host.
"""

from fabric.api import env, run, sudo, get, put, settings, cd, task, execute
import os

# default SciDb dir
DEFAULT_SCIDB_DIR = '/opt/scidb/14.3'
# default data prefix dir
DEFAULT_DATA_PREFIX = '/mnt/sdb/scidb'
# default task dir
DEFAULT_TASK_DIR = '~'
# default sudo password filr
DEFAULT_SUDO_PASSWD_FILE = 'sudo_passwd'

# Some utilities

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
        run('python bin/scidb.py stop_all %s' % db)

def drop_caches():
    """ Drops OS caches on a host """
    sudo('echo 3 > /proc/sys/vm/drop_caches')

def start_scidb(db):
    """ Starts the specified SciDb cluster """
    with cd(DEFAULT_SCIDB_DIR):
        run('python bin/scidb.py start_all %s' % db)

def copy_task(task_file):
    """
    Copies the specified task file to a host.
    
    By default, the task file is copied in the user's home.
    """
    put(task_file, DEFAULT_TASK_DIR)
    
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
    get(os.path.join(DEFAULT_TASK_DIR, task_file), log_dir)

# Main tasks

@task()
def prepare(cluster, task_file='sl.json', user='scidb'):
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
    env.hosts = get_hosts(cluster)
    env.user = user
    env.password = get_sudo_passwd()
    execute(drop_caches)
    execute(copy_task, task_file)
    execute(start_scidb, cluster, hosts=[env.hosts[0]])

@task()
def stop(cluster, log_dir, task_file='sl.json', user='scidb'):
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
    env.hosts = get_hosts(cluster)
    env.user = user
    env.inst_count = 0  # total instance counter
    execute(stop_scidb, cluster, hosts=[env.hosts[0]])
    execute(copy_logs, cluster, log_dir)
    execute(get_remote_task, log_dir, task_file, hosts=[env.hosts[0]])
