# Searchlight example

To run Searchlight we need:
1. Start a SciDB cluster.
2. Generate some data
3. Load the data
4. Run a simple query

## SciDB Cluster
In general we want a bunch of servers and a single instance per server. This can be achieved by the following config.ini section (config.ini is
located in the etc directory of SciDb installation):
```
[sl]
server-0=hades,0
server-1=artemis,1
db_user=sl
db_passwd=sl
install_root=/opt/scidb/14.12
pluginsdir=/opt/scidb/14.12/lib/scidb/plugins
logconf=/opt/scidb/14.12/share/scidb/log4cxx.properties
base-path=/data/sl
base-port=1239
interface=eth0
smgr-cache-size=1024
mem-array-threshold=2048
```

Where "hades" and "artemis" are machine names. smgr-cache-size is the size of the cache for the _local_ portion of the data array. 
mem-array-threshold is the size of the cache for the chunks transferred from _remote_ nodes. This one should be usually larger than
for the local, if the deafault array distribution is used. SciDB uses round-robin chunk distribution. However, validators divide array
in continuous stripes. That means if, for instance, the number of valdators is 4 on a 4 node cluster, each validator will have 1/4 local
chunks and 3/4 remote chunks, so mem-array-threshold should be aaproximately 3 times larger than smgr-cache-size. Validators don't necessarily
read the entire array, so these params require some tuning. On real experiments, they must be set according to the amount on memory
available! Otherwise, the memory is wasted.

The config.ini must be the same at every SciDB node and should reside at the same location.

Additionally, we'll need a running PostgreSQL server. The server should be open to the world, since other nodes will try to connect to
it. The server will contain metadata for the entire SciDB cluster in the database named the same as the cluster (i.e., "sl" in our example).
PostgreSQL neeeds to run only on one node, which will be the initiator.

SciDB has a utility called `scidb-prepare-db.sh` in its binary dir. It will bootstrap the metadata into the PostgreSQL. The db user and
password must match the `db_user` and `db_passwd` fields in the config.ini. The database name must be the same as the cluster name,
which equals to the section name in the config.ini (e.g., `[sl]`).

The SciDB cluster itself can be started by running:
```
python scidb.py initall sl
python scidb.py startall sl
```
where "sl" is the cluster name. scidb.py is located in the SciDB binary directory.

If the start is successful, a command like `iquery -a -q "list('instances')"` should output running instances.

## Generating simple data
The main Searchlight directory contains a number of utils in the `utils` dir. One of those is a toy data generator, which
can be run as `python generators/synth_gen.py`. Let's create a simple 10x10 array:

```
python utils/generators/synth_gen.py --lbs 0 0 --rbs 9 9 --density 0.8 --sgrid 10 10 --mean 50 51 1 -- toy_data.csv
```
lbs and rbs specify the boundaries. sgrid and lgrid define the data configuration. The array is divided into "clusters"
with the steps defined by lgrid. Each cluster is divided into a number of "cells", with the steps defined in sgrid.
Each cell will contain data generated with normal distribution with mean and dev params. mean is taken from the array of means,
which is specified by the `mean` param. The `mean l r s` creates an array of means from l to r (excluded) with step s. The mean
values is reset randomly for every cell. While the density is reset for each cluster.

After the toy data is generated, it can be loaded into SciDB. For this we need to create an array and transform CSV to the
SciDB format. First, the array. The iquery AQL command looks like this:
```
create array t<val: double>[x=0:9,10,0, y=0:9,10,0]
```
which should be self-explanatory. 0:9,10,0 means the indexes go from 0 to 9, with SciDB chunks of size 10 (so onle a single chunk here) and
0 overlap between the chunks.

The next step is to transform the data to SciDB format. For this we can use the csv2scidb located in the SciDB binary dir:
```
/opt/scidb/14.12/bin/csv2scidb -i toy_data.csv -o toy_data.dat -c 1000 -f 0 -p NNN -s 1
```
parameters can be seen with the "--help" option. But the basic idea is that the command creates a 1D array starting from (0, chunk size 1000,
format (Number, Number, Number), skipping first line (the CSV header).

Now we can load the array via a proxy 1D array:
```
create array traw<x:int64,y:int64,val:double>[i=0:*,1000,0];
load(traw, '/searchlight/toy_data.dat');
insert(redimension(traw, t), t);
remove(traw);
```
The redimension command here does the trick by converting the 1D to 2D array.

The next step in running Searchlight is to create at least one synopsis of the data, which is stored in another array with
predefined name. There is a utility called "sampler.py" in the utils Searchlight dir. The synopsis is described in the paper, but what
the command does is divides the data into logical cells and computes the aggregate values for each cell. Then it can dump the result
in the binary SciDB format and create a scipt for loading it. The command below creates a synopsis with cells of size 5x5 for the specified
region of the array ((0,0)-(9,9); the region should cover the whole array):
```
PYTHONPATH=/searchlight/utils/scidb_query:/searchlight/utils/mimic:$PYTHONPATH python utils/sampler/sampler.py --bindir /opt/scidb/14.12/bin/ --chunks 5 5 --outfile toy_data.sample --outscript toy_data.sh --region 0 0 9 9 --port 1239 --binary t val
```
NB: 1239 is the default iquery port.

We need a couple of additional modules, hence the PYTHONPATCH change. The main params from the commands are "chunks", which defines the steps of the synopsis grid, and
"region", which defines the region to sample.

We then can load the sample into SciDB by using the generated script:
```
bash toy_data.sh
```
NB: right now the scipt only works when the output format is binary (the "--binary" option).

## Running Searchlight
To run Searchlight first we need to load the library via the following AFL command:
```
load_library('searchlight_udo')
```
This has to be done only once after cluster init. Although, if the library ever updated, the cluster should be restarted.

Searchlight itself is run by using "searchlight()" operator, which requires to specify a number of parameters:

1. The array to search over. Will be 't' in our example.
2. Synopsis arrays (might be several of them). Will be 't_val_5x5' in our example. Note, that the synopsis names follow a specific format:
	`<array_name>_<attribute_name>_<dim1xdim2x...xdimN>` for a N-dimensional synopsis.
3. The query specification. Query specification has the following format: `<task_name>:<query_name>:<path_to_config>`, where:
	* Task is a collection of queries. This is just a logical coupling, but in fact it corresponds to a shared (SO) library. In our example it will be "searchlight_sw",
		which corresponds to a collection of SemanticWindows (SW) tasks.
	* Query is taken from the corresponsing task. This corresponds to a function in the libary. Will be "SemWindowsAvg" in our example. It searches for windows with
		specified constraints.
	* Config is a JSON-formatted file, that contains both query-specific and global Searchlight params. Will be "toy.json" (attached as a separate file).

The command looks like this:
```
searchlight(depart(t), depart(t_val_5x5), 'searchlight_sw:SemWindowsAvg:/searchlight/toy.json')
```
Note: the JSON config file must be available at the specified location on _every_ SciDB node in the cluster.
Note: each array including synopsis must be accessed with the "depart()" operator. This is a Searchlght-speicic operator facilitating array chunks redistribution
	and caching needed for Searchlight, since Searchlight generally requires different distribution than SciDB provides.

As can be seen after running the command, no results has been found. You can change the values of "avg_l" and "avg_u" to 50 and 60 correspondingly, which should find
a bunch of regions, but it's more interesting to see relaxation in action. Set relax.on to 1 and run again. As you can see, we now get much more results, with relaxation
ranks attached. Note, Searchlight outputs the _running_ results, which contain intermediate ones (i.e., with worse ratings than the final results).

## Appendix A: toy.json query config file
```
{
    "searchlight": {
        "sl": {
            "max_probes": 100,
            "probes_percentage": 0.1,
            "splits": 100,
            "probe_submit": 0,
            "interval_penalty": 0.5,

            "time_strategy": "fin",
            "time_interval": 300,
            "restarts": 0,

            "luby_scale": 0,

            "fails_restart_probes": 1000,
            "fails_restart_thr": -1,

            "val_sync": 0
        },

        "validator": {
            "max_validations": 100000000,
            "restart_period": 1024,
            "send_info_period": 100,
            "max_helpers": 0,
            "helper_workload": 10,
            "zones": 1,
            "low_watermark": 100,
            "high_watermark": 200,
            "rebal_watermark": 1000,
            "sort": "rd"
        },

        "sampler": {
            "preload": 0,
            "cache": "eager",
            "memory_per_attr": 1024,
            "cell_thr": 0.75,
            "mbr_thr": 0.75,
            "cell_limit": 1000
        },

        "load_aux_samples": 0,
        "dynamic_scheduling": 1
    },

    "balance": {
        "solver_balance": 1,
        "validator_balance": "stripes",
        "map_update_frequency": 0,
        "general_low": 0.05,
        "general_high": 0.5,
        "slices_number": 100,
        "load_slices": 1
    },

    "setup": {
        "solvers": {"0": 2, "1": 2, "2": 2, "3": 2},
        "validators": [0, 1, 2, 3]
    },

    "relax": {
        "on": 0,
        "contr_type": 0,
        "heur": "guess",
        "card": 50,
        "replay": "viol",
        "force_replays": 0,
        "sort": "best",
        "save_udfs": 1,
        "replay_rd": 0.3,
        "spec": 0
    },

    "sw": {
        "lx": 0,
        "ux": 9,
        "ly": 0,
        "uy": 9,
        "avg_l": 60,
        "avg_u": 70,
        "avg_relax_l": 50,
        "avg_relax_u": 80,
        "avg_contr": 1,
        "size_l": 9,
        "size_u": 25,
        "len_lx": 3,
        "len_ux": 5,
        "len_ly": 3,
        "len_uy": 5,
        "step_x": 1,
        "step_y": 1,
        "step_lx": 1,
        "step_ly": 1,

        "neighborhood": {
            "size": 0,
            "max_diff": 1,
            "max_relax_diff": 1,
            "max_diff_high": 100,
            "max_diff_contr": 1,

            "min_diff": 1,
            "min_relax_diff": 1,
            "min_diff_high": 100,
            "min_diff_contr": 1
        },

        "db": "split",
        "time_limit": 3600
    }
}
```

## Appendix B: Searchlight logging

Searchlight uses the same logging mechanism SciDB uses: log4cxx. The following lines can be added to the log config file (share/scidb/log4cxx.properties by default):
```
log4j.logger.searchlight=DEBUG, sl
log4j.additivity.searchlight=false

log4j.appender.sl=org.apache.log4j.RollingFileAppender
log4j.appender.sl.File=searchlight.log
log4j.appender.sl.MaxFileSize=10000KB
log4j.appender.sl.MaxBackupIndex=5
log4j.appender.sl.layout=org.apache.log4j.PatternLayout
log4j.appender.sl.layout.ConversionPattern=%d [%t] %c %-5p: %m%n

log4j.logger.searchlight.result=TRACE, results
log4j.additivity.searchlight.result=false

log4j.appender.results=org.apache.log4j.FileAppender
log4j.appender.results.File=searchlight_results.log
log4j.appender.results.Append=false
log4j.appender.results.layout=org.apache.log4j.PatternLayout
log4j.appender.results.layout.ConversionPattern=%d: %m%n
```

This will dump logging information (TRACE level) to the searchligh.log file in the SciDB directory (alongside the scidb.log file). The results the node finds
are dumped at searchlight_results.log.
