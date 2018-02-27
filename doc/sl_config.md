# Searchlight Configuration
So, the config params/args for a query generally correspond to the attached JSON. It consists of several sections (JSON maps).

## The "searchlight" section
"searchlight" section contains the general Searchlght params. Should be pretty much the same between tasks/queries. It has several sub-maps: "sl" contains config for the solvers,
"validator" is for the validator, "sampler" for the sampler. sl and validator correspond to the main components you probably know from the paper. Solvers (sl) do
the search, while validators check the possible solutions (candidates) over real data. The number of Solver threads is configured for each query, but the number of validator
threads is determined dynamically during execution. Sampler, on the other hand is an internal component that handles all the synopses of the data (the synopses are
described in the paper). Each synopsis is contained in a sample array, which you specify when you run the searchlight() statement. So I use synopses and samples
interchangeably.

### The "searchligh.sl" section
This one contains general Solver params impacting the search heristic. The search heuristic (known as Decision Builder in or-tools) is
set in the task-dependent section (see below), but the general values are:

* "impact". This is the default or-tools strategy. It's not that useful. I don't remeber what they use there, to be honest. Can dig it
	up from the sources. It probably has some init phase, where they try to assess which intervals worth exploring based on the number and
	frequency of fails. Or something like that.
* "random". At every search step (search tree node) simply assign a random value to a random var. This is also from or-tools.
* "split". At every step choose the var with the largest _current_ domain and split it in half. Go to the lowest half first. I
	emphasize "current" here, since var domains change (reduce) when the search goes deeper, which is obvious. This strategy supports
	solver balancing described above. I used it as a default most of the time, since it works very well.
* "sl". So this is our custom heuristic, which has two stages. I'll describe it in the next paragraph.

The custom "sl" heuristic is our attempt to play with heuristics. The main idea is that if your search
space is large, you can devide it into hypercubes (e.g., squares for a 2D-matrix). We do "searchligh.sl.splits" (this is the param name)
intervals along each dimensions (e.g., splits=2 will result in 4 squares in a 2D-matrix). Then we probe each square by trying var assignments
from it. The total number of assignments is obviously the product of square's dimensions, say X. But we restrict the number of probes to
be the max("searchligh.sl.probes_percentage" times X, "searchligh.sl.max_probes"). If "searchligh.sl.probe_submit" is 1,
we submit successful probes to the validator immediately in hopes of a quick answer. The percentage of
successfull probes determines which sqaures will explore first (e.g., zoom-in into). After zooming-in into
a square, we just run the random heuristic there, trying random assignments. There's a bunch of config options there, like how long to run it
"searchligh.sl.time_strategy"/"time_interval", if we want to restart search from time to time (searchligh.sl.restarts; common techinque in random search),
and if we do what Luby scale we want to use (searchligh.sl.luby_scale; this is common restart strategy in search as well). The latter also
aplies to the "random" heuristic. Right now we don't need to use this probably, but if in the future you
want to play with it, we can discuss that. These params in the "searchlight.sl" section should not have any impact on anything else.

searchligh.sl.val_sync basically stops the solver in case there are more than searchligh.validator.max_validastions pending
candidates. Should be always 0, since dynamic balancing (described below) works better.

### The "searchligh.validator" section
This one controls the validator:

* "max_validations" defines the number of pending validations. If the number is larger, the solver will block. Solvers usually much faster than validators, but this should be usually set to a high number.
* "max_period" is a memory optimization. Validator uses an or-tools solver to validate solution, but makes it read the data intead of the synopsis. The
	problem is that it creates a left deep tree. So, we restart it every "max_period" assignments to clean up the or-tools solver. The default is usually fine.
* "max_helpers" is the maximum number of additonal validator threads (there is always one). Since I/O is expensive this is not to waste time waiting on it
	and do the CPU part of validation in parallel.
	Validator determines the need in such threads automatically, but this sets the maximum. Note, if the dynamic scheduling is turned on, this param
	is ignored, since we balance the number of threads between the solver and validator. I'll write about this below.
* "helper_workload". Number of solutions each validator thread takes at a time. It's bacthed to avoid queue thrashing.
* "zones". This is to try some locality optimization during validations. As described in the paper, a validator "owns" a stripe of the data
	(e.g., a slice in a 2D array). A slice usually doesn't fit in memory. So, we can try to divied the slice into zones. Then
	we try to stick to the same zones, if possible. The idea is that a zone or two would fit into main memory and will avoid
	disk I/O. So, it's similar to caching, in the sense that it provides space and time locality optimization. But if cancidates
	come from different zones, that might result in thrashing, where we constantly load and evict array blocks from memory.
	So, "zones=0" actually computes the number of zones automatically based on the amount of memory available and assuming we want
	2 zones in memory. "zones=1" should disable zones all together.
* "sort". How the candidates are sorted in the validator queue. "rd" means they are sorted on the relaxation degree (described
	in the paper and relaxation section here). If there is no relaxation active, then it is basically FIFO. "cost" means we
	try to assess the number of chunks required to validate the candidate and sort based on this cost.
* "send_info_period". This is for balancing validators. If the local validator is really overloaded, then we can offload the
	some candidates to a remote validator. To make it more efficient, validators periodically broadcast some information
	about the number of candidates that are pending. So this param says when to broadcast this info, e.g., 100 means
	when the pending candidates queue size changes by >= 100.
* "low/high_watermark". So, this is for the validator balancing. The idea here is that we have threads limit.
	Initially, we start with S threads for the solver (defined in setup.solvers) + 1 fixed thread for the validator. When a solver
	finishes execution, it says so and validator can potentially grab a thread, if needed (when there are some candidates pending).
	If it does so, it generally just runs one workload ("helper_workload" candidates) and releases the thread. In that case it
	can be re-used by both solvers and validators. But if the number of pending candidates exceeds "high_watermark", it actually
	keeps the thread until the number of candidates goes below "low_watermark". The idea here is that too many candidates shifts
	resources to the validator, since the backlog becomes too large. In that case solver might slow down, since it cannot
	reserve threads anymore, which is what we try to achieve (like, "stop I need to catch on"). Note, all this applies only
	if the searchlight.dynamic_scheduling is set to 1. Otherwise, the validator uses no more than max_helpers additional threads.
* "rebal_watermark". This is for deciding when to send candidates to a _remote_ validator. If >= tha this number, we'll first
	try that and only then to span more threads.

### The "searchligh.sampler" section
This one controls the sampler, which handles all the synopses for the query. The params:

* searchligh.sampler.cache. Defines the cache policy for the synopsis array. "eager" means it allocates memory for all
	synopsis cells right away (by cell I mean the unit that contains aggregate info about a sub-array). "lazy" means it takes
	SciDb array chunks configuration into accout. The idea here is to reduce the memory footprint. If a synopsis cell isn't
	required, we won't access that data. This is useful, since Solvers don't generally touch the whole array (the search is divided),
	but it creates additional runtime overhead. "none" means we don't cache at all. Every time we want synopsis information
	we go to disk (sans whatever caching SciDb does for its arrays).
* searchligh.sampler.preload The synopsis data itself isn't loaded immediately at the query start unless this param is set to 1.
	The synopsis array is read lazily, when solvers really need the data. Setting it to 1 in general isn't a good idea, since
	the query warm up time might go up significantly (you need to read all synopses into memory, which might be a lot of data).
	However, the good news is that once loaded the synopses persist _between_ the queries. So, you might try that one and
	don't restart SciDb for subsequent queries. Still might be a fair comparison, if we can argue that synopses should stay in
	memory for a series of queries.
* searchligh.sampler.memory_per_attr. So this is basically to restrict the number of synopsis we cache in memory. When
	we run a query, we can specify multiple synopses per attribute (e.g., different resolutions, from coarser to finer). So,
	we'll cache as much as we can, but the rest will be read from disk. This is in MB.
* searchlight.sampler.cell_thr/mbr_thr. This tuned the estimations via synopses. I think it was described in the paper/talks, but
	not sure. Say, you have a bunch of synopses of different resolutions: from coarse to fine. A query region doesn't have to be aligned with the
	synopsis grid. So, first you actually determine the Minimum-Bounding-Rectangle (MBR) aligned with the synopsis grid that
	contains the query region. If the ratio of query region to MBR area is above the mbr_thr, then only the current synopsis will be used for
	estimations. If below, then we divide the query region in chunks and the current synopsis is used to compute only the chunks
	that overlap with the corresponsing synopsis cell on more than cell_thr (e.g., 0.7 means the synopsis cell must be covered
	70% by the query region). Otherwsie, this query region part will be computed with a finer synopsis. Needless to say that works only
	when: (a) we have multiple synopses; (b) the aggregates are combinable.


## The "balance" section
"balance" section is for redistributing the search space between solvers (including remote ones) and solutions between validators.

* "solver_balance" enables the _dynamic_ solver balancing. That means the solver can offload some of its work to a local or even
	remote solver. Note, only "sl" and "split" search heuristics use that. The others just ignore the param.
* "general_low/high" related to the dynamic balancing. The idea here is that we monitor the search and look at the current search space
	size (product of current variable domains). If the ratio of the current size to the initial search space size is
	between "general_low" and "general_high", then we just cut this subtree from the search and transfer it to another
	solver (local or remote). The current solver just rolls back and goes into another search path.
* "validator_balance" decides how we submit candidates to other validators. When a validator gets a candidate it can either
	validate it locally or submit to another validsator. "stripes" means each validator takes a "horizontal" stripe of the array
	(as described in the paper) and we submit candidates based on that. "dynamic" means we actually take the dynamic array distribution
	into account. The idea here is to send the candidate (which is a query region) to the validator that has most of the data needed
	for the validation locally. Another param, "map_update_frequency", sets the periodicity with which Searclight broadcasts the
	distribution of data from validators.
* "slices_number". Total number of search slices. The search array is divided into slices along the largest dimension. So,
	each solver takes at least one slice. If the number of slices > number of solvers, the slices are distributed round-robin.
* "load_slices". 1 means it will create a single search with all local slices at once. While 0 means it will create
	several seacrhes, one per slice. There were some performance implications. Probably should be 1 always.

## The "setup" section
"setup" basically tells the Searchlight where to put Solvers in Validators.

Each SciDb instance has instance_id. I think they are assigned linearly in the config order. The SciDb config I use usually looks
like this:

server-0=host-1,0
server-1=host-2,1
...
server-N=host-N,1

That puts one SciDb instance on each machine (for server-0 it says 0, but it always puts an instance there, since server-0 is the coordinator).
In that case server-I basically will have instance number I. In that case "searchlight.setup.solvers" basically says how many
solver threads should be put on each instance. It should usually equal to the number of virtual CPU cores. "searchlight.setup.validators"
basically an array specifying which instances will run the validator. Usually should include all instances. But obviously, you can
play with these param.

Note, when you specify S solvers for instance I, that means Searchlight will use S+1 threads when dynamically balancing work
between solvers and validator on that instance. This is decribed in more detail in the validator section.

## Relaxation/contraction
This is for enabling relaxation and contraction mechanism. This is described well in the paper.

* relax.on Set to 1 basically allows it. Note, relaxation or contraction is determined automatically, based on the params you
	specify in the tasks.
* relax.card Desired result cardinality to determine if we want to run it.
* relax.dist_w Defined the weight of the distance part of the relaxation value (remember, it takes into account the distance and
	the number of the failed constraints). The default is 1/2.
* relax.spec The number of speculative solvers to run. The paper describes it. Basically, when we get fails, we run relaxation
	on the corresponding sub-trees immediately instead of waiting for the main search to finish. Eats away CPU obviously. When
	the main search is finished, the speulative solvers turn off and only the main solver threads are used for search.
* relax.heur. This defines the fail registration heuristic. When the search fails, we need to understand if we want to replay
	it later. For this we need the current values for expressions participating in constraints. We can do it 3 ways:
	* "all". That means we estimate every expression via synopses. Might be expensive, since if the search fails, not all
		constraints are necessarily computed (e.g., if one constraint out of five fails, the search fails).
	* "guess". This one tries to remedy the "all" drawback by taking the current values of the expressions. Since this is quite
		approximate, it might result in more fail replays than needed.
	* "guess-all". Combines the two. Initially, we use "guess". But if we have a fail _during_ replay, we try "all" to refine
		our estimations.
	This param is important, oftem I ran experiments with different settings to see which one performs better.
* relax.replay. Determines which constraints to relax when replaying a fail. "viol" means relaxing only violated consraints, which
	is the bare minimum. That might result in fails during replay (re-replay), since other constraints might fail. "all" relaxes all constraints.
	While it avoids re-replays, it might be more expensive, since relaxing too much increases the search space and potential number
	of candidates. I think "viol" is usually better.
* relax.sort Defines the order of replays. Each replay has the best and worst relaxation possible. Worst is basically the maximum one.
	If you relax that much, it will succeed. Best is the minimal relaxation you can try on replay to move forward, but you might get
	another fail quite fast, if the constraint fails again. "best" sorts the replays based on the best relaxation. "worst" based on the
	worst. "prod" will sort based on the product of the two. "time" will sort based on the time the replays were detected (this is useful
	to prove that in-time relaxation is usually inferior to the best/worst). I think "best" is usually better, since such candidates tend
	to be the "best bang for the buck". Cheaper to relax and often result in answers.
* relax.replay_rd As stated in the previous paragraph. We have a range of possible relaxations: from best (smallest) to worst (largest).
	This parameter is from 0.0 to 1.0, and defines the degree of relaxation, which will be "best + (worst - best) x replay_rd". For instance,
	1.0 means we use the worst, which is quite loose, but avoids re-replays. This is quite useful for experimentation. Usually, you want
	to play with it.
* relax.lrd In the paper we introduce the concept of relaxation degree as [0,1]-normalized degree of a constraint relaxation, 0 being the
	best possible (i.e., equivalent to the original constraint). LRD is basically the Lowest Relaxation Degree currently possible. If for the
	current fail, your best relaxation degree is greater than LRD, you can safely get rid of it, since it will never result in meaningful results.
	It's synchronized over the nodes periodically. So this param sets the initial LRD. Should be always 1.0 for correctness. But it's useful to
	see "what if" it was lower in the beginning. The idea is that we want to try to make it as lower as possible, since this is the only
	way to provably get rid of fail replays.
* relax.force_replays This is in case we relax and get immediate fail after that. That might happen if we use "guess" heuristic and
	use "viol" for relax.replay. In that case if this param is 1, we replay it immediately. Should be always 1, I think.
* relax.save_udfs So this is the optimization described in the paper to save functions state when saving a fail. Has a small footprint,
	actually. Should be always 1.
* relax.contr_type. This is to define the _contraction_ mode. 0 means no contraction. 1 means rank-based (ranks are described in the paper).
	2 means the skyline contraction (again, based on ranks).


## Task-dependent sections
The next section's name is dependent on the task. Each task (i.e., the query function in the SO lib) looks for params in the specific
section. You can look up the name in the sources. MIMIC tasks look for "mimic" JSON map, Semantic Windows (SW) look for "sw", etc. I'll call
it the task section. You can actually have a single JSON for multiple tasks, since SW, for instance, won't read the "mimic" session, and so on.
But I usually keep them separate. There's a bunch of common params:
* task_name.db (e.g., sw.db) defines the search heuristic, which I described earlier. For most runs I use "split", since it works very well.
* task_name.time_limit (e.g., sw.time_limit). This is the total time limit for the search. After that it'll be terminated. The search still outputs
	the results it has found so far, since Searchlight is online search.

The rest of params depends on the semantics of the task, but some tasks have common ones.

### Semantic Windows
This task is defined in the "sw" section. The task is basically described in the Semantic Windows paper. In a nutshell, we
define a region (window), possibly specifying ranges for the sides, and the neighborhood (basically, hypercube-like expansion) around it.
Then, we set criteria for the region and the neighborhood. 

The params are:
* lx, ux, ly, uy. Define the total search array in the array. Allows to search in a sub-array, instead of the entire array. The task is 2D only.
* len_lx, len_ux, len_ly, len_uy. Define the ranges of possible lengths for the window (l and u mean lower and upper).
	If len_lx=len_ux, for instance, that means the window has fixed sized over the x dimension.
* size_l, size_u define the range of possible areas for the window (as in length x width).
* step_x, step_y define the step of the window. Semantically, the left-upper corner starts at lx,ly and then moves with this step. Allows
	to reduce the search space size, since otherwise it might be way too large.
* step_lx, step_ly is related to len_l params. Instead of trying all side lengths in the interval, we go with with step. Example:
	len_lx=10, len_ux=50, step_lx=20 will try lengths 10, 30, 40.
* avg_l, avg_u define the range of avg(val) function we're looking for. This is the average constraint. It's always computed over the array
	attribute called "val".
* avg_relax_l, avg_relax_u. This is the maximum relaxation for the avg() constraint. To avoid too huge relaxation space. Should be always defined
	to some meanigful degree. Defaults might be too much loose.
* avg_contr Sets the contraction type for the avg() constraint. If set to 1 we refer larger values. -1 we prefer smaller values. 0 means
	we don't contract it. Has obvious effect on the contraction ranks. Can be set to 0 or 1 almost always.
* neighborhood.size Defines the size of the neigborhood of the region in all directions. So, the extension is always the same in all directions.
* neighborhood.max_diff, neighborhood.min_diff. If not 0, basically create the following constraints: "rmax - nmax >= max_diff",
	"rmin - nmax >= min_diff", where rmax/min is the max(val)/min(val) in the query region, and nmax is the max(val) in the neighborhood (not
	including the region, just the neighborhood itself).
* neighborhood.max_relax_diff, max_diff_high. So this sets the relaxation interval for the "rmax - nmax >= max_diff" constraint. The constraint
	is satisfied with the [max_diff; +inf] interval. max_relax_diff is basically to specify relaxation (must be lower than max_diff), but max_diff_high
	basically sets the maximum value for the constraint function. It's needed to basically compute relaxation ranks properly. Should be set to
	a suitable domain value, not too high.
* neighborhood.min_relax_diff, min_diff_high. The same but for the "rmin - nmax >= min_diff" constraint.
* neighborhood.max_diff_contr, min_diff_contr. Analogous to avg_contr for the constraints above.

### MIMIC
The MIMIC task (mimic section) is kind of similar to Semantic Windows, bu runs over the MIMIC dataset. The MIMIC is also a 2D array, where X dimesion is patient IDs,
Y dimension is normalized time. Every (id, time) point contains a measurement. I usually run MIMIC on the ABP attribute, which is Arterial Blood
Pressure. The task can be run over a single attribute. The window here is basically 1D, over time only. But the search goes over multiple
patients, so the window will be a rectangle of size 1 x len.

The params:
* l_id, u_id, l_time, u_time the same aas lx,ly, ux, uy for the SW task.
* avg_l, avg_u are the same.
* avg_relax_l, avg_relax_u, avg_contr is the same as for SW.
* len_l, len_u is analogous to len_ly, len_uy.
* step_len is analogous to step_ly.
* step_time is analogous to step_y.
* signal defines the measurement to look at (array attribute).
* neighborhood.lsize, neighborhood.rsize define the neighborhood of the time interval left and right respectively. Since the window is 1D, no
	need to specify it for patient ids, only for time.
* neighborhood.left_max_diff, neighborhood.right_max_diff are analogous to max_diff for the sw task. It adds the following constraints:
	rmax - lmax >= left_max_diff and rmax - rmax >= right_max_diff, where rmax is region max, lmax/rmax are the max of the left/right neighborhoods.
* neighborhood.left_relax_diff, left_diff_high. Analogous to the SW; defines relaxation interval for the left_max_diff constraint.
* neighborhood.right_relax_diff, right_diff_high. Analogous to the above, defines relaxation interval for the right_max_diff constraint.

### SDSS
This is the SDSS task (sdss section), which we ran for the first paper. Basically just Semantic Windows, bit over sdss data. No neighborhood.
The parameters are the same as for sw with the addition of: min_l, min_h, which are analogous to avg() but for min(). The avg_l/u, min_l/u
params here are actually arrays, since sdss always assumes 5 attributes: u, g, r, i, z. So all contraints are set for each attribute and the
param arrays are treated in this attribute order.
 
### DIST
This is applicable to the MIMIC data set. The task has section "dist" in the config. The idea here is that you're given a sequence of measurements
of some lenght (e.g., 10 20 10 30), which might represent a sequence of ABP, for instance. And in your data you want to find all sub-sequences of the
same length, which are similar to the given one. The similarity is defined as the Euclidean distance between the two. It needs to be below the specified
threshold. There're some techniques to make the computation faster, like PAA or FFT synopses (which we support), but this is an advanced topic. Not sure, if we want to
run this task at all right now. But we used for SIGMOD demo. The params:

* l_id, u_id, l_time, u_time, step_time, signal. The same as for MIMIC.
* dist Similarity threshold for the Eucledean distance.
* dist_relax relaxed threshold for relaxation purposes.
* query. Filename where the query sequence is stored. The format is: length m1 m2 ... For example, 10 10 20 will be stored as "3 10 10 20" in the file.
