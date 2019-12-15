


#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin/open-graph name doesn't need to match main/open-graph.c
export GAPP               = open-graph

# test name needs to match the file name test/test_accel-graph.c
export GAPP_TEST          = test_open-graph

# dirs Root app 
export APP_DIR           	= .
export BENCHMARKS_DIR    	= ../03_test_graphs

#dir root/managed_folders
export SRC_DIR           	= src
export OBJ_DIR			  	= obj
export INC_DIR			  	= include
export BIN_DIR			  	= bin

#if you want to compile from cmake you need this directory
export BUILD_DIR		  	= build

# relative directories used for managing src/obj files
export STRUCT_DIR		  	= structures
export PREPRO_DIR		  	= preprocessing
export ALGO_DIR		  		= graphalgorithms
export UTIL_DIR		  		= utils

#contains the tests use make run-test to compile what in this directory
export TEST_DIR		  	= tests

#contains the main for the graph processing framework
export MAIN_DIR		  	= main

##################################################
##################################################

#########################################################
#       		 ACCEL RUN GRAPH ARGUMENTS    			#
#########################################################

#small test graphs
export GRAPH_NAME = test
export GRAPH_NAME = v51_e1021
# export GRAPH_NAME = v300_e2730

#gem5-Aladdin small dynamic traces
# export GRAPH_NAME = Gnutella
# export GRAPH_NAME = dblp
# export GRAPH_NAME = amazon
# export GRAPH_NAME = euall

# generates large dynamic traces for gem5-Aladdin
# export GRAPH_NAME = com-youtube
# export GRAPH_NAME = web-BerkStan
# export GRAPH_NAME = web-Google
# export GRAPH_NAME = wiki-Talk

# synthetic graphs
# export GRAPH_NAME = RMAT20
# export GRAPH_NAME = RMAT22

# real world large graphs binary format
# export GRAPH_NAME = orkut
# export GRAPH_NAME = gplus
# export GRAPH_NAME = sk-2005
# export GRAPH_NAME = twitter
# export GRAPH_NAME = livejournal
# export GRAPH_NAME = USA-Road
# export GRAPH_NAME = enwiki-2013
# export GRAPH_NAME = arabic-2005


#UNWEIGHTED
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.bin

#WEIGHTED
export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.wbin



#GRAPH RUN
export SORT_TYPE 		= 0
export REORDER 		    = 0
export DATA_STRUCTURES  = 0
export ALGORITHMS 		= 1

export ROOT 			= 164
export PULL_PUSH 		= 2
export TOLERANCE 		= 1e-8
export DELTA 			= 800

# export NUM_THREADS  	= 64
NUM_THREADS  	= $(shell grep -c ^processor /proc/cpuinfo)
export NUM_ITERATIONS 	= 1
export NUM_TRIALS 		= 1

export FILE_FORMAT 	= 1
export CONVERT_FORMAT 	= 1

#STATS COLLECTION VARIABLES
export BIN_SIZE = 512
export INOUT_STATS = 2

export ARGS = -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA)
##################################################

APP_DIR           	= .
MAKE_DIR      = 00_graph_bench
MAKE_NUM_THREADS  	= $(shell grep -c ^processor /proc/cpuinfo)
MAKE_ARGS = -w -C $(APP_DIR)/$(MAKE_DIR) -j$(MAKE_NUM_THREADS)

##################################################
##################################################

##############################################
#         OPEN GRAPH TOP LEVEL RULES        #
##############################################

.PHONY: help
help:
	$(MAKE) help $(MAKE_ARGS)

.PHONY: run
run:
	$(MAKE) run $(MAKE_ARGS)

.PHONY: run-openmp
run-openmp:
	$(MAKE) run-openmp $(MAKE_ARGS)

.PHONY: convert
convert:
	$(MAKE) convert $(MAKE_ARGS)

.PHONY: stats-openmp
stats-openmp: graph-openmp
	$(MAKE) stats-openmp $(MAKE_ARGS)

.PHONY: debug-openmp
debug-openmp: 
	$(MAKE) debug-openmp $(MAKE_ARGS)

.PHONY: debug-memory-openmp
debug-memory-openmp: 
	$(MAKE) debug-memory-openmp $(MAKE_ARGS)

.PHONY: test-verbose
test-verbose:
	$(MAKE) test-verbose $(MAKE_ARGS)
	
# test files
.PHONY: test
test:
	$(MAKE) test $(MAKE_ARGS)
	
.PHONY: run-test
run-test: 
	$(MAKE) run-test $(MAKE_ARGS)

.PHONY: run-test-openmp
run-test-openmp:
	$(MAKE) run-test-openmp $(MAKE_ARGS)

.PHONY: debug-test-openmp
debug-test-openmp: 
	$(MAKE) debug-test-openmp $(MAKE_ARGS)

.PHONY: debug-memory-test-openmp
debug-memory-test-openmp:	
	$(MAKE) debug-memory-test-openmp $(MAKE_ARGS)
# cache performance
.PHONY: cachegrind-perf-openmp
cachegrind-perf-openmp:
	$(MAKE) cachegrind-perf-openmp $(MAKE_ARGS)

.PHONY: cache-perf
cache-perf-openmp: 
	$(MAKE) cache-perf-openmp $(MAKE_ARGS)

.PHONY: clean
clean: 
	$(MAKE) clean $(MAKE_ARGS)

.PHONY: clean-obj
clean-obj: 
	$(MAKE) clean-obj $(MAKE_ARGS)

##################################################
##################################################

