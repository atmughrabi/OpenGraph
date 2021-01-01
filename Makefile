#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin/open-graph name doesn't need to match main/open-graph.c
export APP                ?= open-graph

# test name needs to match the file name test/test_accel-graph.c
export APP_TEST           ?= test_open-graph-match

# export APP_TEST           ?=  sweep_order-OpenGraph-performance-graph
# export APP_TEST           ?=  sweep_order-PR-performance-graph
# export APP_TEST           ?=  sweep_order-BFS-performance-graph

# dirs Root app
export APP_DIR              ?= .

#dir root/managed_folders
export SRC_DIR           	= src
export OBJ_DIR			  	= obj
export INC_DIR			  	= include
export BIN_DIR			  	= bin
export RES_DIR			  	= results


#if you want to compile from cmake you need this directory
#cd build
#cmake ..
export BUILD_DIR		  	= build

# relative directories used for managing src/obj files
export STRUCT_DIR		  	= structures
export PREPRO_DIR		  	= preprocess
export ALGO_DIR		  		= algorithms
export UTIL_DIR		  		= utils
export CAPI_UTIL_DIR		= capi_utils

#contains the tests use make run-test to compile what in this directory
export TEST_DIR		  	= tests

#contains the main for the graph processing framework
export MAIN_DIR		  	= main

##################################################
##################################################

#########################################################
#       		 ACCEL RUN GRAPH ARGUMENTS    			#
#########################################################

# export BENCHMARKS_DIR    	?= ../../01_GraphDatasets
export BENCHMARKS_DIR    	?= ../01_test_graphs

export GRAPH_SUIT ?= TEST
# export GRAPH_SUIT ?= LAW
# export GRAPH_SUIT ?= GAP
# export GRAPH_SUIT ?= SNAP
# export GRAPH_SUIT ?= KONECT
# export GRAPH_SUIT ?= GONG

# TEST # small test graphs
# export GRAPH_NAME ?= test
# export GRAPH_NAME ?= v51_e1021
# export GRAPH_NAME ?= v300_e2730
export GRAPH_NAME ?= graphbrew

# GONG # https://gonglab.pratt.duke.edu/google-dataset
# export GRAPH_NAME ?= GONG-gplus
# export GRAPH_NAME ?= Gong-gplus

# GAP # https://sparse.tamu.edu/MM/GAP/
# export GRAPH_NAME ?= GAP-twitter
# export GRAPH_NAME ?= GAP-road

# SNAP # https://snap.stanford.edu/data/
# export GRAPH_NAME ?= SNAP-cit-Patents
# export GRAPH_NAME ?= SNAP-com-Orkut
# export GRAPH_NAME ?= SNAP-soc-LiveJournal1
# export GRAPH_NAME ?= SNAP-soc-Pokec
# export GRAPH_NAME ?= SNAP-web-Google

# KONECT # http://konect.cc/networks/wikipedia_link_en/
# export GRAPH_NAME ?= KONECT-wikipedia_link_en

# LAW # https://sparse.tamu.edu/MM/LAW/
# export GRAPH_NAME ?= LAW-amazon-2008
# export GRAPH_NAME ?= LAW-arabic-2005
# export GRAPH_NAME ?= LAW-cnr-2000
# export GRAPH_NAME ?= LAW-dblp-2010
# export GRAPH_NAME ?= LAW-enron
# export GRAPH_NAME ?= LAW-eu-2005
# export GRAPH_NAME ?= LAW-hollywood-2009
# export GRAPH_NAME ?= LAW-in-2004
# export GRAPH_NAME ?= LAW-indochina-2004
# export GRAPH_NAME ?= LAW-it-2004
# export GRAPH_NAME ?= LAW-ljournal-2008
# export GRAPH_NAME ?= LAW-uk-2002
# export GRAPH_NAME ?= LAW-uk-2005
# export GRAPH_NAME ?= LAW-webbase-2001

# export FILE_BIN_TYPE ?= graph
export FILE_BIN_TYPE ?= graph.bin
# export FILE_BIN_TYPE ?= graph.wbin

# export FILE_LABEL_TYPE ?= graph_Gorder.labels
export FILE_LABEL_TYPE ?= graph_Rabbit.labels

#GRAPH file
export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_SUIT)/$(GRAPH_NAME)/$(FILE_BIN_TYPE)
export FILE_LABEL = $(BENCHMARKS_DIR)/$(GRAPH_SUIT)/$(GRAPH_NAME)/$(FILE_LABEL_TYPE)

#ALGORITHM
export PULL_PUSH 		?= 0
export ALGORITHMS 		?= 1

#GRAPH DATA_STRUCTURES
export SORT_TYPE		?= 1
export DATA_STRUCTURES  ?= 0
export REORDER_LAYER1 	?= 4
export REORDER_LAYER2   ?= 0
export REORDER_LAYER3   ?= 0

#ALGORITHM SPECIFIC ARGS
export ROOT 			?= 46050
export TOLERANCE 		?= 1e-8
export DELTA			?= 800
export NUM_ITERATIONS	?= 100

#PERFORMANCE
# export NUM_THREADS_PRE  ?= $(shell grep -c ^processor /proc/cpuinfo)
# export NUM_THREADS_ALGO ?= $(shell grep -c ^processor /proc/cpuinfo)
# export NUM_THREADS_KER  ?= 1

export NUM_THREADS_PRE  ?= 1
export NUM_THREADS_ALGO ?= 1
export NUM_THREADS_KER  ?= 1

#EXPERIMENTS
export NUM_TRIALS 		?= 1

#GRAPH FROMAT EDGELIST
export FILE_FORMAT		?= 1
export CONVERT_FORMAT 	?= 0

#STATS COLLECTION VARIABLES
export BIN_SIZE 		?= 1000
export INOUT_STATS 		?= 0
export MASK_MODE 		?= 0

##################################################

APP_DIR                 = .
MAKE_DIR                = 00_graph_bench

MAKE_NUM_THREADS        = $(shell grep -c ^processor /proc/cpuinfo)
MAKE_ARGS               = -w -C $(APP_DIR)/$(MAKE_DIR) -j$(MAKE_NUM_THREADS)

#########################################################
#                RUN  ARGUMENTS                         #
#########################################################
export ARGS ?= -k -M $(MASK_MODE) -j $(INOUT_STATS) -g $(BIN_SIZE) -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS_PRE) -N $(NUM_THREADS_ALGO) -K $(NUM_THREADS_KER) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -F $(FILE_LABEL) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -b $(DELTA)

##################################################
##################################################

##############################################
#         ACCEL GRAPH TOP LEVEL RULES        #
##############################################

.PHONY: help
help:
	$(MAKE) help $(MAKE_ARGS)

.PHONY: run
run:
	$(MAKE) run $(MAKE_ARGS)

.PHONY: sweep-run
sweep-run:
	$(MAKE) run-test $(MAKE_ARGS)

.PHONY: run-openmp
run-openmp:
	$(MAKE) run-openmp $(MAKE_ARGS)

.PHONY: convert
convert:
	$(MAKE) convert $(MAKE_ARGS)

.PHONY: sweep-convert
sweep-convert:
	$(MAKE) sweep-convert $(MAKE_ARGS)

.PHONY: echo-dir
echo-dir: 
	$(MAKE) echo-dir $(MAKE_ARGS)

.PHONY: convert-w
convert-w:
	$(MAKE) convert-w $(MAKE_ARGS)

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

.PHONY: clean-all
clean-all: clean

.PHONY: scrub
scrub: clean clean-nohup clean-stats

.PHONY: clean-stats
clean-stats:
	$(MAKE) clean-stats $(MAKE_ARGS)

.PHONY: clean-nohup
clean-nohup:
	@rm -f $(APP_DIR)/nohup.out

##################################################
##################################################