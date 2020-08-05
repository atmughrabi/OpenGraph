#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin/accel-graph name doesn't need to match main/accel-graph.c
APP                        = open-graph

# test name needs to match the file name test/test_accel-graph.c
export APP_TEST          ?=  test_open-graph
# export APP_TEST          ?=  pagerRank-accuracy-report
# export APP_TEST          ?=  pagerRank-capi-report
# export APP_TEST          ?=  test_cache

# dirs Root app
export APP_DIR              ?= .


export BENCHMARKS_DIR_LOCAL ?= 01_test_graphs

export BENCHMARKS_DIR    	?= ../$(BENCHMARKS_DIR_LOCAL)
# export BENCHMARKS_DIR    	= ../../$(BENCHMARKS_DIR_LOCAL)

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

# # small test graphs
# export GRAPH_NAME ?= test
export GRAPH_NAME ?= v51_e1021
# export GRAPH_NAME ?= v300_e2730
# export GRAPH_NAME ?= amazon


# GAP https://sparse.tamu.edu/MM/GAP/
# https://gonglab.pratt.duke.edu/google-dataset

# export GRAPH_NAME ?= Gong-gplus
# export GRAPH_NAME ?= GAP-road
# export GRAPH_NAME ?= SNAP-soc-pokec
# export GRAPH_NAME ?= SNAP-cit-Patents
# export GRAPH_NAME ?= SNAP-com-orkut
# export GRAPH_NAME ?= SNAP-soc-LiveJournal1
# export GRAPH_NAME ?= KONECT-wikipedia_link_en

# LAW https://sparse.tamu.edu/MM/LAW/
# export GRAPH_NAME ?= amazon-2008
# export GRAPH_NAME ?= arabic-2005
# export GRAPH_NAME ?= cnr-2000
# export GRAPH_NAME ?= dblp-2010
# export GRAPH_NAME ?= enron
# export GRAPH_NAME ?= eu-2005
# export GRAPH_NAME ?= hollywood-2009
# export GRAPH_NAME ?= in-2004
# export GRAPH_NAME ?= indochina-2004
# export GRAPH_NAME ?= it-2004
# export GRAPH_NAME ?= ljournal-2008
# export GRAPH_NAME ?= sk-2005
# export GRAPH_NAME ?= uk-2002
# export GRAPH_NAME ?= uk-2005
# export GRAPH_NAME ?= webbase-2001

# export LAW = amazon-2008 arabic-2005 cnr-2000 dblp-2010 enron eu-2005 hollywood-2009 in-2004 indochina-2004 it-2004 ljournal-2008 sk-2005 uk-2002 uk-2005 webbase-2001
# export MIX = Gong-gplus GAP-road SNAP-soc-pokec SNAP-cit-Patents SNAP-com-orkut SNAP-soc-LiveJournal1 KONECT-wikipedia_link_en
export LAW ?= amazon-2008 arabic-2005 cnr-2000 dblp-2010 enron eu-2005 hollywood-2009 in-2004 indochina-2004 it-2004 ljournal-2008 sk-2005 uk-2002 uk-2005 webbase-2001 Gong-gplus GAP-road SNAP-soc-pokec SNAP-cit-Patents SNAP-com-orkut SNAP-soc-LiveJournal1 KONECT-wikipedia_link_en gplus USA-Road enwiki-2013 KONECT-wikipedia_link_en twitter

# export GAP = GAP-kron GAP-road GAP-twitter GAP-urand GAP-web
# export CU_CONFIG_MODES ?= 0x00000000 0x00041000 0x00841000 0x10041000 0x10841000
# export CU_CONFIG_MODES ?= 0x10000000 0x00800000 0x00040000 0x00001000
export CU_CONFIG_MODES  ?= 0x00041000 0x00841000
# export PUSHPULL_MODES = 0 2 4 9 10 11 12 13

# TEXT formant
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph

#UNWEIGHTED
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.bin

# export FILE_BIN_TYPE = graph
# export FILE_BIN_TYPE = graph.bin
export FILE_BIN_TYPE ?= graph.wbin

#WEIGHTED
export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/$(FILE_BIN_TYPE)

#Direction
export PULL_PUSH 		?=0

#GRAPH RUN
export SORT_TYPE		?=0
export REORDER 		    ?=0
export DATA_STRUCTURES  ?=0
export ALGORITHMS 		?=2

export ROOT 			?=164
export TOLERANCE 		?=1e-8
export DELTA			?= 800

export START_THREADS	?= 1
export INC_THREADS      ?=1
export NUM_THREADS  	?=25
# NUM_THREADS  	= $(shell grep -c ^processor /proc/cpuinfo)
export NUM_ITERATIONS	?= 1
export NUM_TRIALS 		?=1

export FILE_FORMAT		?= 1
export CONVERT_FORMAT 	?=1

#STATS COLLECTION VARIABLES
export BIN_SIZE 		?=1000
export INOUT_STATS 		?=0

##################################################

APP_DIR                 = .
MAKE_DIR                = 00_graph_bench

MAKE_NUM_THREADS        = $(shell grep -c ^processor /proc/cpuinfo)
MAKE_ARGS               = -w -C $(APP_DIR)/$(MAKE_DIR) -j$(MAKE_NUM_THREADS)

#########################################################
#                RUN  ARGUMENTS                         #
#########################################################
export ARGS ?= -j $(INOUT_STATS) -g $(BIN_SIZE) -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA)
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

.PHONY: run-cache
run-cache:
	$(MAKE) run-cache $(MAKE_ARGS)

.PHONY: debug-cache
debug-cache:
	$(MAKE) debug-cache $(MAKE_ARGS)

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

.PHONY: law
law:
	$(MAKE) law $(MAKE_ARGS)

.PHONY: mix
mix:
	$(MAKE) mix $(MAKE_ARGS)

.PHONY: results
results:
	$(MAKE) results $(MAKE_ARGS)

.PHONY: results-law
results-law:
	$(MAKE) results-law $(MAKE_ARGS)

.PHONY: results-mix
results-mix:
	$(MAKE) results-mix $(MAKE_ARGS)

PHONY: stats
stats:
	$(MAKE) stats $(MAKE_ARGS)

.PHONY: stats-law
stats-law:
	$(MAKE) stats-law $(MAKE_ARGS)

.PHONY: stats-mix
stats-mix:
	$(MAKE) stats-mix $(MAKE_ARGS)

##################################################
##################################################