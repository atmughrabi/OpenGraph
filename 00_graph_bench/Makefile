
#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin/accel-graph name doesn't need to match main/accel-graph.c
APP               ?= open-graph

# test name needs to match the file name test/test_accel-graph.c
APP_TEST          ?= test_open-graph


# dirs Root app
APP_DIR           	?= .
# BENCHMARKS_DIR    	?= ../../01_GraphDatasets
BENCHMARKS_DIR    	?= ../01_test_graphs

#dir root/managed_folders
SRC_DIR           	?= src
OBJ_DIR			  	?= obj
INC_DIR			  	?= include
BIN_DIR			  	?= bin
RES_DIR			  	?= results

#if you want to compile from cmake you need this directory
#cd build
#cmake ..
BUILD_DIR		  	?= build

# relative directories used for managing src/obj files
STRUCT_DIR		  	?= structures
PREPRO_DIR		  	?= preprocess
ALGO_DIR		  	?= algorithms
UTIL_DIR		  	?= utils
CAPI_UTIL_DIR		?= capi_utils

# Folders needed when using CAPI


#contains the tests use make run-test to compile what in this directory
TEST_DIR		  	?= tests

#contains the main for the graph processing framework
MAIN_DIR		  	?= main

##################################################
##################################################

#########################################################
#       		 ACCEL RUN GRAPH ARGUMENTS    			#
#########################################################

#small test graphs
GRAPH_NAME ?= test
# GRAPH_NAME ?= v51_e1021
# GRAPH_NAME  ?= v300_e2730

#gem5-Aladdin small dynamic traces
# GRAPH_NAME ?= Gnutella
# GRAPH_NAME ?= dblp
# GRAPH_NAME ?= amazon
# GRAPH_NAME ?= euall

# generates large dynamic traces for gem5-Aladdin
# GRAPH_NAME ?= com-youtube
# GRAPH_NAME ?= web-BerkStan
# GRAPH_NAME ?= web-Google
# GRAPH_NAME ?= wiki-Talk

# synthetic graphs
# GRAPH_NAME ?= RMAT20
# GRAPH_NAME ?= RMAT22

# real world large graphs binary format
# GRAPH_NAME ?= orkut
# GRAPH_NAME ?= gplus
# GRAPH_NAME ?= sk-2005
# GRAPH_NAME ?= twitter
# GRAPH_NAME ?= livejournal
# GRAPH_NAME ?= USA-Road
# GRAPH_NAME ?= enwiki-2013
# GRAPH_NAME ?= arabic-2005

# real world large graphs binary LAW
# GRAPH_NAME ?= amazon-2008
# GRAPH_NAME ?= arabic-2005
# GRAPH_NAME ?= cnr-2000
# GRAPH_NAME ?= dblp-2010
# GRAPH_NAME ?= enron
# GRAPH_NAME ?= eu-2005
# GRAPH_NAME ?= hollywood-2009
# GRAPH_NAME ?= in-2004
# GRAPH_NAME ?= indochina-2004
# GRAPH_NAME ?= it-2004
# GRAPH_NAME ?= ljournal-2008
# GRAPH_NAME ?= sk-2005
# GRAPH_NAME ?= uk-2002
# GRAPH_NAME ?= uk-2005
# GRAPH_NAME ?= webbase-2001

LAW ?= amazon-2008 arabic-2005 cnr-2000 dblp-2010 enron eu-2005 hollywood-2009 in-2004 indochina-2004 it-2004 ljournal-2008 sk-2005 uk-2002 uk-2005 webbase-2001
MIX ?= GAP-kron GAP-road GAP-twitter GAP-urand GAP-web
CU_CONFIG_MODES ?= 0x00000000 0x00041000 0x00841000 0x10041000 0x10841000

# TEXT formant
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph

#UNWEIGHTED
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.bin

# FILE_BIN_TYPE ?= graph
# FILE_BIN_TYPE ?= graph.bin
FILE_BIN_TYPE ?= graph.wbin

#WEIGHTED
FILE_BIN ?= $(BENCHMARKS_DIR)/$(GRAPH_NAME)/$(FILE_BIN_TYPE)



#GRAPH RUN
SORT_TYPE 		?= 0
REORDER 		?= 0
DATA_STRUCTURES ?= 0
ALGORITHMS 		?= 1

ROOT 			?= 164
PULL_PUSH 		?= 2
TOLERANCE 		?= 1e-8
DELTA 			?= 800

START_THREADS   ?= 1
INC_THREADS     ?= 1
NUM_THREADS  	?= 64
# NUM_THREADS  	= $(shell grep -c ^processor /proc/cpuinfo)
NUM_ITERATIONS 	?= 1
NUM_TRIALS 		?= 1

FILE_FORMAT 	?= 1
CONVERT_FORMAT 	?= 1

#STATS COLLECTION VARIABLES
BIN_SIZE ?= 512
INOUT_STATS ?= 0

ARGS ?= -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA)


##############################################
#      ACCEL GRAPH COMPILATION VARIABLES     #
##############################################


#each implementation is contained in a separate folder. when you compile for openmp for example the files in openmp folder gets included
SRC_FILES_ALGO_OPENMP		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/*.c)
SRC_FILES_ALGO_CACHE		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/cachesim/*.c)

SRC_FILES_UTIL		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/*.c)
SRC_FILES_PREPRO	=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/*.c)
SRC_FILES_STRUCT	=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/*.c)
SRC_FILES_MAIN		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/*.c)
SRC_FILES_TEST		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/*.c)

#create object file rules for the final binary target
OBJ_FILES_ALGO_OPENMP 		= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c,$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o,$(SRC_FILES_ALGO_OPENMP))
OBJ_FILES_ALGO_CACHE 		= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/cachesim/%.c,$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/cachesim/%.o,$(SRC_FILES_ALGO_CACHE))


OBJ_FILES_UTIL 		= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o,$(SRC_FILES_UTIL))
OBJ_FILES_PREPRO 	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o,$(SRC_FILES_PREPRO))
OBJ_FILES_STRUCT 	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o,$(SRC_FILES_STRUCT))
OBJ_FILES_MAIN  	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/%.o,$(SRC_FILES_MAIN))
OBJ_FILES_TEST  	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/%.o,$(SRC_FILES_TEST))

CC				  = gcc

OPENMP_INC = 	-I$(APP_DIR)/$(INC_DIR)/$(STRUCT_DIR)       \
				-I$(APP_DIR)/$(INC_DIR)/$(ALGO_DIR)/$(INTEGRATION) 	\
				-I$(APP_DIR)/$(INC_DIR)/$(PREPRO_DIR)       \
				-I$(APP_DIR)/$(INC_DIR)/$(UTIL_DIR)


# flags
CFLAGS   =  -O3 -Wall -m64 -fopenmp -g
LFLAGS = -lm -lJudy


# TEST PASS OR FAIL
PASS = \033[92mPASS\033[0m
FAIL = \033[91mFAIL\033[0m

##################################################
##################################################

#########################################################
#       		 ACCEL GRAPH GENERATION    				#
#########################################################

.PHONY: app-openmp
app-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(APP)-openmp
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> AccelGraph: ./$(word 3,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: test-openmp
test-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> TEST: ./$(word 3,$^)"
	@echo " ****************************************************************************** \n"

.PHONY: app-cache
app-cache : directories directories-cache $(APP_DIR)/$(BIN_DIR)/$(APP)-cachesim
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> AccelGraph: ./$(word 3,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: test-cache
test-cache : directories directories-cache $(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-cachesim
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> TEST: ./$(word 3,$^)"
	@echo " ****************************************************************************** \n"

.PHONY: directories
directories :
	@mkdir -p $(APP_DIR)/$(BIN_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)

.PHONY: directories-openmp
directories-openmp :
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp

.PHONY: directories-cache
directories-cache :
	@mkdir -p $(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/cachesim

$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(APP).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(BIN_DIR)/$(APP)-cachesim : $(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(APP).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_CACHE) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(APP_TEST).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(APP).o : $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/$(APP).c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(APP_TEST).o : $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/$(APP_TEST).c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o : $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/cachesim/%.o : $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/cachesim/%.c
	$(CC) $(CFLAGS) $(MODE) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

.PHONY: clean
clean:
	@rm -fr $(APP_DIR)/$(OBJ_DIR)
	@rm -fr $(APP_DIR)/$(BIN_DIR)


.PHONY: clean-stats
clean-stats:
	@rm -f $(BENCHMARKS_DIR)/*.perf
	@rm -f $(BENCHMARKS_DIR)/*/*.perf
	@rm -f $(BENCHMARKS_DIR)/*/*.stats


.PHONY: results-law
results-law:
	@mkdir -p $(APP_DIR)/$(RES_DIR)/results-law_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)
	@for graph in $(LAW); do \
			if [ -e $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).perf ]; then \
				cp $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).perf $(APP_DIR)/$(RES_DIR)/results-law_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)/$$graph; \
			fi ; \
	done

.PHONY: results-mix
results-mix:
	@mkdir -p $(APP_DIR)/$(RES_DIR)/results-mix_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)
	@for graph in $(MIX); do \
			if [ -e $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).perf ]; then \
				cp $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).perf $(APP_DIR)/$(RES_DIR)/results-mix_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)/$$graph; \
			fi ; \
	done

.PHONY: stats-law
stats-law:
	@mkdir -p $(APP_DIR)/$(RES_DIR)/stats-law_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)
	@for graph in $(LAW); do \
			if [ -e $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).stats ]; then \
				cp $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).stats $(APP_DIR)/$(RES_DIR)/stats-law_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)/$$graph; \
			fi ; \
	done

.PHONY: stats-mix
stats-mix:
	@mkdir -p $(APP_DIR)/$(RES_DIR)/stats-mix_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)
	@for graph in $(MIX); do \
			if [ -e $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).stats ]; then \
				cp $(BENCHMARKS_DIR)/$$graph/$(FILE_BIN_TYPE)_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH).stats $(APP_DIR)/$(RES_DIR)/stats-mix_$(ALGORITHMS)_$(DATA_STRUCTURES)_$(NUM_TRIALS)_$(PULL_PUSH)/$$graph; \
			fi ; \
	done


.PHONY: results
results: results-law results-mix

.PHONY: stats
stats: stats-law stats-mix

##################################################
##################################################


##################################################
##################################################

##############################################
#         ACCEL GRAPH TOP LEVEL RULES        #
##############################################

.PHONY: help
help: INTEGRATION = openmp
help: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  --help

.PHONY: run
run: run-openmp

.PHONY: run-cache
run-cache: MODE 		= -DCACHE_HARNESS
run-cache: INTEGRATION 	= cachesim
run-cache: app-cache
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-cachesim  -f $(FILE_BIN) $(ARGS)

.PHONY: run-openmp
run-openmp: INTEGRATION = openmp
run-openmp: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -f $(FILE_BIN) $(ARGS)

.PHONY: convert
convert: INTEGRATION = openmp
convert: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -w -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT)

.PHONY: stats-openmp
stats-openmp: INTEGRATION = openmp
stats-openmp: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -x -f $(FILE_BIN) -n $(NUM_THREADS) -i $(BIN_SIZE) -o $(SORT_TYPE) -l $(REORDER) -p $(INOUT_STATS)

.PHONY: debug-openmp
debug-openmp: INTEGRATION = openmp
debug-openmp: app-openmp
	gdb -ex=r --args ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: debug-memory-openmp
debug-memory-openmp: INTEGRATION = openmp
debug-memory-openmp: app-openmp
	valgrind --leak-check=full --show-leak-kinds=all ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: test-verbose
test-verbose: INTEGRATION = openmp
test-verbose:
	$(MAKE) run-test -C $(APP_DIR) | grep -E 'PASS|FAIL'

# test files
.PHONY: test
test: INTEGRATION = openmp
test:
	@if $(MAKE) run-test -C $(APP_DIR) | grep FAIL > /dev/null; \
		then exit 1; \
	fi

.PHONY: run-test
run-test: run-test-openmp

.PHONY: run-test-openmp
run-test-openmp: test-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp

.PHONY: debug-test-openmp
debug-test-openmp: test-openmp
	gdb ./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp

.PHONY: debug-memory-test-openmp
debug-memory-test-openmp: test-openmp
	valgrind --leak-check=full --show-leak-kinds=all ./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp
# cache performance
.PHONY: cachegrind-perf-openmp
cachegrind-perf-openmp: app-openmp
	valgrind --tool=cachegrind ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: cache-perf
cache-perf-openmp: app-openmp
	sudo perf record -d ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: law
law: app-openmp
	@for graph in $(LAW); do \
	    for num in `seq $(START_THREADS) $(INC_THREADS) $(NUM_THREADS)`; do \
	    	echo $$config $$graph $$num ; \
	    	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -j $(INOUT_STATS) -f $(BENCHMARKS_DIR)/$$graph/graph.wbin -m $(AFU_CONFIG_GENERIC) -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $$num  -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA) ;\
	    done ; \
    done

.PHONY: mix
mix: app-openmp
	@for graph in $(MIX); do \
	    for num in `seq $(START_THREADS) $(INC_THREADS) $(NUM_THREADS)`; do \
	    	echo $$config $$graph $$num ; \
	    	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -j $(INOUT_STATS) -f $(BENCHMARKS_DIR)/$$graph/graph.wbin -m $(AFU_CONFIG_GENERIC) -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $$num  -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA) ;\
	    done ; \
    done

##################################################
##################################################