

#########################################################
#       		 GENERAL DIRECTOIRES   	    			#
#########################################################
# globals binaary /bin/open-graph name doesn't need to match main/open-graph.c
GAPP               ?= open-graph

# test name needs to match the file name test/test_accel-graph.c
GAPP_TEST          = test_open-graph

# dirs Root app 
APP_DIR           	?= .
SCRIPT_DIR          ?= ../04_scripts
BENCHMARKS_DIR    	?= ../03_test_graphs

#dir root/managed_folders
SRC_DIR           	?= src
OBJ_DIR			  	?= obj
INC_DIR			  	?= include
BIN_DIR			  	?= bin

#if you want to compile from cmake you need this directory
#cd build
#cmake ..
BUILD_DIR		  	?= build

# relative directories used for managing src/obj files
STRUCT_DIR		  	?= structures
PREPRO_DIR		  	?= preprocessing
ALGO_DIR		  	?= graphalgorithms
UTIL_DIR		  	?= utils

#contains the tests use make run-test to compile what in this directory
TEST_DIR		  	?= tests

#contains the main for the graph processing framework
MAIN_DIR		  	?= main

##################################################
##################################################

#########################################################
#       		 OPEN RUN GRAPH ARGUMENTS    			#
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


#UNWEIGHTED
# FILE_BIN ?= $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.bin

#WEIGHTED
FILE_BIN ?= $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.wbin



#GRAPH RUN
SORT_TYPE 		?= 0
REORDER 		?= 0
DATA_STRUCTURES ?= 0
ALGORITHMS 		?= 1

ROOT 			?= 164
PULL_PUSH 		?= 2
TOLERANCE 		?= 1e-8
DELTA 			?= 800

NUM_THREADS  	?= 4
# NUM_THREADS  	= $(shell grep -c ^processor /proc/cpuinfo)
NUM_ITERATIONS 	?= 1
NUM_TRIALS 		?= 1

FILE_FORMAT 	?= 1
CONVERT_FORMAT 	?= 1

#STATS COLLECTION VARIABLES
BIN_SIZE ?= 512
INOUT_STATS ?= 2

ARGS ?= -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER) -b $(DELTA)
##################################################
##################################################


##############################################
#      OPEN GRAPH COMPILATION VARIABLES     #
##############################################


#each implementation is contained in a separate folder. when you compile for openmp for example the files in openmp folder gets included
SRC_FILES_ALGO_OPENMP		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/*.c)

SRC_FILES_UTIL		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/*.c) 
SRC_FILES_PREPRO	=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/*.c) 
SRC_FILES_STRUCT	=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/*.c)
SRC_FILES_MAIN		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/*.c)
SRC_FILES_TEST		=	$(wildcard $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/*.c)

#create object file rules for the final binary target
OBJ_FILES_ALGO_OPENMP 		= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c,$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o,$(SRC_FILES_ALGO_OPENMP)) 

OBJ_FILES_UTIL 		= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o,$(SRC_FILES_UTIL)) 
OBJ_FILES_PREPRO 	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o,$(SRC_FILES_PREPRO)) 
OBJ_FILES_STRUCT 	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o,$(SRC_FILES_STRUCT)) 
OBJ_FILES_MAIN  	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/%.o,$(SRC_FILES_MAIN)) 
OBJ_FILES_TEST  	= $(patsubst $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/%.o,$(SRC_FILES_TEST)) 

CC				  = gcc 

OPENMP_INC = 	-I$(APP_DIR)/$(INC_DIR)/$(STRUCT_DIR)   \
				-I$(APP_DIR)/$(INC_DIR)/$(ALGO_DIR)/openmp 	\
				-I$(APP_DIR)/$(INC_DIR)/$(PREPRO_DIR)   \
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
#       		 OPEN GRAPH GENERATION    				#
#########################################################

.PHONY: graph-openmp
graph-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> OpenGraph: ./$(word 3,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: test-openmp
test-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(GAPP_TEST)-openmp
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

$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(GAPP).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)
	
$(APP_DIR)/$(BIN_DIR)/$(GAPP_TEST)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(GAPP_TEST).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(GAPP).o : $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/$(GAPP).c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(GAPP_TEST).o : $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/$(GAPP_TEST).c 
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o : $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c
	$(CC) $(CFLAGS) $(OPENMP_INC) -c -o $@ $<

.PHONY: clean
clean:
	@rm -fr $(APP_DIR)/$(OBJ_DIR)
	@rm -fr $(APP_DIR)/$(BIN_DIR)

##################################################
##################################################

##############################################
#         OPEN GRAPH TOP LEVEL RULES        #
##############################################

.PHONY: help
help: graph-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp  --help

.PHONY: run
run: run-openmp

.PHONY: run-openmp
run-openmp: graph-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp  -f $(FILE_BIN) $(ARGS)

.PHONY: convert
convert: graph-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp  -w -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT)

.PHONY: stats-openmp
stats-openmp: graph-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp -x -f $(FILE_BIN) -n $(NUM_THREADS) -i $(BIN_SIZE) -o $(SORT_TYPE) -l $(REORDER) -p $(INOUT_STATS)

.PHONY: debug-openmp
debug-openmp: graph-openmp	
	gdb -ex=r --args ./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: debug-memory-openmp
debug-memory-openmp: graph-openmp	
	valgrind --leak-check=full --show-leak-kinds=all ./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: test-verbose
test-verbose:
	$(MAKE) run-test -C $(APP_DIR) | grep -E 'PASS|FAIL'

# test files
.PHONY: test
test:
	@if $(MAKE) run-test -C $(APP_DIR) | grep FAIL > /dev/null; \
		then exit 1; \
	fi
	
.PHONY: run-test
run-test: run-test-openmp

.PHONY: run-test-openmp
run-test-openmp: test-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(GAPP_TEST)-openmp

.PHONY: debug-test-openmp
debug-test-openmp: test-openmp	
	gdb ./$(APP_DIR)/$(BIN_DIR)/$(GAPP_TEST)-openmp

.PHONY: debug-memory-test-openmp
debug-memory-test-openmp: test-openmp	
	valgrind --leak-check=full --show-leak-kinds=all ./$(APP_DIR)/$(BIN_DIR)/$(GAPP_TEST)-openmp
# cache performance
.PHONY: cachegrind-perf-openmp
cachegrind-perf-openmp: graph-openmp
	valgrind --tool=cachegrind ./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp -f $(FILE_BIN) $(ARGS)
.PHONY: cache-perf
cache-perf-openmp: graph-openmp
	sudo perf record -d ./$(APP_DIR)/$(BIN_DIR)/$(GAPP)-openmp -f $(FILE_BIN) $(ARGS)

##################################################
##################################################
