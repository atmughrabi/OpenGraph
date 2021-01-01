
#########################################################
#                GENERAL DIRECTOIRES                    #
#########################################################
# globals binaary /bin/accel-graph name doesn't need to match main/accel-graph.c
APP               ?= open-graph

# test name needs to match the file name test/test_accel-graph.c
APP_TEST          ?= test_open-graph-match


# dirs Root app
APP_DIR             ?= .
# BENCHMARKS_DIR        ?= ../../01_GraphDatasets
BENCHMARKS_DIR      ?= ../01_test_graphs

#dir root/managed_folders
SRC_DIR             ?= src
OBJ_DIR             ?= obj
INC_DIR             ?= include
BIN_DIR             ?= bin
RES_DIR             ?= results

#if you want to compile from cmake you need this directory
#cd build
#cmake ..
BUILD_DIR           ?= build

# relative directories used for managing src/obj files
STRUCT_DIR          ?= structures
PREPRO_DIR          ?= preprocess
ALGO_DIR            ?= algorithms
UTIL_DIR            ?= utils

# Folders needed when using CAPI


#contains the tests use make run-test to compile what in this directory
TEST_DIR            ?= tests

#contains the main for the graph processing framework
MAIN_DIR            ?= main

##################################################
##################################################

#########################################################
#                ACCEL RUN GRAPH ARGUMENTS              #
#########################################################

#small test graphs
GRAPH_NAME ?= test
# GRAPH_NAME ?= v51_e1021
# GRAPH_NAME  ?= v300_e2730

# TEXT formant
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph

#UNWEIGHTED
# export FILE_BIN = $(BENCHMARKS_DIR)/$(GRAPH_NAME)/graph.bin

# FILE_BIN_TYPE ?= graph
FILE_BIN_TYPE ?= graph.bin
# FILE_BIN_TYPE ?= graph.wbin

#WEIGHTED
FILE_BIN ?= $(BENCHMARKS_DIR)/$(GRAPH_NAME)/$(FILE_BIN_TYPE)


#GRAPH RUN
SORT_TYPE       ?= 0
REORDER_LAYER1  ?= 0
REORDER_LAYER2  ?= 0
REORDER_LAYER3  ?= 0
DATA_STRUCTURES ?= 0
ALGORITHMS      ?= 1

ROOT            ?= 164
PULL_PUSH       ?= 2
TOLERANCE       ?= 1e-8
DELTA           ?= 800

NUM_THREADS_PRE   ?= $(shell grep -c ^processor /proc/cpuinfo)
NUM_THREADS_ALGO  ?= $(shell grep -c ^processor /proc/cpuinfo)
NUM_THREADS_KER   ?= $(shell grep -c ^processor /proc/cpuinfo)
NUM_ITERATIONS  ?= 1
NUM_TRIALS      ?= 1

FILE_FORMAT     ?= 1
CONVERT_FORMAT  ?= 1

#STATS COLLECTION VARIABLES
BIN_SIZE        ?= 512
INOUT_STATS     ?= 0
MASK_MODE       ?= 0

ARGS ?= -w -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS_PRE) -N $(NUM_THREADS_ALGO) -K $(NUM_THREADS_KER) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -b $(DELTA)

##############################################
#      ACCEL GRAPH COMPILATION VARIABLES     #
##############################################


#each implementation is contained in a separate folder. when you compile for openmp for example the files in openmp folder gets included
SRC_FILES_ALGO_OPENMP       =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/*.c)
SRC_FILES_UTIL      =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/*.c)
SRC_FILES_PREPRO    =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/*.c)
SRC_FILES_STRUCT    =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/*.c)
SRC_FILES_MAIN      =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/*.c)
SRC_FILES_TEST      =   $(wildcard $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/*.c)


INC_FILES_ALGO_OPENMP       =   $(wildcard $(APP_DIR)/$(INC_DIR)/$(ALGO_DIR)/openmp/*.h)
INC_FILES_UTIL      =   $(wildcard $(APP_DIR)/$(INC_DIR)/$(UTIL_DIR)/*.h)
INC_FILES_PREPRO    =   $(wildcard $(APP_DIR)/$(INC_DIR)/$(PREPRO_DIR)/*.h)
INC_FILES_STRUCT    =   $(wildcard $(APP_DIR)/$(INC_DIR)/$(STRUCT_DIR)/*.h)


ALL_HEADER_FILES    =   $(INC_FILES_ALGO_OPENMP) $(INC_FILES_UTIL) $(INC_FILES_PREPRO) $(INC_FILES_STRUCT)

#create object file rules for the final binary target
OBJ_FILES_ALGO_OPENMP = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c,$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o,$(SRC_FILES_ALGO_OPENMP))

OBJ_FILES_UTIL      = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o,$(SRC_FILES_UTIL))
OBJ_FILES_PREPRO    = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o,$(SRC_FILES_PREPRO))
OBJ_FILES_STRUCT    = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o,$(SRC_FILES_STRUCT))
OBJ_FILES_MAIN      = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/%.o,$(SRC_FILES_MAIN))
OBJ_FILES_TEST      = $(patsubst $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/%.c,$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/%.o,$(SRC_FILES_TEST))

CC                = gcc

OPENMP_INC =    -I$(APP_DIR)/$(INC_DIR)/$(STRUCT_DIR)               \
				-I$(APP_DIR)/$(INC_DIR)/$(ALGO_DIR)/$(INTEGRATION)  \
				-I$(APP_DIR)/$(INC_DIR)/$(PREPRO_DIR)               \
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
#                ACCEL GRAPH GENERATION                 #
#########################################################

.PHONY: app-openmp
app-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(APP)-openmp
	@echo "\n ******************************************************************************  "
	@echo " * DONE!! NOTHING ELSE TO COMPILE ---> OpenGraph-openmp: ./$(word 3,$^)"
	@echo " ******************************************************************************  \n"

.PHONY: test-openmp
test-openmp : directories directories-openmp $(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp
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
	@mkdir -p $(APP_DIR)/openmp-results


$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(APP).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp : $(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(APP_TEST).o $(OBJ_FILES_UTIL) $(OBJ_FILES_ALGO_OPENMP) $(OBJ_FILES_PREPRO) $(OBJ_FILES_STRUCT)
	@$(CC) $(CFLAGS) -o $@ $^  $(LFLAGS)

$(APP_DIR)/$(OBJ_DIR)/$(MAIN_DIR)/$(APP).o : $(APP_DIR)/$(SRC_DIR)/$(MAIN_DIR)/$(APP).c $(ALL_HEADER_FILES)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(TEST_DIR)/$(APP_TEST).o : $(APP_DIR)/$(SRC_DIR)/$(TEST_DIR)/$(APP_TEST).c $(ALL_HEADER_FILES)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(UTIL_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(UTIL_DIR)/%.c $(INC_FILES_UTIL)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(ALGO_DIR)/openmp/%.o : $(APP_DIR)/$(SRC_DIR)/$(ALGO_DIR)/openmp/%.c $(ALL_HEADER_FILES)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(PREPRO_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(PREPRO_DIR)/%.c $(INC_FILES_PREPRO)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

$(APP_DIR)/$(OBJ_DIR)/$(STRUCT_DIR)/%.o : $(APP_DIR)/$(SRC_DIR)/$(STRUCT_DIR)/%.c $(INC_FILES_STRUCT)
	$(CC) $(CFLAGS) $(OPENMP_INC) $(MODE) -c -o $@ $<

.PHONY: clean
clean:
	@rm -fr $(APP_DIR)/$(OBJ_DIR)
	@rm -fr $(APP_DIR)/$(BIN_DIR)


.PHONY: clean-stats
clean-stats:
	@rm -f $(BENCHMARKS_DIR)/*.perf
	@rm -f $(BENCHMARKS_DIR)/*/*.perf
	@rm -f $(BENCHMARKS_DIR)/*/*.stats
	@rm -fr $(APP_DIR)/openmp-results

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

.PHONY: run-openmp
run-openmp: INTEGRATION = openmp
run-openmp: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -f $(FILE_BIN) $(ARGS)

.PHONY: convert
convert: INTEGRATION = openmp
convert: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3)

.PHONY: sweep-convert
sweep-convert: INTEGRATION = openmp
sweep-convert: app-openmp
	@for dir in ${BENCHMARKS_DIR}/*;      													\
	do                                      												\
		for graph_dir in $${dir}/*;               											\
		do                                  												\
			if [ -e $${graph_dir}/graph.bin ]; then											\
				echo "convert $${graph_dir}/graph.bin";									    \
				./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -x -c 0 -f $${graph_dir}/graph.bin -z 1 -l 10 -L 10 -O 10;\
				mv $${graph_dir}/graph.bin.txt $${graph_dir}/graph.rand;                    \
				mv $${graph_dir}/graph_Gorder.labels $${graph_dir}/graph_Gorder.org.labels; \
				mv $${graph_dir}/graph_Rabbit.labels $${graph_dir}/graph_Rabbit.org.labels; \
				../../06_Gorder/Gorder $${graph_dir}/graph.rand ;                  							\
				../../05_rabbit_order/demo/reorder $${graph_dir}/graph.rand ;      							\
				mv $${graph_dir}/graph_Gorder.labels $${graph_dir}/graph_Gorder.rand.labels; \
				mv graph_Rabbit.labels $${graph_dir}/graph_Rabbit.rand.labels; \
				./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -x -c 1 -f $${graph_dir}/graph.rand -z 0 ;\
				rm $${graph_dir}/graph.rand;                                                \
			fi; 																			\
		done                                												\
	done

.PHONY: convert-w
convert-w: INTEGRATION = openmp
convert-w: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -w -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3)

.PHONY: stats-openmp
stats-openmp: INTEGRATION = openmp
stats-openmp: app-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -x -f $(FILE_BIN) -n 4 $(NUM_THREADS) -i $(BIN_SIZE) -o $(SORT_TYPE) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -j $(INOUT_STATS)

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
run-test: INTEGRATION = openmp
run-test: run-test-openmp

.PHONY: run-test-openmp
run-test-openmp: INTEGRATION = openmp
run-test-openmp: test-openmp
	./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp >> $(APP_DIR)/openmp-results/full_results_$(ALGORITHMS).perf

.PHONY: debug-test-openmp
debug-test-openmp: INTEGRATION = openmp
debug-test-openmp: test-openmp
	gdb ./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp  

.PHONY: debug-memory-test-openmp
debug-memory-test-openmp: INTEGRATION = openmp
debug-memory-test-openmp: test-openmp
	valgrind --leak-check=full --show-leak-kinds=all ./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp
# cache performance
.PHONY: cachegrind-perf-openmp
cachegrind-perf-openmp: INTEGRATION = openmp
cachegrind-perf-openmp: app-openmp
	valgrind --tool=cachegrind ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

.PHONY: cache-perf-openmp
cache-perf-openmp: INTEGRATION = openmp
cache-perf-openmp: app-openmp
	sudo perf record -d ./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp -f $(FILE_BIN) $(ARGS)

