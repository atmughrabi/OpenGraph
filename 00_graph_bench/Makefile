
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
CONFIG_DIR			?= config



#contains the tests use make run-test to compile what in this directory
TEST_DIR            ?= tests

#contains the main for the graph processing framework
MAIN_DIR            ?= main

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
export CACHE_SIZE       ?= 32768

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

export ARGS ?= -k -M $(MASK_MODE) -j $(INOUT_STATS) -g $(BIN_SIZE) -z $(FILE_FORMAT) -d $(DATA_STRUCTURES) -a $(ALGORITHMS) -r $(ROOT) -n $(NUM_THREADS_PRE) -N $(NUM_THREADS_ALGO) -K $(NUM_THREADS_KER) -i $(NUM_ITERATIONS) -o $(SORT_TYPE) -p $(PULL_PUSH) -t $(NUM_TRIALS) -e $(TOLERANCE) -F $(FILE_LABEL) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -b $(DELTA) -C $(CACHE_SIZE)

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
INC_FILES_CONFIG    =   $(wildcard $(APP_DIR)/$(INC_DIR)/$(CONFIG_DIR)/*.h)


ALL_HEADER_FILES    =   $(INC_FILES_ALGO_OPENMP) $(INC_FILES_UTIL) $(INC_FILES_PREPRO) $(INC_FILES_STRUCT) $(INC_FILES_CONFIG)

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
				-I$(APP_DIR)/$(INC_DIR)/$(UTIL_DIR)   				\
				-I$(APP_DIR)/$(INC_DIR)/$(CONFIG_DIR)              


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
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -s -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -F $(FILE_LABEL)

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
	./$(APP_DIR)/$(BIN_DIR)/$(APP)-openmp  -w -x -c $(CONVERT_FORMAT) -f $(FILE_BIN) -z $(FILE_FORMAT) -l $(REORDER_LAYER1) -L $(REORDER_LAYER2) -O $(REORDER_LAYER3) -F $(FILE_LABEL)

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
	./$(APP_DIR)/$(BIN_DIR)/$(APP_TEST)-openmp >> $(APP_DIR)/openmp-results/full_results_$(NUM_THREADS_KER).perf

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

