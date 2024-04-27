[![Build Status](https://travis-ci.com/atmughrabi/OpenGraph.svg?branch=master)](https://travis-ci.com/atmughrabi/OpenGraph)
[<p align="center"><img src="./02_slides/fig/logo.svg" width="200" ></p>](#opengraph-benchmark-suite)

# OpenGraph Benchmark Suite

## Graph Processing Framework with OpenMP

## Overview

![End-to-End Evaluation](./02_slides/fig/theme.png "OpenGraph")

OpenGraph is an open source graph processing framework, designed as a modular benchmarking suite for graph processing algorithms. It provides an end to end evaluation infrastructure which includes the preprocessing stage (forming the graph structure) and the graph algorithm. The OpenMP part of OpenGraph has been developed on Ubuntu 20.04.1, with PowerPC/Intel architecture taken into account.
OpenGraph is coded using C giving the researcher full flexibility with modifying data structures and other algorithmic optimizations.

* Presentations that explains end-to-end graph processing (implementation is inspired from these sources)
  * Preprocessing two steps (third one is optional) :
    1. [[Sorting the edge-list](./02_slides/02_preprocessing_countsort.pdf)], using count-sort or radix-sort.
    2. [[Building the graph structure](./02_slides/03_preprocessing_DataStructures.pdf)]. CSR, Gird, Adjacency-Linked-List, and Adjacency-Array-List.
        * [Ref](https://github.com/thu-pacman/GridGraph): Xiaowei Zhu, Wentao Han and Wenguang Chen. [GridGraph: Large-Scale Graph Processing on a Single Machine Using 2-Level Hierarchical Partitioning](https://www.usenix.org/system/files/conference/atc15/atc15-paper-zhu.pdf). Proceedings of the 2015 USENIX Annual Technical Conference, pages 375-386.
        * [Ref](https://github.com/epfl-labos/EverythingGraph): Malicevic, Jasmina, Baptiste Lepers, and Willy Zwaenepoel. "Everything you always wanted to know about multicore graph processing but were afraid to ask." 2017 USENIX Annual Technical Conference. Proceedings of the 2015 USENIX Annual Technical Conference, pages 375-386.
    3. [[Relabeling the graph](./02_slides/01_algorithm_PR_cache.pdf)], this step achieves better cache locality (better performance) with preprocessing overhead.
        * [Ref](https://github.com/araij/rabbit_order): J. Arai, H. Shiokawa, T. Yamamuro, M. Onizuka, and S. Iwamura. Rabbit Order: Just-in-time Parallel Reordering for Fast Graph Analysis. IEEE International Parallel and Distributed Processing Symposium (IPDPS), 2016.
        * [Ref](https://github.com/faldupriyank/dbg):P. Faldu and J. Diamond and B. Grot, "A Closer Look at Lightweight Graph Reordering," in Proceedings of the International Symposium on Workload Characterization (IISWC), November 2019.
  * Graph Algorithm step depends on the direction of the data (Push/Pull):
    1. [[BFS example](./02_slides/00_algorithm_BFS.pdf)], although it doesn't show direction optimized. But we discusses the Push and Pull approach separately.
        * [[Ref](https://github.com/sbeamer/gapbs)]: Scott Beamer, Krste Asanović, David Patterson. [The GAP Benchmark Suite](http://arxiv.org/abs/1508.03619). arXiv:1508.03619 [cs.DC], 2015.
    2. [[Page-Rank (PR) example](./02_slides/01_algorithm_PR_cache.pdf)]: Discussing PR cache behavior.
       * [Ref](https://github.com/araij/rabbit_order): J. Arai, H. Shiokawa, T. Yamamuro, M. Onizuka, and S. Iwamura. Rabbit Order: Just-in-time Parallel Reordering for Fast Graph Analysis. IEEE International Parallel and Distributed Processing Symposium (IPDPS), 2016.

<!-- ## Details -->
<!-- ### OpenGraph Supported Algorithms -->

# Installation and Dependencies

## CPU Mode

[<img src="./02_slides/fig/openmp_logo.png" height="45" align="right" >](https://www.openmp.org/)
### OpenMP

1. Judy Arrays
```console
open@graph:~$ sudo apt-get install libjudy-dev
```
2. OpenMP is already a feature of the compiler, so this step is not necessary.
```console
open@graph:~$ sudo apt-get install libomp-dev
```

#### Setting up the source code

1. Clone OpenGraph.
```console
open@graph:~$ git clone https://github.com/atmughrabi/OpenGraph.git
```
2. From the home directory go to the OpenGraph directory:
```console
open@graph:~$ cd OpenGraph/
```
3. Make the source code
```console
open@graph:~OpenGraph$ make
```

# Running OpenGraph

## CPU Mode

[<img src="./02_slides/fig/openmp_logo.png" height="45" align="right" >](https://www.openmp.org/)
### Initial compilation for the Graph framework with OpenMP

1. The default compilation is `openmp` mode:
```console
open@graph:~OpenGraph$ make
```
2. From the root directory you can modify the Makefile with the [(parameters)](#OpenGraph-options) you need for OpenMP:
```console
open@graph:~OpenGraph$ make run
```
* OR
```console
open@graph:~OpenGraph$ make run-openmp
```

* You can pass parameters modifying `Makefile` parameters (easiest way) - cross reference with [(parameters)](#OpenGraph-options) to pass the correct values.


| PARAMETER  | FUNCTION | 
| :--- | :--- |
| ARGS  | arguments passed to open-graph |

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *Graph Files Directory* |
| FILE_BIN  | graph edge-list location |
| FILE_LABEL  | graph edge-list reorder list | 

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *Graph Structures PreProcessing* |
| SORT_TYPE  | graph edge-list sort (count/radix) |
| DATA_STRUCTURES  | CSR,GRID,LinkedList,ArrayList |
| REORDER_LAYER1  | Reorder graph for cache optimization |

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *Algorithms General* |
| ALGORITHMS  | BFS, PR, DFS, etc |
| PULL_PUSH  | Direction push,pull,hybrid |

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *Algorithms Specific* |
| ROOT  | source node for BFS, etc |
| TOLERANCE  | PR tolerance for convergence |
| NUM_ITERATIONS  | PR iterations or convergence |
| DELTA  | SSSP delta step |

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *General Performance* |
| NUM_THREADS_PRE  | number of threads for the preprocess step (graph sorting, generation) |
| NUM_THREADS_ALGO  | number of threads for the algorithm step (BFS,PR, etc) |
| NUM_THREADS_KER  | (Optional) number of threads for the algorithm kernel (BFS,PR, etc) |
| NUM_TRIALS  | number of trials for the same algorithms | 


# Graph structure Input (Edge list)

* If you open the Makefile you will see the convention for graph directories : `BENCHMARKS_DIR/GRAPH_NAME/graph.wbin`.
* `.bin` stands to unweighted edge list, `.wbin` stands for wighted, `In binary format`. (This is only a convention you don't have to use it)
* The reason behind converting the edge-list from text to binary, it is simply takes less space on the drive for large graphs, and easier to use with the `mmap` function.

| Source  | Dest | Weight (Optional) |
| :---: | :---: | :---: |
| 30  | 3  |  1 |
| 3  | 4  |  1 |

* Example:
* INPUT: (unweighted textual edge-list)
* ../BENCHMARKS_DIR/GRAPH_NAME/graph
 ```
  30    3
  3     4
  25    5
  25    7
  6     3
  4     2
  6     12
  6     8
  6     11
  8     22
  9     27

 ```
* convert to binary format and add random weights, for this example all the weights are `1`.
* `--graph-file-format` is the type of graph you are reading, `--convert-format` is the type of format you are converting to.
* NOTE: you can read the file from text format without the convert step. By adding `--graph-file-format 0` to the argument list. The default is `1` assuming it is binary. please check `--help` for better explanation.
* `--stats` is a flag that enables conversion. It used also for collecting stats about the graph (but this feature is on hold for now).
* (unweighted graph)
```console
open@graph:~OpenGraph/00_graph_bench$ make convert
```
* OR (weighted graph)
```console
open@graph:~OpenGraph/00_graph_bench$ make convert-w
```
* OR (weighted graph)
```console
open@graph:~OpenGraph/00_graph_bench$ ./bin/open-graph-openmp  --generate-weights --stats --graph-file-format=0 --convert-format=1 --graph-file=../BENCHMARKS_DIR/GRAPH_NAME/graph
```

* `Makefile` parameters

| PARAMETER  | FUNCTION | 
| :--- | :--- |
| *File Formats* |
| FILE_FORMAT  | the type of graph read |
| CONVERT_FORMAT  | the type of graph converted |


* OUTPUT: (weighted binary edge-list)
*  ../BENCHMARKS_DIR/GRAPH_NAME/graph.wbin
```
1e00 0000 0300 0000 0100 0000 0300 0000
0400 0000 0100 0000 1900 0000 0500 0000
0100 0000 1900 0000 0700 0000 0100 0000
0600 0000 0300 0000 0100 0000 0400 0000
0200 0000 0100 0000 0600 0000 0c00 0000
0100 0000 0600 0000 0800 0000 0100 0000
0600 0000 0b00 0000 0100 0000 0800 0000
1600 0000 0100 0000 0900 0000 1b00 0000
0100 0000
```

# Graph Structure Preprocessing:
OpenGraph can handle multiple representations of the graph structure in memory, each has their own theoretical benefits and shortcomings.

## Regular unsorted Edge-list as input.

<p align="center"><img src="./02_slides/fig/datastructures/edgelist-file.png" width ="500" ></p>

## CSR (Compressed Sparse Row)
<p align="center"><img src="./02_slides/fig/datastructures/csr.png" width ="450" ></p>

## Grid
<p align="center"><img src="./02_slides/fig/datastructures/grid.png" width ="700" ></p>

## Array-List
<p align="center"><img src="./02_slides/fig/datastructures/arraylist.png" width="550" ></p>

## Linked-List
<p align="center"><img src="./02_slides/fig/datastructures/linkedlist.png" width="600" ></p>


# OpenGraph Options

```
Usage: open-graph-openmp [OPTION...]
            -f <graph file> -d [data structure] -a [algorithm] -r [root] -n
            [num threads] [-h -c -s -w]

OpenGraph is an open source graph processing framework, it is designed to be a
benchmarking suite for various graph processing algorithms using pure C.

   -a, --algorithm=[DEFAULT:[0]-BFS]

                             [0]-BFS, 
                             [1]-Page-rank, 
                             [2]-SSSP-DeltaStepping,
                             [3]-SSSP-BellmanFord, 
                             [4]-DFS,
                             [5]-SPMV,
                             [6]-Connected-Components,
                             [7]-Betweenness-Centrality, 
                             [8]-Triangle Counting,
                             [9-BUGGY]-IncrementalAggregation.

  -b, --delta=[DEFAULT:1]    
                             SSSP Delta value [Default:1].

  -c, --convert-format=[DEFAULT:[1]-binary-edgeList]

                             [serialize flag must be on --serialize to write]
                             Serialize graph text format (edge list format) to
                             binary graph file on load example:-f <graph file>
                             -c this is specifically useful if you have Graph
                             CSR/Grid structure and want to save in a binary
                             file format to skip the preprocessing step for
                             future runs. 
                             [0]-text-edgeList, 
                             [1]-binary-edgeList,
                             [2]-graphCSR-binary.

  -C, --cache-size=<LLC>     
                             LLC cache size for MASK vertex reodering


  -d, --data-structure=[DEFAULT:[0]-CSR]

                             [0]-CSR, 
                             [1]-Grid, 
                             [2]-Adj LinkedList, 
                             [3]-Adj ArrayList 
                             [4-5] same order bitmap frontiers.

  -e, --tolerance=[EPSILON:0.0001]

                             Tolerance value of for page rank
                             [default:0.0001].

  -f, --graph-file=<FILE>    

                             Edge list represents the graph binary format to
                             run the algorithm textual format change
                             graph-file-format.

  -F, --labels-file=<FILE>   
                             Read and reorder vertex labels from a text file,
                             Specify the file name for the new graph reorder,
                             generated from Gorder, Rabbit-order, etc.

  -g, --bin-size=[SIZE:512]  
                             You bin vertices's histogram according to this
                             parameter, if you have a large graph you want to
                             illustrate.

  -i, --num-iterations=[DEFAULT:20]

                             Number of iterations for page rank to converge
                             [default:20] SSSP-BellmanFord [default:V-1].

  -j, --verbosity=[DEFAULT:[0:no stats output]

                             For now it controls the output of .perf file and
                             PageRank .stats (needs --stats enabled)
                             filesPageRank .stat [1:top-k results] [2:top-k
                             results and top-k ranked vertices listed.

  -k, --remove-duplicate     
                             Removers duplicate edges and self loops from the
                             graph.

  -K, --Kernel-num-threads=[DEFAULT:algo-num-threads]

                             Number of threads for graph processing kernel
                             (critical-path) (graph algorithm)

  -l, --light-reorder-l1=[DEFAULT:[0]-no-reordering]

                             Relabels the graph for better cache performance
                             (first layer). 
                             [0]-no-reordering, 
                             [1]-out-degree,
                             [2]-in-degree, 
                             [3]-(in+out)-degree, 
                             [4]-DBG-out,
                             [5]-DBG-in, 
                             [6]-HUBSort-out, 
                             [7]-HUBSort-in,
                             [8]-HUBCluster-out, 
                             [9]-HUBCluster-in,
                             [10]-(random)-degree,  
                             [11]-LoadFromFile (used for Rabbit order).

  -L, --light-reorder-l2=[DEFAULT:[0]-no-reordering]

                             Relabels the graph for better cache performance
                             (second layer). 
                             [0]-no-reordering, 
                             [1]-out-degree,
                             [2]-in-degree, 
                             [3]-(in+out)-degree, 
                             [4]-DBG-out,
                             [5]-DBG-in, 
                             [6]-HUBSort-out, 
                             [7]-HUBSort-in,
                             [8]-HUBCluster-out, 
                             [9]-HUBCluster-in,
                             [10]-(random)-degree,  
                             [11]-LoadFromFile (used for Rabbit order).

 -O, --light-reorder-l3=[DEFAULT:[0]-no-reordering]

                             Relabels the graph for better cache performance
                             (third layer). 
                             [0]-no-reordering, 
                             [1]-out-degree,
                             [2]-in-degree, 
                             [3]-(in+out)-degree, 
                             [4]-DBG-out,
                             [5]-DBG-in, 
                             [6]-HUBSort-out, 
                             [7]-HUBSort-in,
                             [8]-HUBCluster-out, 
                             [9]-HUBCluster-in,
                             [10]-(random)-degree,  
                             [11]-LoadFromFile (used for Rabbit order).

  -M, --mask-mode=[DEFAULT:[0:disabled]]

                             Encodes [0:disabled] the last two bits of
                             [1:out-degree]-Edgelist-labels
                             [2:in-degree]-Edgelist-labels or
                             [3:out-degree]-vertex-property-data
                             [4:in-degree]-vertex-property-data with hot/cold
                             hints [11:HOT]|[10:WARM]|[01:LUKEWARM]|[00:COLD]
                             to specialize caching. The algorithm needs to
                             support value unmask to work.

  -n, --pre-num-threads=[DEFAULT:MAX]

                             Number of threads for preprocessing (graph
                             structure) step 

  -N, --algo-num-threads=[DEFAULT:MAX]

                             Number of threads for graph processing (graph
                             algorithm)

  -o, --sort=[DEFAULT:[0]-radix-src]

                             [0]-radix-src, 
                             [1]-radix-src-dest, 
                             [2]-count-src,
                             [3]-count-src-dst.



  -p, --direction=[DEFAULT:[0]-PULL]

                             [0]-PULL, 
                             [1]-PUSH,
                             [2]-HYBRID. 

                             NOTE: Please consult the function switch table for each
                             algorithm.

  -r, --root=[DEFAULT:0]     
                             BFS, DFS, SSSP root

  -s, --symmetrize           
                             Symmetric graph, create a set of incoming edges.

  -S, --stats                
                             Write algorithm stats to file. same directory as
                             the graph.PageRank: Dumps top-k ranks matching
                             using QPR similarity metrics.

  -t, --num-trials=[DEFAULT:[1 Trial]]

                             Number of trials for whole run (graph algorithm
                             run) [default:1].

  -w, --generate-weights     
                             Load or Generate weights. Check ->graphConfig.h
                             #define WEIGHTED 1 beforehand then recompile using
                             this option.

  -x, --serialize            
                             Enable file conversion/serialization use with
                             --convert-format.

  -z, --graph-file-format=[DEFAULT:[1]-binary-edgeList]

                             Specify file format to be read, is it textual edge
                             list, or a binary file edge list. This is
                             specifically useful if you have Graph CSR/Grid
                             structure already saved in a binary file format to
                             skip the preprocessing step. 
                             [0]-text edgeList,
                             [1]-binary edgeList, 
                             [2]-graphCSR binary.

  -?, --help                 Give this help list
      --usage                Give a short usage message
  -V, --version              Print program version


```


# Organization

* `00_graph_bench`
  * `include` - Major function headers
    * `graphalgorithms` - supported Graph algorithms
      * `openmp`  - OpenMP integration
        * `BFS.h`   - Breadth First Search
        * `DFS.h`   - Depth First Search
        * `SSSP.h`  - Single Source Shortest Path
        * `bellmanFord.h` - Single Source Shortest Path using Bellman Ford
        * `incrementalAgreggation.h` - Incremental Aggregation for clustering
        * `pageRank.h` - Page Rank Algorithm
        * `SPMV.h` - Sparse Matrix Vector Multiplication
    * `preprocessing` - preprocessing graph structure
      * `countsort.h` - sort edge list using count sort
      * `radixsort.h` - sort edge list using radix sort
      * `reorder.h` - cluster reorder the graph for better cache locality
      * `sortRun.h` - chose which sorting algorithm to use
    * `structures` - structures that hold the graph in memory
      * `graphAdjArrayList.h` - graph using adjacency list array with arrays
      * `graphAdjLinkeList.h` - graph using adjacency list array with linked lists
      * `graphCSR.h` - graph using compressed sparse matrix
      * `graphGrid.h` - graph using Grid
  * `src` - Major function Source files
    * `graphalgorithms` - supported Graph algorithms
      * `openmp`  - OpenMP integration
        * `BFS.c`   - Breadth First Search
        * `DFS.c`   - Depth First Search
        * `SSSP.c`  - Single Source Shortest Path
        * `bellmanFord.c` - Single Source Shortest Path using Bellman Ford
        * `incrementalAgreggation.c` - Incremental Aggregation for clustering
        * `pageRank.c` - Page Rank Algorithm
        * `SPMV.c` - Sparse Matrix Vector Multiplication
    * `preprocessing` - preprocessing graph structure
      * `countsort.c` - sort edge list using count sort
      * `radixsort.c` - sort edge list using radix sort
      * `reorder.c` - cluster reorder the graph for better cache locality
      * `sortRun.c` - chose which sorting algorithm to use
    * `structures` - structures that hold the graph in memory
      * `graphAdjArrayList.c` - graph using adjacency list array with arrays
      * `graphAdjLinkeList.c` - graph using adjacency list array with linked lists
      * `graphCSR.c` - graph using compressed sparse matrix
      * `graphGrid.c` - graph using Grid

* *`Makefile`* - Global makefile

# Tasks TODO CSR Graphs Only:

- [x] Finish preprocessing sort
  - [x] Radix sort
  - [x] Count sort
  - [ ] Bitonic sort
- [x] Finish preprocessing Graph Data-structures
  - [x] CSR   (Compressed Sparse Row)
  - [x] Grid
  - [x] Adjacency Linked List
  - [x] Adjacency Array List
- [x] Add Light weight reordering
- [x] Finish graph algorithms suite OpenMP
  - [x] BFS   (Breadth First Search)
  - [x] PR    (Page-Rank)
  - [x] DFS   (Depth First Search)
  - [x] IA    (Incremental Aggregation) BUGGY*
  - [x] SSSP  (BellmanFord)
  - [x] SSSP  (Delta Stepping)
  - [x] SPMV  (Sparse Matrix Vector Multiplication)
  - [x] CC    (Connected Components)
  - [x] TC    (Triangle Counting)
  - [x] BC    (Betweenness Centrality)
- [x] Support testing


Report bugs to:
- <atmughrabi@gmail.com>
- <atmughra@ncsu.edu>
[<p align="right"> <img src="./02_slides/fig/logo.svg" width="100" ></p>](#OpenGraph-benchmark-suite)