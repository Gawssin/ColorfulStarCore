# ColorfulStarCore

> Source code for the Paper:


## Structure
- **The Colorful h-star Core Decomposition**
  - *ColorfulStarCoreDecomposition*
- **The h-clique Densest Subgraph Problem**
  - *hCliquePeel*  
  - *OptimizedCliqueCore* 
- **Header Files**
  - *header*

### 1. ColorfulStarCoreDecomposition
> The colorful h-star core decompositon algorithms, containing the basic and advanced version.
#### To compile
```
$ cd ./ColorfulStarCoreDecomposition/
$ g++ -O3 -o HStar ColorfulStarCoreDecompositionWithHeap.cpp
```

#### To run
```
$ ./HStar h filepath [basic]
```

| arguments  | Description |
| :-----| :---- |
| **HStar** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |
| **[basic]** | *Optional*, indicate which algorithm will be performed. "basic" for the basic version **HStarDP**, i.e. recomputing the colorful h-star degree when a neighbor are remusioved. Default for the advanced version **HStarCD**, i.e. using the proposed updating technique to calculate the colorful h-star degree |


### 2. hCliquePeel
> The naive approach to approximate the h-clique densest subgraph following the peeling paradigm, that is to remove the node with minimum h-clique degree and return the subgraph which achieves the highest h-clique density among all subgraphs.
#### To compile
```
$ cd ./hCliquePeel/
$ g++ -O3 -o hCliquePeel hCliquePeeling.cpp
```

#### To run
```
$ ./hCliquePeel h filepath
```

| arguments  | Description |
| :-----| :---- |
| **hCliquePeel** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |

### 3. OptimizedCliqueCore

#### 3.1 The colorful h-star core based algorithm
> *The colorful h-star core based algorithm*: using the colorful h-star θ core as a reduction, and performing the **CoreApp** algorithm or the peeling algorithm on the reduced graph.
#### To compile
```
$ cd ./OptimizedCliqueCore/
$ g++ -O3 -o OptimizedCliqueCore OptimizedCliqueCore.cpp
```

#### To run
```
$ ./OptimizedCliqueCore h filepath
```

| arguments  | Description |
| :-----| :---- |
| **OptimizedCliqueCore** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |


#### 3.2 The heuristic approach (**HStarMPP**)
> using the colorful h-star Kmax core as an approximation solution, or performing the **CoreApp** algorithm or the peeling algorithm on the colorful h-star Kmax core.

#### To compile
```
$ cd ./OptimizedCliqueCore/
$ g++ -O3 -o ColorfulStarKmaxCore ColorfulStarKmaxCore.cpp
```

#### To run
```
$ ./ColorfulStarKmaxCore h filepath alg
```

| arguments  | Description |
| :-----| :---- |
| **OptimizedCliqueCore** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |
| **alg** | "MaxCore" for computing the colorful h-star Kmax core and its h-clique density; "MaxCorePeel" for running the peeling algorithm on the colorful h-star Kmax core |



