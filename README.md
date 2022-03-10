# ColorfulStarCore

> Source code for the Paper: "Colorful *h*-star Core Decomposition"

> The default graph coloring algorithm: **Degree**

## Structure
- **The Colorful h-star Core Decomposition**
  - *ColorfulStarCoreDecomposition*
- **The h-clique Densest Subgraph Problem**
  - *hCliquePeel*  
  - *OptimizedCliqueCore* 
- **Header Files**
  - *header*

## 1. ColorfulStarCoreDecomposition
> The colorful h-star core decompositon algorithms, including the basic and advanced version.
### To compile
```
$ cd ./ColorfulStarCoreDecomposition/
$ g++ -O3 -o HStar ColorfulStarCoreDecomposition.cpp
```

### To run
```
$ ./HStar h filepath colorAlgo [basic]
```

| arguments  | Description |
| :-----| :---- |
| **HStar** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |
| **colorAlgo** | graph coloring algorithm: <br>"0" for **Degree**; <br>"1" for **Degen**; <br>"2" for **FF**; <br>"3" for **SD**. |
| **[basic]** | *Optional*, indicate which algorithm will be performed. <br>"basic" for the basic version **HStarDP**, i.e. recomputing the colorful h-star degree when a neighbor is removed; <br>the default is the advanced version **HStarCD**, i.e. using the proposed updating technique to calculate the colorful h-star degree. |


## 2. hCliquePeel
> The naive approach to approximate the h-clique densest subgraph following the peeling paradigm, that is to remove the node with minimum h-clique degree and return the subgraph which achieves the highest h-clique density among all subgraphs.
### To compile
```
$ cd ./hCliquePeel/
$ g++ -O3 -o hCliquePeel hCliquePeeling.cpp
```

### To run
```
$ ./hCliquePeel h filepath
```

| arguments  | Description |
| :-----| :---- |
| **hCliquePeel** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |

## 3. OptimizedCliqueCore

### 3.1 The colorful h-star core based algorithm (**HStarPP**)
> using the colorful h-star θ core as a reduction, and performing the **CoreApp** algorithm or the peeling algorithm on the reduced subgraph.
### To compile
```
$ cd ./OptimizedCliqueCore/
$ g++ -O3 -o HStarPP OptimizedCliqueCore.cpp
```

### To run
```
$ ./HStarPP h filepath
```

| arguments  | Description |
| :-----| :---- |
| **HStarPP** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |


### 3.2 The heuristic approach (**HStarMPP**)
> using the colorful h-star Kmax core as an approximation solution, or performing the **CoreApp** algorithm or the peeling algorithm on the colorful h-star Kmax core.

### To compile
```
$ cd ./OptimizedCliqueCore/
$ g++ -O3 -o HStarMPP ColorfulStarKmaxCore.cpp
```

### To run
```
$ ./HStarMPP h filepath alg
```

| arguments  | Description |
| :-----| :---- |
| **HStarMPP** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |
| **alg** | "MaxCore" for computing the colorful h-star Kmax core and its h-clique density; <br> "MaxCorePeel" for running the peeling algorithm on the colorful h-star Kmax core. |



### 3.3 Binary search for the colorful h-star Kmax core (**HStarMB**)
> compute the the colorful h-star Kmax core directly using binary search method, return the Kmax and the Kmax core

### To compile
```
$ cd ./OptimizedCliqueCore/
$ g++ -O3 -o HStarMB ColorfulStarKmaxCore-BS.cpp
```

### To run
```
$ ./HStarMB h filepath
```

| arguments  | Description |
| :-----| :---- |
| **HStarMB** | executable file |
| **h** | the size of stars |
| **filepath** | input file path |


