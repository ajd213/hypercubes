# Hypergraphs

A Python module written in C to analyse percolation problems on hypercubes and other related graphs, such as the Fibonacci cube. 

# Introduction



A hypercube is a generalisation of a square (dimension $2$) or a cube (dimension $3$) to any dimension $N$. A Fibonacci cube (also called the PXP graph) is a subgraph of the Hypercube (see [wikipedia](https://en.wikipedia.org/wiki/Fibonacci_cube)).

Percolation problems are simply stated. First, we choose a probability $p$ (where $0< p < 1$). We then take a graph and "activate" each of the edges independently with the probability $p$. After doing this for all edges, we are interested in whether the graph is split up into many disconnected clusters of sites ("non-percolating") or is still more-or-less one big connected graph, but with some edges missing ("percolating"). If the probability $p$ is small, then it is more likely that there will be many small clusters. Values of $p$ close to $1$, on the other hand, are more likely to leave the graph intact. There are some subtleties, but this is the big picture. 

Percolation theory gets its name from the physical processes of a fluid moving through a porous material, like water through coffee grounds in a percolating coffee machine, or oil in a porous rock bed. 

When a graph is split into many disjoint clusters, we are primarily interested in how large those clusters are. Many properties of percolation, such as the value of $p$ at which the graph switches from the non-percolating to the percolating phase (the so-called "percolation transition"), can be computed if you know how likely clusters of different sizes are to form. Obtaining these probabilities, however, can be computationally expensive, especially for complex graphs such as the hypercube lattice. The most efficient method is often to 'grow' clusters from a starting site, generating the graph and evaluating the probabilities on-the-fly. By growing many clusters, we can estimate the cluster probabilities. 

This process of growing clusters is what the hypergraphs module is for. At its core is an efficient depth-first-search-like algorithm for growing clusters, written in C. The output of this algorithm is a list of many cluster sizes, which can be analysed by the provided functions to compute percolation properties. The module also contains functions for creating Hamiltonians (adjacency matrices) for the two graphs, as well as for running Dijkstra's algorithm.

For more information on percolation, I recommend the book *Introduction to Percolation Theory, D. Stauffer & A. Aharony, Taylor & Francis (2003).*

## Installation

Navigate to pymodule/ and run `pip install .`. 

Alternatively, you can avoid using pip by running `python3 setup.py build` from the pymodule/ directory. This will create a new directory, called build/, in one of the subdirectories, you will find a file with a .so extension, which should be copied to your working directory. The module can then be imported as usual.

Included also is a `main.c` file and associated Makefile, which allows for the complication of the code into a standalone executable run via command line arguments.


### Prerequisites

GNU C, including

* GNU Scientific Library: https://www.gnu.org/software/gsl/ - used for random number generation

Python 3.9 and

* NumPy
* Matplotlib

## Using the module

In the directory analysis/, we include a file distributions.py, which contains functions for analysing cluster sizes and saving/loading data. When analysing the data, be sure to import this file.

### Getting cluster sizes

Let's say we want a hypercube of dimension 16 with percolation concentration $p = 0.5$. To generate the sizes of $10000$ clusters, we run:

```
N = 16
NR = 10000
p = 0.5

cs = hypergraphs.hypercube_clusters(N, NR, p)


```

To compute the mean cluster size S, for example, we then run

```
mean_size = distributions.S(cs)
```

The useful functions `get_clusters_hypercube(N, NR, p, <data_path>)` and `get_clusters_PXP(N, NR, p, <data_path>)` are defined in distributions.py. Given values of N, NR and p, as well as a directory used for data storage, the function will attempt to load the required cluster sizes. If data cannot be found, then it will call the hypergraphs module to generate it, and then save it for next time.

Note that the number of nodes in a hypercube of dimension $N$ is $N_\mathcal{H} = 2^N$, and hence the complexity of the hypergraphs algorithm grows exponentially with $N$. We impose a limit of $N=32$, though computational power will likely limit studies to around $N=20$.

### Creating Hamiltonians

The functions `H_hypercube(N, p)` and `H_PXP(N, p)` generate Hamiltonians (adjacency matrices) for the Hypercube and Fibonacci cube. The RNG is seeded only upon importing the module, so it is safe to create many instances of a Hamiltonian with a particular $p$ value in a short period of time.

```
N = 12
p = 0.5

# create an instance of a hypercube Hamiltonian
H_hyp = hypergraphs.H_hypercube(N, p)

# create an instance of a PXP Hamiltonian
H_p = hypergraphs.H_PXP(N, p)


```

In the directory `analysis/`, the `test_distributions.py` file includes many unit tests pertaining to properties of the Hamiltonians.

### Running Dijkstra's algorithm

The function `hypercube_dijkstra(N, p)` runs an instance of Dijkstra's algorithm from the $0$ site of the hypercube, returning the shortest distances to all of the sites in the same cluster as the $0$ site. By running this function many times, we can plot probability distributions of the shortest paths, which can also be used to diagnose the percolation transition.

```
N = 22
p = 2/N

# run Dijkstra's algorithm
shortest_paths = hypergraphs.hypercube_dijkstra(N, p)

...


```

## An example calculation: locating the percolation transition

It is known analytically that the location of the percolation transition is $p_c = 1/N$ in the limit $N\to\infty$. In this example, we show that the transition is apparent even for modest $N$. The code for this example is located in analysis/percolation_transition.py.

Let's briefly explain the important bits of code. We first create a grid of $p$-values, and choose how many realisations (NR) to use, as well as what values of $N$:
```
p_min = 0
p_max = 1
N_p = 101
plist = np.linspace(p_min, p_max, N_p)

# number of realisations to use
NR = 100000

# system sizes to use
Nlist = np.arange(8, 16, 2)
```

Next, we loop over $N$ and generate and plot the mean cluster size

```
  for Ni, N in enumerate(Nlist):

        NH = 2**N # total number of nodes in the graph

        S, max = s_with_p(N, NR, plist)

        # Plot data ...
```
The function `s_with_p(N, NR, plist)` simply fetches the mean cluster size S and the maximum cluster size for each value of $N$, returning two Numpy arrays. 

In the [sample output](https://github.com/ajd213/hypercubes/blob/master/analysis/example_hypercube_percolation.pdf), we plot the mean size $S$ and the maximum value of $s$ against $p$, and also against $pN$. The plot of $S/N_\mathcal{H}$ against $pN$ (with $N_\mathcal{H}=2^N$ the number of nodes in the graph) reveals most clearly the location of the transition at $pN=1$. We find that for small $p$, $S/N_\mathcal{H}$ is very small, as the average cluster size is a vanishing fraction of the total number of nodes. For $pN>1$, on the other hand, we are in the percolating phase and hence the largest cluster contains a finite fraction of the total number of nodes. $S/N_\mathcal{H}$ is therefore non-zero and increases to $1$ as $p\to 1$. The transition point between these two regimes is represented by a black dashed line at $pN=1$, and becomes sharper as $N$ is increased.

Note that this plot will take some time to reproduce, as many data points must be generated. By lowering the value of NR and the maximum value of $N$, results can be obtained much quicker (due to the exponential scaling with $N$)!

This is just one example: the D. Stauffer & A. Aharony textbook contains many more. 

### Unit tests

Included also in analysis/ is a set of unit tests designed to test both the core hypergraphs code, and the functions included in distributions.py. This includes checking properties of the outputs, as well as comparing the code to some analytic results.

In the case of the Hamiltonians, we check output against third-party code. For the hypercube, this is NetworkX, and for the PXP model a colleague's Python code.

## Full list of functions

### Hypercube

`hypercube_clusters(N, NR, p)`: generate the sizes of NR clusters on a hypercube with percolation concentration $p$, returned as a NumPy array of type NPY_ULONG. Because the sites of the hypercube are all equivalent, we always choose the site `0` as the starting site.

 *  N (int): the dimension of the hypercube
 *  NR (int): the Number of Realisations: number of clusters to grow
 *  p (float): the percolation concentration. $0 <= p <= 1$

`hypercube_H(N, p)`: generate the Hamiltonian (adjacency matrix) for the hypercube, as a NumPy array.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `hypercube_dijkstra(N, p)`: run Dijkstra's algorithm for the hypercube, returning the distances to all of the sites in the same cluster as the $0$ site as a NumPy array. That is, the length of the returned array depends on how many sites are in the cluster.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$


### PXP graph/Fibonacci cube


 `PXP_clusters(N, NR, p)`: generate the sizes of NR clusters on a PXP graph (Fibonacci cube) with percolation concentration $p$, returned as a NumPy array of type NPY_ULONG. Because the sites of the PXP graph are NOT equivalent, the starting site is randomly chosen for each realisation.

 *  N (int): the dimension of the hypercube
 *  NR (int): the Number of Realisations: number of clusters to grow
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `PXP_H(N, p)`: generate the Hamiltonian (adjacency matrix) for the hypercube, as a NumPy array of type NPY_ULONG.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `PXP_Sites(N)`: generate the (ordered) array of nodes for the PXP graph, as a NumPy array of type NPY_ULONGLONG.

  *  N: the dimension of the Fibonacci cube


### Helper functions

`Hamming_distance(state1, state2)`: Compute the Hamming distance between two basis sites (integers). This is the number of bits different between the two sites, i.e., the minimum number of spin-flips required to go from one to the other. 

 *  state1 (int): the first basis site
 *  state2 (int): the second basis site

## Authors

* **Alexander Duthie** [GitHub](https://github.com/ajd213)
