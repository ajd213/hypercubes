# Hypergraphs

A Python module written in C to analyse percolation problems on hypercubes and other related graphs, such as the Fibonacci cube. 

# Introduction



A hypercube is a generalisation of a square (dimension $2$) or a cube (dimension $3$) to any dimension $N$. A Fibonacci cube (also called the PXP graph) is a subgraph of the Hypercube (see [wikipedia](https://en.wikipedia.org/wiki/Fibonacci_cube)).

Percolation problems are simply stated. First, we choose a probability $p$ (where $0< p < 1$). We then take a graph and "activate" each of the edges independently with the probability $p$. After doing this for all edges, we are interested in whether the graph is split up into many disconnected clusters of sites ("non-percolating") or is still more-or-less one big connected graph, but with some edges missing ("percolating"). If the probability $p$ is small, then it is more likely that there will be many small clusters. Values of $p$ close to $1$, on the other hand, are more likely to leave the graph intact.

When a graph is split into many disjoint clusters, we are interested in how large those clusters are. Many properties of percolation, such as the value of $p$ at which the percolation transition takes place, can be computed if you know how likely clusters of different sizes are to form. Obtaining these probabilities, however, can be computationally expensive, especially for complex graphs such as the hypercube lattice. The most efficient method is often to 'grow' clusters from a starting site, generating the graph and evaluating the probabilities on-the-fly. By growing many clusters, we can estimate the cluster probabilities. 

This process of growing clusters is what this module is for. At its core is an efficient depth-first-search-like algorithm for growing clusters, written in C. The output of this algorithm is a list of many cluster sizes, which can be analysed by the provided functions to compute percolation properties. The module has also been extended to generate Hamiltonians (adjacency matrices) and to compute distributions of shortest paths within individual clusters.

For more information on percolation, see the textbook *Introduction to Percolation Theory, D. Stauffer & A. Aharony, Taylor & Francis (2003).*

## Installation

Navigate to src/ and run `pip install .`. 

### Prerequisites

GNU C, including

* GNU Scientific Library: https://www.gnu.org/software/gsl/ - used for random number generation

Python 3.9 and

* NumPy
* Matplotlib

You may need to edit your CPATH and LIBRARY_PATH to include the Python, NumPy and GSL headers. For Python, you can find these by running:

```
from distutils.sysconfig import get_python_inc
get_python_inc()
```

and for NumPy with:
```
import numpy
numpy.get_include()
```

The include path of GSL will depend on where it was installed.

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

Note that the number of nodes in a hypercube of dimension $N$ is $N_\mathcal{H} = 2^N$, and hence the complexity of the hypergraphs algorithm grows exponentially with $N$. We impose a limit of $N=32$, though computational power will likely limit studies before this cutoff is reached.

### Creating Hamiltonians

The functions `hypercube_H(N, p)` and `PXP_H(N, p)` generate Hamiltonians (adjacency matrices) for the Hypercube and Fibonacci cube. The RNG is seeded only upon importing the module, so it is safe to create many instances of a Hamiltonian with a particular $p$ value in a short period of time.

```
N = 12
p = 0.5

# create an instance of a hypercube Hamiltonian
H_hyp = hypergraphs.hypercube_H(N, p)

# create an instance of a PXP Hamiltonian
H_p = hypergraphs.PXP_H(N, p)


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

` hypercube_clusters(N, NR, p)`: generate the sizes of NR clusters on a hypercube with percolation concentration $p$, returned as a NumPy array of type NPY_ULONG.

 *  N (int): the dimension of the hypercube
 *  NR (int): the Number of Realisations: number of clusters to grow
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `hypercube_dijkstra(N, p)`: run Dijksta's algorithm from the root site of the hypercube, and compute the distances
 from the root site to all the other sites in the same cluster. The length of the returned array is equal to the size of
 the cluster containing the $0$ site.

  `hypercube_dijkstra_LC(N, p)`: run Dijksta's algorithm within the largest cluster of a realisation. Return a list of distances from one site within that cluster to all the others. This is more expensive than the standard Dijkstra algorithm, as we must run the algorithm many times over Fock space, until we are sure we have found the largest cluster.

`hypercube_H(N, p)`: generate the Hamiltonian (adjacency matrix) for the hypercube, as a NumPy array.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `hypercube_H_SC(N, p)`: generate the Hamiltonian for a single cluster of the hypercube. Return the Hamiltonian and the cluster size, as a tuple.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 Returns $(H, s)$, where $H$ is the Hamiltonian, and $s$ the size of the cluster.

  `hypercube_H_LC(N, p)`: generate the Hamiltonian for the largest cluster of the hypercube. Return the Hamiltonian and the cluster size, as a tuple. This is more expensive than generating the Hamiltonian for just a single cluster, as we must enumerate all (or most of) Hilbert space to be sure that we have the largest cluster.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 Returns $(H, s)$, where $H$ is the Hamiltonian, and $s$ the size of the cluster.



### PXP graph/Fibonacci cube


 ` PXP_clusters(N, NR, p)`: generate the sizes of NR clusters on a PXP graph (Fibonacci cube) with percolation concentration $p$, returned as a NumPy array of type NPY_ULONG.

 *  N (int): the dimension of the hypercube
 *  NR (int): the Number of Realisations: number of clusters to grow
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `PXP_H(N, p)`: generate the Hamiltonian (adjacency matrix) for the hypercube, as a NumPy array of type NPY_ULONG.

 *  N (int): the dimension of the hypercube
 *  p (float): the percolation concentration. $0 <= p <= 1$

 `PXP_sites(N)`: generate the (ordered) array of nodes for the PXP graph, as a NumPy array of type NPY_ULONGLONG.

  *  N: the dimension of the Fibonacci cube


### Helper functions

`Hamming_distance(state1, state2)`: Compute the Hamming distance between two basis sites (integers). This is the number of bits different between the two sites, i.e., the minimum number of spin-flips required to go from one to the other. 

 *  state1 (int): the first basis site
 *  state2 (int): the second basis site


### Functions in `distributions.py`

Included in `analysis/distributions.py` are functions both for analysing arrays of cluster sizes, as well as generating and saving/loading data. 


## Authors

* **Alexander Duthie** [GitHub](https://github.com/ajd213)
