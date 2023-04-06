# Hypercubes

A Python module written in C to analyse percolation problems on N-dimensional hypercubes. 

# Introduction

A hypercube is a generalisation of a square (dimension $2$) or a cube (dimension $3$) to any dimension $N$. 

Percolation problems are simply stated: we start with a graph (in this case the hypercube) and "activate" each of the edges independently with some probability p. After doing this for all edges, we are interested in whether the hypercube is split up into many small, disconnected clusters of sites ("non-percolating") or is still more-or-less one big connected graph, but with some edges missing ("percolating"). If the probability $p$ is small, then it is more likely that the graph will fragment. Values of $p$ close to $1$, on the other hand, are more likely to leave the graph intact. There are some subtleties, but this is the big picture. 

Percolation theory gets its name from the physical processes of a fluid moving through a porous material, like water through coffee grounds in a percolating coffee machine, or oil in a porous rock bed. 

When a graph is split into many disjoint clusters, we are primarily interested in how large those clusters are. Many properties of percolation, such as the value of $p$ at which the graph switches from the non-percolating to the percolating phase (the so-called "percolation transition"), can be obtained only from the probabilities that clusters of different sizes will form. Obtaining these probabilities, however, can be computationally expensive, especially for complex graphs such as the hypercube lattice. The most efficient method is often to 'grow' clusters from a starting site, generating the graph and evaluating the probabilities on-the-fly. By growing many clusters, we can estimate the cluster probabilities. 

This process of growing clusters is what the hypercubes module is for. At its core is an efficient depth-first-search-like algorithm for growing clusters, written in C. The output of this algorithm is a list of many cluster sizes, which can be analysed by the provided functions to compute percolation properties.

For more information on percolation, I recommend the book *Introduction to Percolation Theory, D. Stauffer & A. Aharony, Taylor & Francis (2003).*

## Installation

Provided you have a local installation of Python, the module can be installed by navigating to pymodule/ and running `pip install .`. The module can then be imported by including `import hypercubes` at the top of your .py file. Alternatively, you can avoid using pip by running `python3 setup.py build` from the pymodule/ directory. This will create a new directory, called build/, in which you will find several subdirectories labelled according to your machine and Python installation. For example, "lib.macosx-11.0-arm64-cpython-39". In one of these folders, you will find a file with a .so extension (for example, "hypercubes.cpython-39-darwin.so"). If you copy this file to the same directory as your Python code, you can import it in the same way, by including the line `import hypercubes`. 


### Prerequisites

C:

* GNU Scientific Library: https://www.gnu.org/software/gsl/ - used for random number generation

Python:

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

cs = hypercubes.clusters(N, NR, p)


```

To compute the mean cluster size S, for example, we then run

```
mean_size = distributions.S(cs)
```



## An example calculation: locating the percolation transition

It is known analytically that the location of the percolation transition is $p_c = 1/N$ in the limit $N\to\infty$. In this example, we show that the transition is apparent even for modest $N$. The code for this example is located in analysis/percolation_transition.py.

We first create a grid of $p$-values, and choose how many realisations (NR) to use, as well as what values of $N$.
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




### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

