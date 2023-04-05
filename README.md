This is the README


Installation:

go to pymodule/ and run "pip install ."

then use "import hypercubes" 


You can also go to pymodule/build and from one of the folders, e.g. "lib.macosx-11.0-arm64-cpython-39/" extract the hypercubes.XXX.so file
and place this in the directory where you would like to import hypercubes

REQUIREMENTS:

C:

* GNU Scientific Library: https://www.gnu.org/software/gsl/ 
  Used for random number generation

Python:

* NumPy
* Matplotlib





# Hypercubes

A Python module written in C to analyse percolation problems on N-dimensional hypercubes. 

A hypercube is a generalisation of a square (dimension 2) or a cube (dimension 3) to any dimension N. A percolation problem is simply stated: we start with a graph (in this case the hypercube) and "activate" each of the edges with some probability p. The question is then whether the hypercube is split up into many small, disconnected pieces ("non-percolating") or is still more-or-less one big connected graph, but with some edges missing ("percolating"). There are some subtleties, but this is the big picture. For more information on percolation, I recommend the book *Introduction to Percolation Theory, D. Stauffer & A. Aharony, Taylor & Francis (2003).*

## Getting Started

Provided you have a local installation of Python, the module can be installed by navigating to pymodule/ and running `pip install .`. The module can then be imported by including `import hypercubes` at the top of your .py file. Alternatively, you can avoid using pip by running `python3 setup.py build` from the pymodule/ directory. This will create a new directory, called build/, in which you will find several subdirectories labelled according to your machine and Python installation. For example, "lib.macosx-11.0-arm64-cpython-39". In one of these folders, you will find a file with a .so extension (for example, "hypercubes.cpython-39-darwin.so"). If you copy this file to the same directory as your Python code, you can import it in the same way, by including the line `import hypercubes`. 

The module is very simple to use, and was designed to perform one task very efficiently: computing the sizes of percolation clusters on the hypercube many times. The core algorithm is a depth-first search which generates the nodes of the hypercube 'on the fly' and enumerates only those nodes which are present in the cluster. The reason for this is that the total number of nodes in the hypercube scales exponentially with the hypercube dimension N. An 18-dimensional hypercube, for example, has over two million edges. Representing the whole graph as an adjacency matrix quickly becomes impractial.

### Prerequisites

C:

* GNU Scientific Library: https://www.gnu.org/software/gsl/ - used for random number generation

Python:

* NumPy
* Matplotlib

## An example calculation

We've included a directory analysis/, in which some key functions relating to percolation are defined, as well as some unit tests. We've also included an example calculation to locate the percolation threshold.




The file `distributions.py` includes functions to calculate key statistical properties



an example calculation to locate the so-called "percolation threshold" is found. This is also 



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

