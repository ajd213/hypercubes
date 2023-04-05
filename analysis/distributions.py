"""Functions to compute statistical properties of percolation clusters. 

For a good introduction to the theory of percolation, see:
Introduction to Percolation Theory, D. Stauffer & A. Aharony, Taylor & Francis (2003). 
"""

import numpy as np
import hypercubes


def cluster_numbers(cs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """From an array of cluster sizes, compute the cluster numbers.
    
    The cluster number n_s for a particular cluster size s is the number of s-clusters
    per lattice site. As a larger cluster is more likely to be grown (as it contains more
    potential starting sites), we must account for this weighting. Hence:

    n_s = #s-clusters / (s * #clusters)

    Parameters:

    cs : a numpy array of cluster sizes.

    Returns:

    s : an array of cluster sizes
    n_s : an array of cluster numbers. s[i] is the size of a cluster with cluster number n_s[i].
    """
    N_CLUSTERS = len(cs)

    # how many (counts) of each cluster size s are there?
    s, counts = np.unique(cs, return_counts=True) 
    n_s = counts/(s*N_CLUSTERS)

    return s, n_s


def S(cs: np.ndarray) -> np.floating:
    """Compute the mean cluster size.
    
    This mean is the mean cluster size obtained by repeatedly picking a site
    at random and computing the size of the cluster to which it belongs. This
    is exactly what the code in hypercubes does.
    
    Parameters:

    cs : a numpy array of cluster sizes.

    Returns:

    S : the mean cluster size
    """
    return np.mean(cs)


def S_prime(cs: np.ndarray) -> np.floating:
    """Compute the alternative mean cluster size.
    
    This mean is the meav cluster size obtained by averaging the sizes
    of all clusters present in the graph. It is therefore the average 
    of n_s weighted by s.
    
    Parameters:

    cs : a numpy array of cluster sizes.

    Returns:

    S' : the alternative mean cluster size
    """
    s, n_s = cluster_numbers(cs)
    return np.sum(s*n_s)/np.sum(n_s)



def w_s(cs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """The probability of a randomly-chosen occupied site belonging to as s-cluster.
    
    w_s is the probability that if we choose a site at random, it belongs to a cluster
    of size s. This is straightforward to obtain from the list of cluster sizes, as it
    is the number of s-clusters divided by the total number of clusters.

    Parameters:

    cs : a numpy array of cluster sizes.

    Returns:

    s : an array of cluster sizes
    w_s : an array of w_s values, which are the probabilities that a randomly-chosen site
          is part of an s-cluster. s[i] is the size of a cluster with w_s value w_s[i].
    """
    N_CLUSTERS = len(cs)
    s, counts = np.unique(cs, return_counts=True)
    w_s = counts/N_CLUSTERS

    return s, w_s


def get_clusters(N : int, NR : int, p : float, data_path : str) -> np.ndarray:
    """ Attempt to load cluster data, else generate and save it."""
    name = f"clusters_N{N}_NR{NR}_p{p:.4f}.npy"

    try:
        clusters = np.load(data_path + name)

    except FileNotFoundError:
        print(f"Generating data: {name}")
        clusters = hypercubes.clusters(int(N), int(NR), p)
        np.save(data_path + name, clusters)
    
    return np.array(clusters)
