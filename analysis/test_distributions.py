import unittest
import distributions
import hypergraphs
import numpy as np
import os, shutil
from scipy.special import comb
import networkx as nx

# where to save the temp data for testing
DATA_PATH = "./data/testing/"

class Testdistributions(unittest.TestCase):

    def setUp(self):
        if not os.path.exists(DATA_PATH):
            os.makedirs(DATA_PATH)

    def tearDown(self):
        # delete the contents of the data/testing directory
        for filename in os.listdir(DATA_PATH):
            file_path = os.path.join(DATA_PATH, filename)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print('Failed to delete %s. Reason: %s' % (file_path, e))

    def test_cluster_properties(self):

        # test parameters
        N = 10
        NH = 2**N
        NR = 111

        # first test an intermediate value of p. Check that the right
        # number of clusters are generated.
        p = 0.5
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        self.assertEqual(len(cs), NR)

        # for p=0, check that all the clusters are of size 1.
        p = 0
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        np.testing.assert_equal(cs, 1)

        # for p=1, check that all clusters span the whole graph.
        p = 1
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        np.testing.assert_equal(cs, NH)

    def test_cluster_numbers(self):

        # test parameters
        N = 8
        NH = 2**N
        NR = 97

        # test p = 0: minimally connected
        p = 0
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, n_s = distributions.cluster_numbers(cs)

        # for p=0, there should only be one cluster number
        self.assertEqual(len(n_s), 1)

        # its value should be 1, as there is one 1-cluster per lattice site
        self.assertEqual(n_s[0], 1)


        # test p = 1: maximally connected
        p = 1
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, n_s = distributions.cluster_numbers(cs)

        # again, there should only be one cluster number, as every cluster
        # is equal in size to NH, the total number of lattice nodes
        self.assertEqual(len(n_s), 1)

        # its value should be 1/NH
        self.assertEqual(n_s[0], 1/NH)


        # test p = 0.5
        p = 0.5
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, n_s = distributions.cluster_numbers(cs)

        # for intermediate p, we check normalisation. This essentially says
        # that every site is in a cluster of some size.
        # for bond percolation, the rule is: \sum_s s*n_s = 1
        np.testing.assert_almost_equal(np.sum([s*n_s]), 1, 4)


        # For N = 2, test the cluster numbers against analytic results
        # We need to use many realisations to achieve convergence.

        NR = 1000000
        N = 2
        NH = 2**N
        p = 0.8
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, n_s = distributions.cluster_numbers(cs)


        # analytic results for cluster numbers on a square
        n_1 = (1-p)**2
        n_2 = p*((1-p)**2)
        n_3 = (p**2)*((1-p)**2)
        n_4 = p**3 - (3/4)*(p**4)

        # aim for numerical calculations to be within one percent of the analytics.
        # as this is a random process, it might sometimes fail by a small amount
        np.testing.assert_allclose(np.array([n_1, n_2, n_3, n_4]), np.array([n_s[0], n_s[1], n_s[2], n_s[3]]), rtol = 0.01)

    def test_w_s(self):

        # test parameters
        N = 12
        NR = 100
        p = 0.49

        # test for normalisation. w_s is a probability, and so its 
        # values should sum to 1.
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, w_s = distributions.w_s(cs)
        self.assertAlmostEqual(1, sum(w_s), places=5)


        # test p = 0
        p = 0
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, w_s = distributions.w_s(cs)

        # when p=0, there should only be the value s=1 returned
        self.assertEqual(s[0], 1)
        self.assertEqual(len(s), 1)

        # when p=0, the probability of ending up in a 1-cluster is 1.
        self.assertEqual(w_s[0], 1)
        self.assertEqual(len(w_s), 1)


        # test p = 1
        p = 1
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        s, w_s = distributions.w_s(cs)

        # when p=1, the probability of ending up in a 2**N-cluster is 1.
        self.assertEqual(s[0], 2**N)
        self.assertEqual(len(s), 1)
        self.assertEqual(w_s[0], 1)
        self.assertEqual(len(w_s), 1)

    def test_S(self):

        N = 7
        NR = 99

        # check the result that the mean cluster size is equal to 1
        # when p=0.
        p = 0
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        S = distributions.S(cs)
        self.assertEqual(S, 1)


        # when p=1, the mean cluster size should be NH = 2**N, the total number
        # of nodes
        p = 1
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        S = distributions.S(cs)
        self.assertEqual(S, 2**N)

    def test_PXP_clusters(self):
        
        # check that for p=1, the sizes of the clusters for a given N are equal to F(N+2)
        NR = 7
        p = 1

        # the first 25 Fibonacci numbers
        fibs = np.array([0,1,1,2,3,5,8,13,21,34,55,89,144,233,377,610,987,1597,2584,4181,6765,10946,17711,28657,46368])

        for index in range(1, 23):
            N = index
            cs = distributions.get_clusters_PXP(N, NR, p, DATA_PATH)
            np.testing.assert_equal(cs, fibs[N+2])
        
        # check that for p=0, the cluster size should be 1
        NR = 16
        p = 0
        N = 12
        cs = distributions.get_clusters_PXP(N, NR, p, DATA_PATH)
        np.testing.assert_equal(cs, 1)

        # check that the right number of clusters are produced for intermediate p
        NR = 111
        N = 7
        p = 0.25
        cs = distributions.get_clusters_PXP(N, NR, p, DATA_PATH)
        self.assertEqual(len(cs), NR)


        # test w_s for intermediate p
        N = 12
        NR = 100
        p = 0.49

        # test for normalisation. w_s is a probability, and so its 
        # values should sum to 1.
        cs = distributions.get_clusters_PXP(N, NR, p, DATA_PATH)
        s, w_s = distributions.w_s(cs)
        self.assertAlmostEqual(1, sum(w_s), places=5)

    def test_H(self):

        Nlist = range(1, 11)
        HS_DIM_LIST = np.power(2, Nlist)


        # first test p = 1. The Hamiltonian should represent a complete hypercube
        p = 1
        for j in range(len(Nlist)):
            N = Nlist[j]
            NH = HS_DIM_LIST[j]
            H = hypergraphs.H_hypercube(N, p)

            # ensure that the Hamiltonian matrix has dim 2**N by 2**N
            np.testing.assert_equal(H.shape, (NH, NH))

            # ensure that all rows and all columns sum to N
            np.testing.assert_equal(np.fromiter((np.sum(row) for row in H), dtype=int), N)
            np.testing.assert_equal(np.fromiter((np.sum(row) for row in H.T), dtype=int), N)

        
        # next, let's use NetworkX to check if we get a hypercube for large N
        # this is INREDIBLY slow, so use it once

        N = 10
        nx_hypercube = nx.hypercube_graph(N)
        H = hypergraphs.H_hypercube(N, 1)
        self.assertTrue(nx.is_isomorphic(nx_hypercube, nx.from_numpy_array(H)))



        # test p = 0: in this limit, each node is disconnected
        p = 0

        for j in range(len(Nlist)):
            N = Nlist[j]
            NH = HS_DIM_LIST[j]
            H = hypergraphs.H_hypercube(N, p)

            # ensure that the Hamiltonian matrix has dim 2**N by 2**N
            np.testing.assert_equal(H.shape, (NH, NH))

            # ensure that everything is just zeros!
            np.testing.assert_equal(H, 0)

        
        # a couple of tests for p = 0.5
        p = 0.5
        N = 10

        NH = 2**N
        H = hypergraphs.H_hypercube(N, p)

        # ensure that the Hamiltonian matrix has dim 2**N by 2**N
        np.testing.assert_equal(H.shape, (NH, NH))

        # the number of edges should be between the two limits
        # (almost certainly for N = 10)
        self.assertTrue(0 < np.sum(H) < N*NH)

    
    def test_H_PXP(self):
        
        """ Test the Fibonacci cube against analytic results from the 
        mathematics literature. """


        # source: http://fare.tunes.org/files/fun/fibonacci.lisp
        fib = lambda n:pow(2<<n,n+1,(4<<2*n)-(2<<n)-1)%(2<<n)

        def number_edges(n):
            """ The number of edges a Fibonacci cube of dim N has.
            see: https://www.sciencedirect.com/science/article/pii/S0012365X05002748"""

            fiblist1 = np.array([fib(i) for i in range(1, n-1)])
            fiblist2 = np.array([fib(n+1-i) for i in range(1, n-1)])
            num_e = fib(n+1) + np.sum(fiblist1 * fiblist2)

            return num_e

        Nlist = range(1, 11)
        HS_DIM_LIST = np.array([fib(N+2) for N in Nlist])

        # first test p = 1. The Hamiltonian should represent a complete PXP graph
        p = 1
        for Ni, N in enumerate(Nlist):
            NH = HS_DIM_LIST[Ni]

            H = hypergraphs.H_PXP(N, p)

            # ensure that the Hamiltonian matrix has correct dimensions
            np.testing.assert_equal(H.shape, (NH, NH))

            # ensure that the total number of edges is correct
            self.assertEqual(np.sum(H)/2, number_edges(N))

            # ensure that the matrix is Hermitian
            np.testing.assert_array_equal(H, H.T)




            # for i in range(len(connect)):
            #     degree = connect[i]
            #     number = count[i]

            #     self.assertEqual(number_nodes_coord_k(N, degree), number)














if __name__ == "__main__":
    unittest.main()

