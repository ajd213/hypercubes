""" An example analysis script. Here, we locate the percolation transition at p_c = 1/N using a scale collapse.

The essence of this plot is that if the percolation transition occurs at 1/N in the thermodynamic limit,
then if we plot a measure of percolation (such as the mean cluster size, S) against pN, then the data
should collapse around a common transition point of pN = 1.

"""

import hypergraphs
import numpy as np
import matplotlib.pyplot as plt
import distributions
import os

# set font size and enable LaTeX
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})

DATA_PATH = "./plotting_data/"
if not os.path.exists(DATA_PATH):
    os.makedirs(DATA_PATH)


def main():
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 8], sharey=True)
    ax = ax.flatten()

    # grid of percolation concentrations, between 0 (non-percolating) and 1 (percolating).
    p_min = 0
    p_max = 1
    N_p = 101
    plist = np.linspace(p_min, p_max, N_p)

    # number of realisations to use
    NR = 100000

    # system sizes to use
    Nlist = np.arange(8, 16, 2)

    # define colours for the plots
    Nlist_numbers = (Nlist - Nlist[0]) / (Nlist[-1] - Nlist[0])
    clist = plt.cm.viridis_r(Nlist_numbers)


    #####################
    ### PLOT THE DATA ###
    #####################

    for Ni, N in enumerate(Nlist):

        NH = 2**N # total number of nodes in the graph

        S, max = s_with_p(N, NR, plist)

        ax[0].plot(plist, max/NH, '^-', c=clist[Ni], markersize='3', label=f"${N}$")
        ax[1].plot(plist, S/NH, 'o-', c=clist[Ni], markersize='3')
        ax[2].plot(plist*N, max/NH, '^-', c=clist[Ni], markersize='3')
        ax[3].plot(plist*N, S/NH, 'o-', c=clist[Ni], markersize='3')


    #####################
    ###  FORMAT PLOT  ###
    #####################

    ax[2].set_xlim([0,5])
    ax[3].set_xlim([0,10])

    # the locations of the transitions 
    ax[2].axvline(1,linestyle='--', c='black')
    ax[3].axvline(1,linestyle='--', c='black')

    ax[0].set_ylabel("$\mathrm{max}(s)/N_\mathcal{H}$")
    ax[1].set_ylabel("$S/N_\mathcal{H}$")
    ax[2].set_ylabel("$\mathrm{max}(s)/N_\mathcal{H}$")
    ax[3].set_ylabel("$S/N_\mathcal{H}$")

    ax[0].set_xlabel("$p$")
    ax[1].set_xlabel("$p$")
    ax[2].set_xlabel("$pN$")
    ax[3].set_xlabel("$pN$")

    ax[0].legend(loc=[0.8, 0.4], handlelength=0.2, columnspacing=0.5, handletextpad=0.2, framealpha=0)

    plt.subplots_adjust(right=0.95, top=0.99, left=0.11, bottom=0.1, wspace=0.2, hspace=0.25)
    plt.savefig("hypercube_percolation.pdf")



def s_with_p(N, NR, plist):
    """ Create lists of the mean cluster size (Slist) and maximum cluster size."""
    N_p = len(plist)

    Slist = np.zeros(N_p)
    max_sizes = np.zeros(N_p)

    for pi, p in enumerate(plist):
        cs = distributions.get_clusters_hypercube(N, NR, p, DATA_PATH)
        Slist[pi] = distributions.S(cs)
        max_sizes[pi] = np.max(cs)

    return Slist, max_sizes


if __name__ == "__main__":
    main()

