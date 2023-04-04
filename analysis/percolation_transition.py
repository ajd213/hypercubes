# here write code to produce collapsed plot demonstrating transiton at 1/N

import hypercubes
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

"""Plot the mean/maximum cluster size against p."""


def main():
    fig, ax = plt.subplots(nrows=2, ncols=2, figsize=[8, 8], sharey=True)
    ax = ax.flatten()

    p_min = 0
    p_max = 1
    N_p = 101
    NR = 100

    # which system sizes to use
    Nlist = np.array([8, 10, 12])

    # define colours for the plots
    Nlist_numbers = (Nlist - Nlist[0]) / (Nlist[-1] - Nlist[0])
    clist = plt.cm.viridis_r(Nlist_numbers)


    for Ni, N in enumerate(Nlist):
        S, max = s_with_p(N, NR, p_min, p_max, N_p)
        plist = np.linspace(p_min, p_max, N_p)

        ax[0].plot(plist, max/(2**N), '^-',c=clist[Ni], markersize='3', label=f"${N}$")
        ax[1].plot(plist, S/(2**N), 'o-',c=clist[Ni], markersize='3')
        # ax[1].plot(plist, Sp/(2**N), '--', c=clist[Ni])

        ax[2].plot(plist*N, max/(2**N), '^-',c=clist[Ni], markersize='3')
        ax[3].plot(plist*N, S/(2**N), 'o-',c=clist[Ni], markersize='3')
        # ax[3].plot(plist*N, Sp/(2**N), '--', c=clist[Ni])

    ax[2].set_xlim([0,5])
    ax[3].set_xlim([0,10])
    ax[2].axvline(1,linestyle='--', c='black')
    ax[3].axvline(1,linestyle='--', c='black')

    ax[0].set_ylabel("$\mathrm{max}(s)/N_\mathcal{H}$")
    ax[1].set_ylabel("$(S/S')/N_\mathcal{H}$")
    ax[2].set_ylabel("$\mathrm{max}(s)/N_\mathcal{H}$")
    ax[3].set_ylabel("$(S/S')/N_\mathcal{H}$")

    ax[0].set_xlabel("$p$")
    ax[1].set_xlabel("$p$")

    ax[2].set_xlabel("$pN$")
    ax[3].set_xlabel("$pN$")



    props = dict(boxstyle='round', facecolor='white', alpha=0)

    ax[1].text(0.1, 0.9, "$S$",  fontsize=18, verticalalignment='top', bbox=props)
    ax[1].text(0.63, 0.9, "$S'$",  fontsize=18, verticalalignment='top', bbox=props)

    ax[3].text(1.75, 0.9, "$S$",  fontsize=18, verticalalignment='top', bbox=props)
    ax[3].text(8.3, 0.9, "$S'$",  fontsize=18, verticalalignment='top', bbox=props)



    ax[0].legend(loc=[0.8, 0.4], handlelength=0.2, columnspacing=0.5, handletextpad=0.2, framealpha=0)


    plt.subplots_adjust(right=0.95, top=0.99, left=0.11, bottom=0.1, wspace=0.2, hspace=0.25)
    plt.savefig("hypercube_percolation.pdf")




def s_with_p(N, NR, p_min, p_max, N_p):
    plist = np.linspace(p_min, p_max, N_p)

    Slist = np.zeros(N_p)
    max_sizes = np.zeros(N_p)

    for pi, p in enumerate(plist):
        cs = get_clusters(N, int(NR), p, DATA_PATH)

        Slist[pi] = distributions.S(cs)
        max_sizes[pi] = np.max(cs)

    return Slist, max_sizes


def get_clusters(N, NR, p, data_path):
    name = f"clusters_N{N}_NR{NR}_p{p:.4f}.npy"

    try:
        clusters = np.load(data_path + name)
    except FileNotFoundError:
        clusters = hypercubes.clusters(N, NR, p)
        np.save(data_path + name, clusters)
    return clusters


if __name__ == "__main__":
    main()

