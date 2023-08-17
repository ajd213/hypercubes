import os
import sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import distributions
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = "/Users/alex/Code/Hypercubes/analysis/data/"

# set font size and enable LaTeX
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})

def main():
    plot_path_lengths()
    transition_path_lengths()


def transition_path_lengths():

    fig, ax = plt.subplots(nrows=1, ncols=2, figsize=[8, 4], sharey=True)
    ax = ax.flatten()
    props = dict(boxstyle='round', facecolor='white', alpha=0)

    # System sizes to use
    Nlist = np.array([10, 12, 14, 16, 18])#, 20, 22, 24])
    clist = (Nlist - Nlist[0])/(Nlist[-1] - Nlist[0])
    colours = plt.cm.viridis_r(clist)

    ax[0].set_ylim(top=10**(-0.8), bottom = 10**(-6))


    # Defines which values of p are used, as p = N_coeff / N
    N_coeff_list = [1.05, 1]
    NRlist = [100000, 100000]

    ###################################
    ########## PLOTS FOR p<1 ##########
    ###################################

    for Nci, N_coeff in enumerate(N_coeff_list):
        axis = ax[Nci]
        axis.semilogy()

        # Add a label
        axis.text(0.65, 0.95, f"$p = {N_coeff}/N$", transform=axis.transAxes, fontsize=18,verticalalignment='top', bbox=props)

        NR = NRlist[Nci]

        for Ni, N in enumerate(Nlist):
            p = N_coeff_list[Nci]/N
            path_lengths = np.concatenate(distributions.get_path_lengths_hypercube_LC(N, NR, p, DATA_PATH))

            # bins at every integer
            bins = np.arange(0, np.max(path_lengths) + 2) -0.5

            axis.hist(path_lengths, bins=bins, histtype='step', density=True, color=colours[Ni], label=f"${N}$")

    ax[0].set_ylabel(r"$P_\ell(\ell)$")
    ax[0].set_xlabel(r"$\ell$")
    ax[1].set_xlabel(r"$\ell$")

    ax[0].legend(loc=[0.65, 0.45], ncol=2, handlelength=0.2, columnspacing=0.5, handletextpad=0.2, framealpha=0)

    plt.subplots_adjust(right=0.99, top=0.99, left=0.12, bottom=0.15, wspace=0.05, hspace=0.1)
    plt.savefig("path_lengths_transition.pdf")

def plot_path_lengths():

    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=[12, 8], sharey=True)
    ax = ax.flatten()
    props = dict(boxstyle='round', facecolor='white', alpha=0)

    # System sizes to use
    Nlist = np.array([8, 10, 12, 14])
    clist = (Nlist - Nlist[0])/(Nlist[-1] - Nlist[0])
    colours = plt.cm.viridis_r(clist)


    # Defines which values of p are used, as p = N_coeff / N
    N_coeff_list = [7, 3, 1.5, 1.1,  1.01, 1]
    NRlist = [2000, 2000, 2000, 4000, 10000, 10000]

    ###################################
    ########## PLOTS FOR p<1 ##########
    ###################################

    for Nci, N_coeff in enumerate(N_coeff_list):
        axis = ax[Nci]
        axis.semilogy()

        # Add a label
        axis.text(0.65, 0.95, f"$p = {N_coeff}/N$", transform=axis.transAxes, fontsize=18,verticalalignment='top', bbox=props)

        NR = NRlist[Nci]

        for Ni, N in enumerate(Nlist):
            p = N_coeff_list[Nci]/N
            
            path_lengths = np.concatenate(distributions.get_path_lengths_hypercube_LC(N, NR, p, DATA_PATH))

            # bins at every integer
            bins = np.arange(0, np.max(path_lengths) + 2) -0.5

            axis.hist(path_lengths, bins=bins, histtype='step', density=True, color=colours[Ni], label=f"${N}$")

    # TODO: add labels etc

    ax[3].set_xlabel(r"$\ell$")
    ax[4].set_xlabel(r"$\ell$")
    ax[5].set_xlabel(r"$\ell$")
    ax[0].set_ylabel(r"$P_\ell(\ell)$")
    ax[0].set_ylabel(r"$P_\ell(\ell)$")
    ax[3].set_ylabel(r"$P_\ell(\ell)$")
    

    ax[3].legend(loc=[0.65, 0.5], ncol=2, handlelength=0.2, columnspacing=0.5, handletextpad=0.2, framealpha=0)

    plt.subplots_adjust(right=0.99, top=0.99, left=0.08, bottom=0.08, wspace=0.05, hspace=0.1)

    plt.savefig("path_lengths.pdf")



if __name__ == "__main__":
    main()
