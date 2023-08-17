import time_evolution as te
import numpy as np
import matplotlib.pyplot as plt

DATA_PATH = "/Users/alex/Code/Hypercubes/analysis/data/"

# set font size and enable LaTeX
plt.rcParams.update({'font.size': 18})
plt.rcParams.update({'text.usetex': True})

def main():
    plot_MHD()


def plot_MHD():
    fig, ax = plt.subplots(nrows=2, ncols=3, figsize=[12, 8], sharex=True, sharey=True)
    ax = ax.flatten()
    props = dict(boxstyle='round', facecolor='white', alpha=0)

    # System sizes to use
    Nlist = np.array([8, 10, 12])
    clist = (Nlist - Nlist[0])/(Nlist[-1] - Nlist[0])
    colours = plt.cm.viridis_r(clist)

    # Defines which values of p are used, as p = N_coeff / N
    N_coeff_list = [7, 3, 1.5, 1.1,  1.01, 1]

    # Times to evaluate MSD at
    tmax = 10**6
    NT = 100

    # Should we use a linear space or a log space for timepoints?
    log = True
    if log:
        tlist = np.logspace(0, np.log10(tmax), NT)
    else:
        tlist = np.linspace(0, tmax, NT)

    # Do the plots for different values of p

    for Nci, N_coeff in enumerate(N_coeff_list):
        axis = ax[Nci]
        axis.semilogx()
        # axis.semilogy()
        # axis.set_xlim(right=10**(3.5))
        
        # Add a label
        axis.text(0.5, 0.92, f"$p = {N_coeff}/N$", transform=axis.transAxes, fontsize=18,verticalalignment='top', bbox=props)

        for Ni, N in enumerate(Nlist):
            if N == 8: NR = 1000
            elif N == 10: NR = 500
            elif N == 12: NR = 100

            mean_hd = te.MHD(N, N_coeff_list[Nci], tmax, NT, NR, log)
            axis.plot(tlist, mean_hd, 'o-', c=colours[Ni], label=f"${N}$")

    # labels etc

    ax[0].set_yticks([2,3,4,6])
    ax[0].set_yticklabels(["$2$", " ", "$4$", "$6$"])

    ax[3].set_xlabel(r"$t$")
    ax[4].set_xlabel(r"$t$")
    ax[5].set_xlabel(r"$t$")
    ax[0].set_ylabel(r"$\mathrm{MHD}$")
    ax[3].set_ylabel(r"$\mathrm{MHD}$")

    ax[0].legend(loc=[0.6, 0.2], ncol=2, handlelength=0.2, columnspacing=0.5, handletextpad=0.2, framealpha=0)

    plt.subplots_adjust(right=0.99, top=0.99, left=0.05, bottom=0.08, wspace=0.05, hspace=0.05)

    plt.savefig("time_evolution.pdf")


if __name__ == "__main__":
    main()
