#! /bin/python

import numpy as np
from glob import glob

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
plt.rcParams.update({'font.size': 14})





def save_image(filename, pad_inches = 0, tight = False):
	# PdfPages is a wrapper around pdf
	# file so there is no clash and create
	# files with no error.
    with PdfPages(filename) as p:

        # get_fignums Return list of existing
        # figure numbers
        fig_nums = plt.get_fignums()
        figs = [plt.figure(n) for n in fig_nums]

        # iterating over the numbers in list
        for fig in figs:

            # and saving the files
            if tight:
                fig.savefig(p, format='pdf', bbox_inches='tight', pad_inches = pad_inches)
            else:
                fig.savefig(p, format='pdf', pad_inches = pad_inches)

            plt.close(fig)


if __name__ == "__main__":
    files = glob("*.dat")

    for file in files:
        
        T, mu, n, rho, p = np.loadtxt(file).T
        
        fig, ax = plt.subplots(3, 1, sharex = True, figsize = (7, 8))

        fig.suptitle("${}$Fermions macroparameters")

        ax[0].plot(T, n)
        ax[1].plot(T, rho)
        ax[2].plot(T, p)

        ax[0].set_ylabel(r"$n$, cm$^{-3}$")
        ax[1].set_ylabel(r"$\rho$, erg$\cdot$cm$^{-3}$")
        ax[2].set_ylabel(r"$p$, dyn$\cdot$cm$^{-2}$")

        ax[2].set_xlabel(r"T, K")

        save_image(f"{file[:-4]}.pdf")

