import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

data_folder = Path("C:/Users/matte/PycharmProjects/Fortran-Based-CFD-Solver/Fortran_Solver/")
file_to_open = data_folder / "residuals.csv"

columns_name = ['Iteration','Energy','Mass','MomentumX','MomentumY']
plt_axs = ['# Iterations','']
plt_title = 'Norm-2 Residulas'

def plot_res(file_to_open,  columns_name, plt_axs, plt_title):
    res = pd.read_csv(file_to_open, sep=" , ", header=None, engine='python')
    res.columns = columns_name

    x = res.Iteration

    plt.figure(1)
    for i in range(1, len(res.columns)):
        plt.plot(x, res[res.columns[i]])
        plt.xlabel(plt_axs[0])
        plt.ylabel(plt_axs[1])
        plt.yscale('linear')
        plt.title(plt_title)
        # plt.ylabel(res.columns[i])
        plt.grid(True, which='minor')  # minor doesn't work
        plt.legend(res.columns[1::])
    plt.show()

plot_res(file_to_open, columns_name, plt_axs, plt_title)


