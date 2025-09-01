from postprocessing import plot_res
from matplotlib import pyplot as plt
import pandas as pd
from pathlib import Path


data_folder = Path("C:/Users/matte/PycharmProjects/Fortran-Based-CFD-Solver/Fortran_Solver/")
file_to_open = data_folder / "residuals.csv"

columns_name = ['Iteration','Energy','Mass','MomentumX','MomentumY']
plt_axs = ['# Iterations','']
plt_title = 'Norm-2 Residulas'

plot_res(file_to_open, columns_name, plt_axs, plt_title)
