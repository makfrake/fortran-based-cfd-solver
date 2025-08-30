import numpy as np
import matplotlib as plt
import pandas as pd
from pathlib import Path

data_folder = Path("C:/Users/matte/PycharmProjects/Fortran-Based-CFD-Solver/Fortran_Solver/")
file_to_open = data_folder / "residuals.csv"

x = pd.read_csv(file_to_open, sep=" , ", header=None, engine='python')
print(x[0])

