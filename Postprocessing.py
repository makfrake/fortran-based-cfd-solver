import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

data_folder = Path("C:/Users/matte/PycharmProjects/Fortran-Based-CFD-Solver/Fortran_Solver/")
file_to_open = data_folder / "residuals.csv"

res = pd.read_csv(file_to_open, sep=" , ", header=None, engine='python')
res.columns = ['Energy','Mass','MomentumX','MomentumY']
print(res['Energy'])

x = np.linspace(0,len(res)-1,len(res))

plt.figure(1)
for i in range(0,len(res.columns)-1):
    plt.plot(x, res[res.columns[i]])
    plt.xlabel('x')
    plt.show()




