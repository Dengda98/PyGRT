""" Compare two nc files """

import sys
import numpy as np
from scipy.io import netcdf_file

ncpath1 = sys.argv[1]
ncpath2 = sys.argv[2]

with netcdf_file(ncpath1, mmap=False) as f1, netcdf_file(ncpath2, mmap=False) as f2:
    # compare keys
    keys1 = set(f1.variables.keys())
    keys2 = set(f2.variables.keys())

    if keys1 != keys2:
        raise ValueError(f"Different keys in {ncpath1} and {ncpath2}")
    
    # compare data
    for k in keys1:
        arr1 = f1.variables[k][:]
        arr2 = f2.variables[k][:]

        if np.all(arr1==0.0) and np.all(arr2==0.0):
            continue

        tol = 0.01
        err = np.sum(np.abs(arr1 - arr2)) / np.mean(np.abs(arr1))
        if err > tol:
            raise ValueError(f"Error({err}) > {tol} from keys={k} in {ncpath1} and {ncpath2}")
