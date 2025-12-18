""" Compare sac files """

import sys
import numpy as np
from obspy import read

sacpath1 = sys.argv[1]
sacpath2 = sys.argv[2]

st1 = read(sacpath1)
st2 = read(sacpath2)

tol = 0.01
errLst = []

for tr1 in st1:
    ch = tr1.stats.channel
    tr2 = st2.select(channel=ch)[0]

    if tr1.stats.delta != tr2.stats.delta:
        raise ValueError(f"Different delta in {tr1} and {tr2} from {sacpath1} and {sacpath2}!")
    if tr1.stats.npts != tr2.stats.npts:
        raise ValueError(f"Different npts in {tr1} and {tr2} from {sacpath1} and {sacpath2}!")
    
    data1 = tr1.data.copy()
    data2 = tr2.data.copy()

    # All zeros
    if np.all(data1==0.0) and np.all(data2==0.0):
        errLst.append(0.0)
        continue
    
    err = np.sum(np.abs(data2 - data1))/np.mean(np.abs(data1))
    errLst.append(err)
    if np.any(err > tol):
        print(f"Error({err}) > {tol} in {tr1} and {tr2} from {sacpath1} and {sacpath2}!!")

errLst = np.array(errLst)

print('-'*100)
print(sacpath1, sacpath2)
print(errLst)
print('-'*100)
if np.any(errLst > tol):
    raise ValueError(f"Error!!")
