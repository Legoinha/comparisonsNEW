from ROOT import TTree,TFile,TChain
import root_pandas
import numpy as np
import pandas as pd
import root_numpy
from array import array
pd.options.mode.chained_assignment = None  # default='warn'
import sys

q2Bin = int(sys.argv[1])
parity = int(sys.argv[2])
year = int(sys.argv[3])

rMC = TChain("ntuple")
rMC.Add("/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}.root".format(q2Bin,year,parity))
rw = TChain("wTree")
rw.Add("./data_MC_weights_{}_b{}p{}_gitv2_onlyRECO_weight.root".format(year,q2Bin,parity))


ds = pd.DataFrame(
    root_numpy.tree2array(
        tree = rMC
    )
)
dw = root_pandas.read_root("./data_MC_weights_{}_b{}p{}_gitv2_onlyRECO_weight.root".format(year,q2Bin,parity))
dw['MCw']=dw['MCw'].astype(float)
data_final = pd.merge(ds,dw,left_index=True,right_index=True)
ofile = "/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}_rw.root".format(q2Bin,year,parity)
data_final.to_root(ofile, key='ntuple', store_index=False)