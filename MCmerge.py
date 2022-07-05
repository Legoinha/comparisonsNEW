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

if (q2Bin==4):
    rMC.Add("/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_*_100.root".format(q2Bin,year))

elif q2Bin==6:
    rMC.Add("/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv1/recoMCDataset_b{}_{}_1_1_49_57.root".format(q2Bin,year))

print ("MC events:" , rMC.GetEntries())
ds = pd.DataFrame(
    root_numpy.tree2array(
        tree = rMC,
        selection="eventN%2=={}".format(parity)
    )
)

ofile = "/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}.root".format(q2Bin,year,parity)
print ('\t...done. n events: ', len(ds))
ds.to_root(ofile, key='ntuple', store_index=False)