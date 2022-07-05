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
data = int(sys.argv[4])

mc_sigma = 0.040
mc_mass  = 5.27783 
JPsiMass_ = 3.096916
nSigma_psiRej = 3.
selData = '( ( (bMass*tagB0 + (1-tagB0)*bBarMass)   > {M}-3*{S}   && \
            (bMass*tagB0 + (1-tagB0)*bBarMass)   < {M}+3*{S} ) && \
            ( pass_preselection ==1 ) && \
            (abs(mumuMass - {JPSIM}) < {CUT}*mumuMassE))'\
            .format(M=mc_mass,S=mc_sigma,  JPSIM=JPsiMass_, CUT=nSigma_psiRej)
selMC = selData + ' && (trig==1)'

print selData
print selMC

r = TChain("ntuple")

if (data==1):
    r.Add("/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoDATADataset_b{}_{}_p{}.root".format(q2Bin,year,parity))
else:
    r.Add("/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}.root".format(q2Bin,year,parity))

print ("Events before sel:" , r.GetEntries())
if (data==0):
    ds = pd.DataFrame(
        root_numpy.tree2array(
            tree = r,
            selection=selMC
        )
    )

else:
    ds = pd.DataFrame(
        root_numpy.tree2array(
            tree = r,
            selection=selData
        )
    )

if (data==1):
    ofile = "/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoDATADataset_b{}_{}_p{}_aftersel.root".format(q2Bin,year,parity)
    print ('\t...done. n events: ', len(ds))
    ds.to_root(ofile, key='ntuple', store_index=False)
else:
    ofile = "/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}_aftersel.root".format(q2Bin,year,parity)
    print ('\t...done. n events: ', len(ds))
    ds.to_root(ofile, key='ntuple', store_index=False)
