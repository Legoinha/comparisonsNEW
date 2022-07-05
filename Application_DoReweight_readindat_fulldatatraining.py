#!/usr/bin/env python
# coding: utf-8

# # Demonstration of distribution reweighting
# 
#     requirements:
#     xgboost
#     numpy
#     matplotlib
#     sklearn
#     pandas
#     ROOT 6.12
#     root_numpy (using python 2.7, download root_numpy from official web and install by hand)

# In[1]:


from __future__ import division
#get_ipython().run_line_magic('pylab', 'inline')
#figsize(16, 8)
import xgboost as xgb 
from xgboost import plot_importance
from xgboost import plot_tree
import numpy as np
import matplotlib.pyplot as plt 
from sklearn import metrics
import pandas as pd
import root_numpy
from array import array
import ROOT as r
from ROOT import TCut
import time as timer
import sys

r.gROOT.SetBatch(True)

time_start=timer.time()

#number = 100000
q2Bin = int(sys.argv[1])
parity = int(sys.argv[2])
year = int(sys.argv[3])
# In[2]:


variables = {} #dictionary containing all the variables

class variable:

    def __init__(self, name, title, binning=None):
        if binning is not None and (not isinstance(binning, list) or len(binning)!=3):
            raise Exception("Error in declaration of variable {0}. Binning must be a list like [nbins, min, max]".format(name))
        self.title = title
        self.binning = binning
        variables[name] = self #add yourself to dictionary

    def get_nbins(self):
        if self.binning is not None:
            return self.binning[0]
        else: return 50

    def get_xmin(self):
        if self.binning is not None:
            return self.binning[1]
        else: return 0

    def get_xmax(self):
        if self.binning is not None:
            return self.binning[2]
        else: return 0

#declare variables

#variable of BDT
variable("kstarmass", "kstar tagged mass", [100, 0.74, 1.05 ])
variable("bVtxCL","B decay vertex CL",[100, 0.0, 1.0])
variable("bLBS"," bLBS ",[100, 0, 1.5])
variable("bLBSE"," bLBSE ",[100, 0, 0.04]) #significance
variable("bLBSsig"," bLBSsig ",[100, 0, 200])

variable("bCosAlphaBS"," bCosAlphaBS ",[100, 0.994, 1.0])

variable("bDCABS","bDCABS",[100, -0.015, 0.015])
variable("bDCABSE","bDCABSE",[100, 0, 0.01]) #significance
variable("bDCABSsig","bDCABSsig",[100, -5, 5])

variable("kstTrk1DCABS","kstTrk1DCABS",[100, -1, 1])
variable("kstTrk1DCABSE"," kstTrk1DCABSE",[100, 0, 0.025]) #significance
variable("kstTrk1DCABSsig"," kstTrk1DCABSsig",[100, -50,50])

variable("kstTrk2DCABS","kstTrk2DCABS",[100, -1, 1])
variable("kstTrk2DCABSE"," kstTrk2DCABSE",[100, 0, 0.025]) #significance
variable("kstTrk2DCABSsig"," kstTrk2DCABSsig",[100, -50,50])

variable("kstTrk2MinIP2D","kstTrk2MinIP2D ",[100, 0, 1])
variable("kstTrk1MinIP2D","kstTrk1MinIP2D ",[100, 0, 1])

variable("sum_isopt_04","sum_isopt_04 ",[100, 0, 40])


#others draw picture from chuqiao

variable("mumPt", "mumPT", [100,0,40])
variable("mupPt", "mupPT", [100,0,40])

variable("mumEta", "mumEta", [100,-2.5,2.5])
variable("mupEta", "mupEta", [100,-2.5,2.5])

variable("mumPhi", "mumPhi", [100,-3.14159,3.14159])
variable("mupPhi", "mupPhi", [100,-3.14159,3.14159])

variable("mumuLBS", "mumuLBS", [100,0,6.4])
variable("mumuLBSE", "mumuLBSE", [100,0,0.12])

variable("mumuDCA", "mumuDCA", [100,0,1])

variable("mumuCosAlphaBS", "mumuCosAlphaBS", [100,0.93,1])

variable("kstTrk1Pt", "kstTrk1Pt", [100,0,20])
variable("kstTrk2Pt", "kstTrk2Pt", [100,0,20])

variable("kstTrk1Eta", "kstTrk1Eta", [100,-2.5,2.5])
variable("kstTrk2Eta", "kstTrk2Eta", [100,-2.5,2.5])

variable("kstTrk1Phi", "kstTrk1Phi", [100,-3.14159,3.14159])
variable("kstTrk2Phi", "kstTrk2Phi", [100,-3.14159,3.14159])

variable("bPt","bPt",[100,0,80])
variable("bEta", "bEta", [100,-2.5,2.5])
variable("bPhi", "bPhi", [100,-3.14159,3.14159])

variable("tagged_mass", "tagged_mass", [100,5.0,5.6])

variable("cos_theta_l", "cos_theta_l", [100,-1,1])
variable("cos_theta_k", "cos_theta_k", [100,-1,1])
variable("phi_kst_mumu", "phi_kst_mumu", [100,-3.14159,3.14159])


#trkdca significance
#bvtxcl, blbssig, bcosalphaBS, 

###variables in note

variable("kstVtxCL","kstVtxCL",[100,0,1])


variable("weight","PUweight",[100,-1,3])


# In[3]:


def get_needed_branches_from_var(tree, var): #tree must be a TTree, cut must be a string or TCut
    if isinstance(var, TCut): var = var.GetTitle()
    list_of_branches = [x.GetName() for x in tree.GetListOfBranches()]
    list_of_needed_branches = [x for x in list_of_branches if x in var]
    return list_of_needed_branches

def only_activate_variables(tree, variables):
    tree.SetBranchStatus("*", 0)
    for variable in variables:
        branches = get_needed_branches_from_var(tree, variable)
        for branch in branches:
            tree.SetBranchStatus(branch, 1)

def plot_var(varname, cut="1", reweight=False):
    #only_activate_variables(rdata, [varname, cut]) #this speeds things up
    #only_activate_variables(MC, ["weight", "tagged_mass","weight_reco_pteta", varname, cut]) #this speeds things up
    rdata.SetBranchStatus("nsig_sw") #this ensure the introduction of SWeight of real data for signal extraction from background
    if reweight: 
        MC.SetBranchStatus("MCw", 1) #this needs to be done manually because it's on the friend tree
    
    hdata = r.TH1F("hdata", "hdata", variables[varname].get_nbins(), variables[varname].get_xmin(), variables[varname].get_xmax())
    swname = "nsig_sw*"
    #rdata.Draw(varname+">>hdata", cut, "goff")
    #rdata.Draw(varname+">>hdata", swname+"{0}".format(cut), "goff")
    rdata.Draw(varname+">>hdata", swname+"{0}".format(cut), "goff")#, number, 0)
    if hdata.Integral()==0:
        print("Empty histogram!")
        return
    hdata.Scale(1./hdata.Integral())
    
    hMC = r.TH1F("hMC", "hMC", hdata.GetNbinsX(), hdata.GetXaxis().GetXmin(), hdata.GetXaxis().GetXmax())
    #wname = "MCw*" if reweight else ""
    #MC.Draw(varname+">>hMC", "PUweight*"+wname+"({0})".format(cut), "same goff")
    #MC.Draw(varname+">>hMC", wname+"({0})".format(cut), "same goff")
    wname = ""
    MC.Draw(varname+">>hMC", "weight*(tagged_mass>5.0&&tagged_mass<5.6)", "same goff")#, number, 0)
    if hMC.Integral()==0:
        print("Empty histogram!")
        return
    hMC.Scale(1./hMC.Integral())
    hMC.SetLineWidth(3)
    hMC.SetLineColor(r.kRed)
    hMC.SetFillColor(0)
    hMC.SetFillStyle(0)
    print("========================================================")
    print("Chi2 test of {} (MC vs data)".format(varname))
    print("<< {} >>".format(varname)); print(hMC.Chi2Test(hdata,"WWP")); print( "<< {} >>".format(varname))
    print("--------------------------------------------------------")
    
    hWMC = r.TH1F("hWMC", "hWMC", hdata.GetNbinsX(), hdata.GetXaxis().GetXmin(), hdata.GetXaxis().GetXmax())
    wname = "MCw*"
    MC.Draw(varname+">>hWMC", "weight*MCw*(tagged_mass>5.0&&tagged_mass<5.6)", "same goff")#, number, 0)
    if hWMC.Integral()==0:
        print("Empty histogram!")
        return
    hWMC.Scale(1./hWMC.Integral())
    hWMC.SetLineWidth(3)
    hWMC.SetLineColor(r.kRed)
    hWMC.SetFillColor(0)
    hWMC.SetFillStyle(0)    
    print("Chi2 test of weighted {} (MC vs data)".format(varname))
    print("<< {} >>".format(varname)); print(hWMC.Chi2Test(hdata,"WWP")); print( "<< {} >>".format(varname))
    print("========================================================\n\n")
    
    c = r.TCanvas("canvas","canvas",1600,600)
    c.Divide(2,1)
    
    ############################################################################
    c.cd(1)
    hs = r.THStack()
    hs.Add(hMC, "hist")
    hdata.SetMarkerStyle(9)
    hs.Add(hdata, "P")
    #Upper plot will be in pad1
    pad1 = r.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1.SetTopMargin(0.05)
    pad1.SetBottomMargin(0.05)
    pad1.SetGridx()         #Vertical grid
    pad1.Draw()             #Draw the upper pad: pad1
    pad1.cd()               #pad1 becomes the current pad
    hs.Draw("nostack")
    hs.SetTitle("")
    #Y axis h1 plot settings
    hs.GetYaxis().SetTitleSize(30)
    hs.GetYaxis().SetTitleFont(43)
    hs.GetYaxis().SetTitleOffset(0.5)
    hs.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hs.GetYaxis().SetLabelSize(25)
    hs.GetXaxis().SetTitleSize(0)
    hs.GetXaxis().SetLabelSize(20)
    
    c.cd(1)          #Go back to the main canvas before defining pad2
    pad2 = r.TPad("pad2", "pad2", 0, 0.0, 1, 0.25)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.5)
    pad2.SetGridx() #vertical grid
    pad2.Draw()
    pad2.cd()       #pad2 becomes the current pad
    
    #Define the ratio plot
    hratio = hdata.Clone("hratio")
    hratio.SetLineColor(r.kBlack)
    hratio.SetMinimum(0.35)  #Define Y
    hratio.SetMaximum(1.65) #range
    hratio.GetXaxis().SetTitle(variables[varname].title)
    #hratio.Sumw2()
    hratio.SetStats(0)      #No statistics on lower plot
    hratio.Divide(hMC)
    hratio.Draw("ep")
    hratio.SetTitle("")
    
    #Y axis ratio plot settings
    hratio.GetYaxis().SetTitle("ratio")
    hratio.GetYaxis().SetNdivisions(505)
    hratio.GetYaxis().SetTitleSize(30)
    hratio.GetYaxis().SetTitleFont(43)
    hratio.GetYaxis().SetTitleOffset(1.0)
    hratio.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hratio.GetYaxis().SetLabelSize(25)
    
    #X axis ratio plot settings
    hratio.GetXaxis().SetTitleSize(25)
    hratio.GetXaxis().SetTitleFont(43)
    hratio.GetXaxis().SetTitleOffset(4.8)
    hratio.GetXaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hratio.GetXaxis().SetLabelSize(20)
    hratio.GetXaxis().SetLabelOffset(0.01)
    
    line = r.TLine(hratio.GetXaxis().GetXmin(),1,hratio.GetXaxis().GetXmax(), 1)
    line.SetLineStyle(3)
    line.Draw()
    
    ############################################################################
    c.cd(2)
    whs = r.THStack()
    whs.Add(hWMC, "hist")
    hdata.SetMarkerStyle(9)
    whs.Add(hdata, "P")
    #Upper plot will be in pad1
    wpad1 = r.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    wpad1.SetTopMargin(0.05)
    wpad1.SetBottomMargin(0.05)
    wpad1.SetGridx()         #Vertical grid
    wpad1.Draw()             #Draw the upper pad: pad1
    wpad1.cd()               #pad1 becomes the current pad
    whs.Draw("nostack")
    whs.SetTitle("")
    #Y axis h1 plot settings
    whs.GetYaxis().SetTitleSize(30)
    whs.GetYaxis().SetTitleFont(43)
    whs.GetYaxis().SetTitleOffset(0.5)
    whs.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    whs.GetYaxis().SetLabelSize(25)
    whs.GetXaxis().SetTitleSize(0)
    whs.GetXaxis().SetLabelSize(20)
    
    c.cd(2)          #Go back to the main canvas before defining pad2
    wpad2 = r.TPad("pad2", "pad2", 0, 0.0, 1, 0.25)
    wpad2.SetTopMargin(0)
    wpad2.SetBottomMargin(0.5)
    wpad2.SetGridx() #vertical grid
    wpad2.Draw()
    wpad2.cd()       #pad2 becomes the current pad
    
    #Define the ratio plot
    whratio = hdata.Clone("hratio_W")
    whratio.SetLineColor(r.kBlack)
    whratio.SetMinimum(0.35)  #Define Y
    whratio.SetMaximum(1.65) #range
    whratio.GetXaxis().SetTitle(variables[varname].title)
    #whratio.Sumw2()
    whratio.SetStats(0)      #No statistics on lower plot
    whratio.Divide(hWMC)
    whratio.Draw("ep")
    whratio.SetTitle("")
    
    #Y axis ratio plot settings
    whratio.GetYaxis().SetTitle("ratio")
    whratio.GetYaxis().SetNdivisions(505)
    whratio.GetYaxis().SetTitleSize(30)
    whratio.GetYaxis().SetTitleFont(43)
    whratio.GetYaxis().SetTitleOffset(1.0)
    whratio.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    whratio.GetYaxis().SetLabelSize(25)
    
    #X axis ratio plot settings
    whratio.GetXaxis().SetTitleSize(25)
    whratio.GetXaxis().SetTitleFont(43)
    whratio.GetXaxis().SetTitleOffset(4.8)
    whratio.GetXaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    whratio.GetXaxis().SetLabelSize(20)
    whratio.GetXaxis().SetLabelOffset(0.01)
    
    wline = r.TLine(hratio.GetXaxis().GetXmin(),1,hratio.GetXaxis().GetXmax(), 1)
    wline.SetLineStyle(3)
    wline.Draw()

    ############################################################################
    #c.SaveAs("Psi2SK_plots_fulldataapplication/data_vs_mc/"+wname+varname+".png")
    #c.SaveAs("Psi2SK_plots_fulldataapplication/data_vs_mc/"+wname+varname+".pdf")
    c.SaveAs("plots/bin{}/{}_gitv2_addPtEtaPhi/".format(q2Bin,year)+varname+".png")
    
    #at the moment I'm saving the legend in a separate canvas EACH TIME. This is a quick workaraound and it is inefficient. It needs to be fixed later.
    cleg = r.TCanvas("legend","legend",500,350)
    leg = r.TLegend(0.1,0.1,0.9,0.9)
    leg.AddEntry(hdata, "data", "ep")
    leg.AddEntry(hMC, "simulation", "f")
    leg.Draw()
    cleg.SaveAs("plots/bin{}/{}_gitv2_addPtEtaPhi/".format(q2Bin,year)+"/legend"+".png")
    #cleg.SaveAs("Psi2SK_plots_fulldataapplication/data_vs_mc/legend"+".pdf")
    

    #================= Draw 3 histograms in one plot
    c_3In1 = r.TCanvas("canvas_3In1","canvas_3In1",800,600)
    hs_3In1 = r.THStack()
    hdata.SetMarkerStyle(9)
    hdata.SetMarkerColor(r.kBlack)
    hs_3In1.Add(hdata, "P")
    hMC.SetLineColor(r.kRed)
    hWMC.SetLineColor(r.kGreen)
    hs_3In1.Add(hMC, "hist")
    hs_3In1.Add(hWMC, "hist")
    hs_3In1.Add(hdata, "P")
    #pad1
    pad1_3In1 = r.TPad("pad1", "pad1", 0, 0.25, 1, 1.0)
    pad1_3In1.SetTopMargin(0.05)
    pad1_3In1.SetBottomMargin(0.05)
    pad1_3In1.SetGridx()         #Vertical grid
    pad1_3In1.Draw()             #Draw the upper pad: pad1
    pad1_3In1.cd()               #pad1 becomes the current pad
    hs_3In1.Draw("nostack")
    hs_3In1.SetTitle("")
    #Y axis h1 plot settings
    hs_3In1.GetYaxis().SetTitleSize(30)
    hs_3In1.GetYaxis().SetTitleFont(43)
    hs_3In1.GetYaxis().SetTitleOffset(0.5)
    hs_3In1.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hs_3In1.GetYaxis().SetLabelSize(25)
    hs_3In1.GetXaxis().SetTitleSize(0)
    hs_3In1.GetXaxis().SetLabelSize(20)
    #pad2
    c_3In1.cd()
    pad2_3In1 = r.TPad("pad2", "pad2", 0, 0.0, 1, 0.25)
    pad2_3In1.SetTopMargin(0)
    pad2_3In1.SetBottomMargin(0.5)
    pad2_3In1.SetGridx() #vertical grid
    pad2_3In1.Draw()
    pad2_3In1.cd()       #pad2 becomes the current pad
    hratio_3In1 = r.THStack()
    hratio_RDvsMC = hdata.Clone("hratio_RDvsMC")
    hratio_RDvsMC.SetLineColor(r.kRed)
    hratio_RDvsMC.SetMarkerColor(r.kRed)
    hratio_RDvsMC.SetMarkerStyle(20)
    hratio_RDvsMC.Divide(hMC)
    hratio_RDvsWMC = hdata.Clone("hratio_RDvsWMC")
    hratio_RDvsWMC.SetLineColor(r.kGreen)
    hratio_RDvsWMC.SetMarkerColor(r.kGreen)
    hratio_RDvsWMC.SetMarkerStyle(21)
    hratio_RDvsWMC.Divide(hWMC)
    hratio_RDvsMC.SetStats(0)      #No statistics on lower plot
    hratio_RDvsWMC.SetStats(0)      #No statistics on lower plot
    hratio_3In1.Add(hratio_RDvsMC, "P")
    hratio_3In1.Add(hratio_RDvsWMC, "P")
    hratio_3In1.Draw("nostack")
    hratio_3In1.SetMinimum(0.35)  #Define Y
    hratio_3In1.SetMaximum(1.65) #range
    hratio_3In1.GetXaxis().SetTitle(variables[varname].title)
    hratio_3In1.SetTitle("")
    #Y axis ratio plot settings
    hratio_3In1.GetYaxis().SetTitle("ratio")
    hratio_3In1.GetYaxis().SetNdivisions(505)
    hratio_3In1.GetYaxis().SetTitleSize(30)
    hratio_3In1.GetYaxis().SetTitleFont(43)
    hratio_3In1.GetYaxis().SetTitleOffset(1.0)
    hratio_3In1.GetYaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hratio_3In1.GetYaxis().SetLabelSize(25)
    #X axis ratio plot settings
    hratio_3In1.GetXaxis().SetTitleSize(25)
    hratio_3In1.GetXaxis().SetTitleFont(43)
    hratio_3In1.GetXaxis().SetTitleOffset(4.8)
    hratio_3In1.GetXaxis().SetLabelFont(43) #Absolute font size in pixel (precision 3)
    hratio_3In1.GetXaxis().SetLabelSize(20)
    hratio_3In1.GetXaxis().SetLabelOffset(0.01)
    line_3In1 = r.TLine(hratio_3In1.GetXaxis().GetXmin(),1,hratio_3In1.GetXaxis().GetXmax(), 1)
    line_3In1.SetLineStyle(3)
    line_3In1.Draw()
    c_3In1.SaveAs("plots/bin{}/{}_gitv2_addPtEtaPhi/".format(q2Bin,year)+varname+"_3In1.png")
    cleg_3In1 = r.TCanvas("legend_3In1","legend_3In1",500,350)
    leg_3In1 = r.TLegend(0.1,0.1,0.9,0.9)
    leg_3In1.AddEntry(hdata, "data", "ep")
    leg_3In1.AddEntry(hMC, "simulation", "lf")
    leg_3In1.AddEntry(hWMC, "RW simulation", "lf")
    leg_3In1.Draw()
    cleg_3In1.SaveAs("plots/bin{}/{}_gitv2_addPtEtaPhi/".format(q2Bin,year)+"legend"+"_3In1.png")
    cleg_ratio_3In1 = r.TCanvas("legend_3In1","legend_3In1",500,350)
    leg_ratio_3In1 = r.TLegend(0.1,0.1,0.9,0.9)
    leg_ratio_3In1.AddEntry(hratio_RDvsMC, "RD/MC", "ep")
    leg_ratio_3In1.AddEntry(hratio_RDvsWMC, "RD/RWMC", "ep")
    leg_ratio_3In1.Draw()
    cleg_ratio_3In1.SaveAs("plots/bin{}/{}_gitv2_addPtEtaPhi/".format(q2Bin,year)+"legend"+"_ratio_3In1.png")
    #=================
    #c_3In1.Delete()
    #hs_3In1.Delete()
    #pad1_3In1.Delete()
    #pad2_3In1.Delete()
    #hratio_3In1.Delete()
    #line_3In1.Delete()
    #cleg_3In1.Delete()
    #leg_3In1.Delete()
    #cleg_ratio_3In1.Delete()

    #cleanup
    hs.Delete()
    #hratio.Delete()
    hMC.Delete()
    hWMC.Delete()
    hdata.Delete()
    pad1.Delete()
    pad2.Delete()


# ## Prepare for data
#     

# In[4]:

print("Start------------------------------------------")

#This column is used to draw picture
columns_draw = [
                "bPt","bEta","tagged_mass","bVtxCL","bDCABS","bDCABSE","bDCABSsig","bLBS","bLBSE","bLBSsig",
                "kstVtxCL","kstarmass",
                "mumPt","mumEta","mumPhi",
                "mupPt","mupEta","mupPhi",
                #"kstTrk1Pt","kstTrk1Eta","kstTrk1Phi",
                #"kstTrk2Pt","kstTrk2Eta","kstTrk2Phi",
                "kstTrk2Pt","kstTrk2Eta","kstTrk2Phi",
                "kstTrk1Pt","kstTrk1Eta","kstTrk1Phi",
                "cos_theta_l","cos_theta_k","phi_kst_mumu",
                "kstTrk1DCABS","kstTrk1DCABSE","kstTrk1DCABSsig","kstTrk1MinIP2D",
                "kstTrk2DCABS","kstTrk2DCABSE","kstTrk2DCABSsig","kstTrk2MinIP2D",
                #"kstTrkmDCABSsig","kstTrkmMinIP2D",
                #"kstTrkpDCABSsig","kstTrkpMinIP2D",
                "sum_isopt_04"

                ]

# In[5]:

print("TTrees preparation------------------------------------------")

#cut_bpeak = "(Bmass>5.0 && Bmass<5.5)"
cut_bpeak = "(tagged_mass>5.0 && tagged_mass<5.6)"


rdata = r.TChain("ntuple")
rdata.Add('/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoDATADataset_b{}_{}_p{}_aftersel.root'.format(q2Bin,year,parity))

#rdata.Add('/eos/user/c/cjiang/selected/data/2018/Psi2SK/after_preselection_after_plainized/plainized_data_cutPsip0_all_2018_UL_MINIAODv1_all_aftercutPsip0.root')
#rdata.AddFriend("tree_sw","/afs/cern.ch/work/c/cjiang/selectUL/CMSSW_10_6_20/src/sel/data/2018/copytree_all_cutPsip0/splot_for_psi2sk_data/tree_sw.root")

MC = r.TChain("ntuple")
MC.Add('/afs/cern.ch/user/x/xuqin/cernbox/workdir/B0KstMuMu/reweight/Tree/final/gitv2/recoMCDataset_b{}_{}_p{}_aftersel.root'.format(q2Bin,year,parity))
#MC.Add('/eos/user/c/cjiang/selected/MC/2018/Psi2SK/after_preselection_after_plainized/plainized_PileupRweight_BuToMuMuK_SIM_2018_UL_MINIAODv1_all_aftercutPsip0.root')

MC_friend = r.TTree("wTree", "weights tree")
leafValues = array("f", [0.0])
weight_branch = MC_friend.Branch("MCw", leafValues,"MCw[1]/F")
print("Wtree builded")

print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

# In[6]:

print("Training samples preparation------------------------------------------")

#This column is used as input to the XGBoost, branchs in your tree

#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','Bcosalphabs','Mumdcabserr']
#columns = ['Trkdcabserr','Blsbs','Blsbserr','Mumdcabserr','Mupdcabserr','MumEta']
#columns = ['Trkdcabserr','Blsbs','Blsbserr','Mumdcabserr','Mupdcabserr','MumEta','Bcosalphabs']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','Bcosalphabs']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','Mumdcabserr','Mupdcabserr','Trkdcabserr','Blsbserr']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','MumPhi','Mumdcabserr','Mupdcabserr','Mumulsbserr']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','MumPhi','MupPhi','Mumdcabserr','Mupdcabserr','Mumulsbserr']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','MumEta','MumPhi','MupPhi','MumPt','MupPt']
#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','Bcosalphabs','MumEta','MumPhi','MumudR']
#columns = ['Trkdcabs','Trkdcabserr','Bvtxcl','Blsbs','Blsbserr','TrkPt','Bcosalphabs']
#columns = ['Trkdcabs','Bvtxcl','Blsbs','TrkPt','Bcosalphabs','MumEta','MumPhi','MumPt','MupPt','BPt','Mumdcabserr','Mupdcabserr','Trkdcabserr','Blsbserr']

#columns = ['Trkdcasigbs','Bvtxcl','Blxysig','TrkPt','Bcosalphabs','MumEta','MumPhi']
#columns = ['bVtxCL', 'bLBSsig', 'bCosAlphaBS', 'bDCABSsig', 'kstTrkmDCABSsig', 'kstTrkpDCABSsig','kstTrkmMinIP2D','kstTrkpMinIP2D','sum_isopt_04','bPt','kstTrkmPt','kstTrkpPt','mumPt','mupPt','bEta','kstTrkmEta','kstTrkpEta','mumEta','mupEta','bPhi','kstTrkmPhi','kstTrkpPhi','mumPhi','mupPhi']
#columns = ['bVtxCL', 'bLBSsig', 'bCosAlphaBS', 'bDCABSsig', 'kstTrk1DCABSsig', 'kstTrk2DCABSsig','kstTrk1MinIP2D','kstTrk2MinIP2D','sum_isopt_04','bPt','bEta','bPhi']#,'kstTrk1Pt','kstTrk2Pt','mumPt','mupPt']
columns = ['bVtxCL', 'bLBS', 'bLBSE' ,'bPt','bEta','bPhi','bCosAlphaBS', 'bDCABS','bDCABSE', 'kstTrk1DCABS','kstTrk1DCABSE','kstTrk2DCABS','kstTrk2DCABSE','kstTrk1MinIP2D','kstTrk2MinIP2D','sum_isopt_04']
sw_branch = ['nsig_sw']

#Read data using root_numpy, define the tree name, branch list

phsp_ori = root_numpy.tree2array(MC,branches=columns)#,start=0, stop=number)
print("MC sample readed------------------------------------------")

#Translate data into pandas formation
phsp_only_X=pd.DataFrame(phsp_ori)

#Make labels for each data, data marks as 1, MC marks as 0 
phsp_only_Y=np.zeros(len(phsp_only_X))

phsp=pd.DataFrame(phsp_ori)

# In[7]:

phsp_only_a=np.array(phsp)

# In[8]:
# ## Initial the xgboost

#Prepare the input data for XGBoost
xg_phsp_only = xgb.DMatrix(phsp_only_X, label=phsp_only_Y)


# In[9]:

print("Machine learning------------------------------------------")

Save_Dir = './fulldata_trained_model_{}_gitv2_onlyRECO_addPhi.json'.format(year)

# In[10]:

#trained_bst = xgb.Booster()
#trained_bst = trained_bst.load_model(Save_Dir)
trained_bst = xgb.Booster(model_file=Save_Dir)
pr_phsp=np.array(trained_bst.predict(xg_phsp_only,validate_features=False).reshape(phsp_only_Y.shape[0], 2))

#weight_test=pr_test[:,1]/pr_test[:,0]
weight_phsp=pr_phsp[:,1]/pr_phsp[:,0]
print("MC_weights: weight_phsp=pr_phsp")
print(weight_phsp)

# In[11]:

print("MC weights------------------------------------------")

MCwFile = r.TFile("./data_MC_weights_{}_b{}p{}_gitv2_onlyRECO_forplot_addPhi.root".format(year,q2Bin,parity),"RECREATE")

for val in weight_phsp:
    leafValues[0] = val
    MC_friend.Fill()

MCwFile.cd()
MC_friend.Write()

MC.AddFriend(MC_friend)
for v in variables.keys(): plot_var(v, cut_bpeak, True)

#MCwFile.Close()


# In[12]:


time_end=timer.time()
print('totally cost',time_end-time_start)

