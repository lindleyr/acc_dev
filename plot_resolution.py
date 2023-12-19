from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH1F
import math
from array import array

data = pd.read_csv("outputfiles/outputroad_phi0305_eta0103_nopileup.txt",sep=' ',header=None)
data.apply(pd.to_numeric)

event = pd.Series(data[0]).values
recod0 = pd.Series(data[1]).values
truthd0 = pd.Series(data[2]).values
resod0 = pd.Series(data[3]).values
recoqpT = pd.Series(data[4]).values
truthqpT = pd.Series(data[5]).values
resoqpT = pd.Series(data[6]).values
truthphi = pd.Series(data[7]).values

resolution_pT = ROOT.TH1F("resolution_pT","Resolution q/pT;(Truth q/pT - Reco q/pT)/Truth q/pT;Events",25,0.,2.5)

resolution_d0 = ROOT.TH1F("resolution_d0","Resolution d0;(Truth d0 - Reco d0)/Truth d0;Events",25,0.,2.5)

ngoodbarcode = 0
efficiency = 0

for i in range(len(event)):
    resolution_pT.Fill(resoqpT[i]/truthqpT[i])
    resolution_d0.Fill(resod0[i]/truthd0[i])
    if(resod0[i]/truthd0[i] > 0.9 and resod0[i]/truthd0[i] < 1.1):
        print(event[i])

c_pT= ROOT.TCanvas("Canvas_reso_pT","Canvas_reso_pT",800,800)
resolution_pT.Draw()

c_d0 = ROOT.TCanvas("Canvas_reso_d0","Canvas_reso_d0",800,800)
resolution_d0.Draw()

out_file = ROOT.TFile("outputfiles/road_resolution_nopileup_phi0305_eta0103.root","RECREATE")
resolution_pT.Write()
resolution_d0.Write()
