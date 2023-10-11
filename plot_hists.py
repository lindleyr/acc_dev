from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH1F
import math
from array import array

data = pd.read_csv("outputarr_withOuterLayers.txt",sep=' ',header=None)
data.apply(pd.to_numeric)

print(data[11])
print(type(data[11]))

layer = pd.Series(data[2]).values
r = pd.Series(data[3]).values
x = pd.Series(data[4]).values
y = pd.Series(data[5]).values
z = pd.Series(data[6]).values
roadeff = pd.Series(data[11]).values
roadtot = pd.Series(data[12]).values

#truth info
barcode = pd.Series(data[7]).values
charge = pd.Series(data[8]).values
pT = pd.Series(data[9]).values
d0 = pd.Series(data[10]).values
print(roadeff)

arrbins_lowpT = array('d',(0,5,10,15,20,30,40,50,100))
arrbins = array('d',(0,10,20,30,50,100,200,400,800))
eff_avg_pT = array('d',(0,0,0,0,0,0,0,0))
eff_norm_pT = array('d',(0,0,0,0,0,0,0,0))
nroadevts_pT = array('d',(0,0,0,0,0,0,0,0))
ntotevts_pT = array('d',(0,0,0,0,0,0,0,0))
arrbins_d0 = array('d',(range(-100,100,10)))
eff_avg_d0 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_avg_d0_5 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_avg_d0_10 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_avg_d0_15 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_norm_d0 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_norm_d0_5 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_norm_d0_10 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_norm_d0_15 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
print(arrbins_d0)
print(eff_avg_d0)

eff_all = ROOT.TH1F("eff_all","Road Efficiency Overall;Eff;Events",50,0.,1.2)
nroads = ROOT.TH1F("nroads","Roads Per Event;nRoads;Events",800,0,800)
eff_pT = ROOT.TH1F("eff_pT","Road Efficiency vs. pT;pT;Eff",len(arrbins)-1,arrbins)
eff_crude_pT = ROOT.TH1F("eff_crude_pT","Road Efficiency vs. pT;pT;Eff",len(arrbins_lowpT)-1,arrbins_lowpT)
eff_d0 = ROOT.TH1F("eff_d0","Road Efficiency vs. d0;d0;Eff",len(arrbins_d0)-1,arrbins_d0)
eff_d0_5 = ROOT.TH1F("eff_d0_5","Road Efficiency vs. d0;d0;Eff",len(arrbins_d0)-1,arrbins_d0)
eff_d0_10 = ROOT.TH1F("eff_d0_10","Road Efficiency vs. d0;d0;Eff",len(arrbins_d0)-1,arrbins_d0)
eff_d0_15 = ROOT.TH1F("eff_d0_15","Road Efficiency vs. d0;d0;Eff",len(arrbins_d0)-1,arrbins_d0)

nroadevts = 0

for i in range(len(layer)):
    eff_all.Fill(roadeff[i])
    nroads.Fill(roadtot[i])
    if roadtot[i] > 0:
        nroadevts += 1

print(nroadevts)
print(len(roadtot))
nroadevts = nroadevts/len(roadtot)
print(nroadevts)
#ugh i'm sure there's a better way to do this
for i in range(len(layer)):
    p = pT[i]/1000.
    for j in range(len(eff_avg_pT)):
        if p >= arrbins[j] and p < arrbins[j+1]:
            eff_norm_pT[j]+=1
            eff_avg_pT[j]+=roadeff[i]
        if p >= arrbins_lowpT[j] and p < arrbins_lowpT[j+1]:
            ntotevts_pT[j]+=1
            if roadtot[i] > 0:
                nroadevts_pT[j]+=1
    for k in range(len(eff_avg_d0)):
        if d0[i] >= arrbins_d0[k] and d0[i] < arrbins_d0[k+1]:
            eff_norm_d0[k]+=1
            eff_avg_d0[k]+=roadeff[i]
            if p > 5:
                eff_norm_d0_5[k]+=1
                eff_avg_d0_5[k]+=roadeff[i]
            if p > 10:
                eff_norm_d0_10[k]+=1
                eff_avg_d0_10[k]+=roadeff[i]
            if p > 15:
                eff_norm_d0_15[k]+=1
                eff_avg_d0_15[k]+=roadeff[i]

print(eff_avg_d0_5[0])

for i in range(len(eff_avg_pT)):
    eff_avg_pT[i] = eff_avg_pT[i]/eff_norm_pT[i]
    if ntotevts_pT[i] > 0:
        nroadevts_pT[i] = nroadevts_pT[i]/ntotevts_pT[i]
    else:
        nroadevts_pT[i] = 0
    eff_pT.SetBinContent(i+1,eff_avg_pT[i])
    eff_crude_pT.SetBinContent(i+1,nroadevts_pT[i])

print(nroadevts_pT)

for i in range(len(eff_avg_d0)):
    eff_avg_d0[i] = eff_avg_d0[i]/eff_norm_d0[i]
    eff_d0.SetBinContent(i+1,eff_avg_d0[i])
    eff_avg_d0_5[i] = eff_avg_d0_5[i]/eff_norm_d0_5[i]
    eff_avg_d0_10[i] = eff_avg_d0_10[i]/eff_norm_d0_10[i]
    eff_avg_d0_15[i] = eff_avg_d0_15[i]/eff_norm_d0_15[i]
    eff_d0_5.SetBinContent(i+1,eff_avg_d0_5[i])
    eff_d0_10.SetBinContent(i+1,eff_avg_d0_10[i])
    eff_d0_15.SetBinContent(i+1,eff_avg_d0_15[i])

c= ROOT.TCanvas("Canvas","Canvas",800,800)
eff_all.Draw()

c2 = ROOT.TCanvas("Canvas2","Canvas2",800,800)
nroads.Draw()

c3 = ROOT.TCanvas("Canvas3","Canvas3",800,800)
eff_pT.SetMinimum(0)
eff_pT.SetMaximum(1)
eff_pT.Draw()

c4 = ROOT.TCanvas("Canvas4","Canvas4",800,800)
eff_d0.SetMinimum(0)
eff_d0.SetMinimum(1)
eff_d0.Draw()

c5 = ROOT.TCanvas("Canvas5","Canvas5",800,800)
eff_d0_5.SetMinimum(0)
eff_d0_5.SetMaximum(1)
eff_d0_5.Draw()
eff_d0_10.SetLineColor(6)
eff_d0_10.SetLineColor(2)
eff_d0_10.Draw("Same")
eff_d0_15.Draw("Same")
leg = ROOT.TLegend(.1,.7,.3,.9,"Road Efficiency vs. d0")
leg.AddEntry(eff_d0_5,"pT > 5 GeV")
leg.AddEntry(eff_d0_10, "pT > 10 GeV")
leg.AddEntry(eff_d0_15, "pT > 15 GeV")
leg.Draw("Same")

c6 = ROOT.TCanvas("Canvas6","Canvas6",800,800)
eff_crude_pT.SetMinimum(0.9)
eff_crude_pT.SetMaximum(1.01)
eff_crude_pT.Draw()

out_file = ROOT.TFile("road_efficiency_all.root","RECREATE")
eff_all.Write()
nroads.Write()
eff_pT.Write()
eff_d0.Write()
c5.Write()
eff_crude_pT.Write()
