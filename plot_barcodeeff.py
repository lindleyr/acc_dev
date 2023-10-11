from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH1F
import math
from array import array

data = pd.read_csv("outputarr_nopileup_2around_avg.txt",sep=' ',header=None)
data.apply(pd.to_numeric)

print(data[11])
print(type(data[11]))

event = pd.Series(data[0]).values
layer = pd.Series(data[2]).values
r = pd.Series(data[3]).values
x = pd.Series(data[4]).values
y = pd.Series(data[5]).values
z = pd.Series(data[6]).values

roadbfrac = pd.Series(data[12]).values
roadtot = pd.Series(data[13]).values
eventtype = pd.Series(data[14]).values
#print(roadbfrac)

#truth info
hbarcode = pd.Series(data[7]).values
barcode = pd.Series(data[8]).values
charge = pd.Series(data[9]).values
pT = pd.Series(data[10]).values
d0 = pd.Series(data[11]).values

arrbins_lowpT = array('d',(0,5,10,15,20,30,40,50,100))
arrbins = array('d',(0,10,20,30,50,100,200,400,800))
eff_barcode_pT = array('d',(0,0,0,0,0,0,0,0))
eff_norm_pT = array('d',(0,0,0,0,0,0,0,0))
nroadevts_pT = array('d',(0,0,0,0,0,0,0,0))
ntotevts_pT = array('d',(0,0,0,0,0,0,0,0))
arrbins_d0 = array('d',(range(-100,100,10)))
eff_barcode_d0 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
eff_norm_d0 = array('d',(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))


barcodeeff_pT = ROOT.TH1F("eff_pT","Barcode Efficiency vs. pT;pT;Frac. of Events with >0.5 Roads",len(arrbins)-1,arrbins)

barcodeeff_d0 = ROOT.TH1F("eff_d0","Barcode Efficiency vs. d0;d0;Frac. of Events with >0.5 Roads",len(arrbins_d0)-1,arrbins_d0)

barcodefrac_type0 = ROOT.TH1F("eff_type0","Barcode Fraction 1 Hit/Layer;Frac of Roads with >0.5;Events",10,0.,1.)
barcodefrac_type1 = ROOT.TH1F("eff_type1","Barcode Fraction Missing Layers Only;Frac of Roads with >0.5;Events",10,0.,1.)
barcodefrac_type2 = ROOT.TH1F("eff_type2","Barcode Fraction Extra Hits Only;Frac of Roads with >0.5;Events",10,0.,1.)
barcodefrac_type3 = ROOT.TH1F("eff_type3","Barcode Fraction Both Extra + Missing;Frac of Roads with >0.5;Events",10,0.,1.)

roadtot_all = ROOT.TH1F("roadtot_all","Number of Roads Overall;Number of Roads;Events",100,0.,100.)

roadtot_type0 = ROOT.TH1F("roadtot_type0","Number of Roads 1 Hit/Layer;Number of Roads;Events",100,0.,100.)
roadtot_type1 = ROOT.TH1F("roadtot_type1","Number of Roads Missing Layers Only;Number of Roads;Events",100,0.,100.)
roadtot_type2 = ROOT.TH1F("roadtot_type2","Number of Roads Extra Hits Only;Number of Roads;Events",100,0.,100.)
roadtot_type3 = ROOT.TH1F("roadtot_type3","Number of Roads Both Extra+Missing;Number of Roads;Events",100,0.,100.)

roadfrac_all = ROOT.TH1F("roadfrac_all","Number of >0.5 Roads Overall;Number of Roads;Events",100,0.,100.)

roadfrac_type0 = ROOT.TH1F("roadfrac_type0","Number of >0.5 Roads 1 Hit/Layer;Number of Roads;Events",100,0.,100.)
roadfrac_type1 = ROOT.TH1F("roadfrac_type1","Number of >0.5 Roads Missing Layers Only;Number of Roads;Events",100,0.,100.)
roadfrac_type2 = ROOT.TH1F("roadfrac_type2","Number of >0.5 Roads Extra Hits Only;Number of Roads;Events",100,0.,100.)
roadfrac_type3 = ROOT.TH1F("roadfrac_type3","Number of >0.5 Roads Both Extra+Missing;Number of Roads;Events",100,0.,100.)

ngoodbarcode = 0
efficiency = 0

for i in range(len(event)):
    if(roadbfrac[i] > 0):
        ngoodbarcode += 1

print(len(event))
print(ngoodbarcode)

efficiency = float(ngoodbarcode)/float(len(event))
print(efficiency)

for i in range(len(event)):
    p = pT[i]/1000.
    for j in range(len(eff_barcode_pT)):
        if p >= arrbins[j] and p < arrbins[j+1]:
            eff_norm_pT[j] += 1
            if(roadbfrac[i] > 0):
                eff_barcode_pT[j] += 1
    for k in range(len(eff_barcode_d0)):
        if d0[i] >= arrbins_d0[k] and d0[i] < arrbins_d0[k+1]:
            eff_norm_d0[k] += 1
            if(roadbfrac[i] > 0):
                eff_barcode_d0[k] +=1

for i in range(len(eff_barcode_pT)):
    print(eff_barcode_pT[i])
    if(eff_norm_pT[i] != 0):
        eff_barcode_pT[i] = float(eff_barcode_pT[i])/float(eff_norm_pT[i])
        print(eff_barcode_pT[i])
    barcodeeff_pT.SetBinContent(i+1,eff_barcode_pT[i])

for i in range(len(eff_barcode_d0)):
    if(eff_norm_d0[i] != 0):
        eff_barcode_d0[i] = float(eff_barcode_d0[i])/float(eff_norm_d0[i])
    barcodeeff_d0.SetBinContent(i+1,eff_barcode_d0[i])

c_pT= ROOT.TCanvas("Canvas_pT","Canvas_pT",800,800)
barcodeeff_pT.Draw()

c_d0 = ROOT.TCanvas("Canvas_d0","Canvas_d0",800,800)
barcodeeff_d0.Draw()

out_file = ROOT.TFile("road_barcodeeff_nopileup_2around_avg.root","RECREATE")
barcodeeff_pT.Write()
barcodeeff_d0.Write()

