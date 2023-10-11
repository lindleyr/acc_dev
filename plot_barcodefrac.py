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
barcode = pd.Series(data[7]).values
charge = pd.Series(data[8]).values
pT = pd.Series(data[9]).values
d0 = pd.Series(data[10]).values


barcodefrac_all = ROOT.TH1F("eff_all","Barcode Fraction Overall;Frac of Roads with >0.5;Events",10,0.,1.)

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

nroadevts = 0

for i in range(len(event)):
    #print('Event: '+str(event[i])+' | i: '+str(i)
    if(i == 0):
        roadtot_all.Fill(roadtot[i])
        roadfrac_all.Fill(roadbfrac[i])
        if(eventtype[i] == 0):
            roadtot_type0.Fill(roadtot[i])
            roadfrac_type0.Fill(roadbfrac[i])
        if(eventtype[i] == 1):
            roadtot_type1.Fill(roadtot[i])
            roadfrac_type1.Fill(roadbfrac[i])
        if(eventtype[i] == 2):
            roadtot_type2.Fill(roadtot[i])
            roadfrac_type2.Fill(roadbfrac[i])
        if(eventtype[i] == 3):
            roadtot_type3.Fill(roadtot[i])
            roadfrac_type3.Fill(roadbfrac[i])
    else:
        if(event[i-1] == event[i]):
           continue
        if(event[i-1] != event[i]):
           roadtot_all.Fill(roadtot[i])
           roadfrac_all.Fill(roadbfrac[i])
           if(eventtype[i] == 0):
               roadtot_type0.Fill(roadtot[i])
               roadfrac_type0.Fill(roadbfrac[i])
           if(eventtype[i] == 1):
               roadtot_type1.Fill(roadtot[i])
               roadfrac_type1.Fill(roadbfrac[i])
           if(eventtype[i] == 2):
               roadtot_type2.Fill(roadtot[i])
               roadfrac_type2.Fill(roadbfrac[i])
           if(eventtype[i] == 3):
               roadtot_type3.Fill(roadtot[i])
               roadfrac_type3.Fill(roadbfrac[i])

for i in range(len(layer)):
    if(roadtot[i] != 0):
#        roadtot_all.Fill(float(roadtot[i]))
        barcodefrac_all.Fill(float(roadbfrac[i])/float(roadtot[i]))
#        print(float(roadbfrac[i])/float(roadtot[i]))
        if eventtype[i] == 0:
            barcodefrac_type0.Fill(float(roadbfrac[i])/float(roadtot[i]))
        if eventtype[i] == 1:
            barcodefrac_type1.Fill(float(roadbfrac[i])/float(roadtot[i]))
        if eventtype[i] == 2:
            barcodefrac_type2.Fill(float(roadbfrac[i])/float(roadtot[i]))
        if eventtype[i] == 3:
            barcodefrac_type3.Fill(float(roadbfrac[i])/float(roadtot[i]))
    else:
        barcodefrac_all.Fill(0)
        if eventtype[i] == 0:
            barcodefrac_type0.Fill(0)
        if eventtype[i] == 1:
            barcodefrac_type1.Fill(0)
        if eventtype[i] == 2:
            barcodefrac_type2.Fill(0)
        if eventtype[i] == 3:
            barcodefrac_type3.Fill(0)

roadfrac_all.SetLineColor(2) #2 is red
roadfrac_type0.SetLineColor(2)
roadfrac_type1.SetLineColor(2)
roadfrac_type2.SetLineColor(2)
roadfrac_type3.SetLineColor(2)

c= ROOT.TCanvas("Canvas","Canvas",800,800)
barcodefrac_all.Draw()

ctot = ROOT.TCanvas("Canvastot","Canvastot",800,800)
roadtot_all.Draw()
roadfrac_all.Draw("Same")
leg_all = ROOT.TLegend(.3,.7,.5,.9,"Number of Roads")
leg_all.AddEntry(roadtot_all,"Total")
leg_all.AddEntry(roadfrac_all,">0.5 Barcode Frac")
leg_all.Draw("Same")

ctot2 = ROOT.TCanvas("ctot2","ctot2",800,800)
roadtot_all.Draw()

c_roadfrac = ROOT.TCanvas("Canvas","Canvas",800,800)
roadfrac_all.Draw()

c0 = ROOT.TCanvas("Canvas","Canvas",800,800)
barcodefrac_type0.Draw()

c0_tot = ROOT.TCanvas("Canvastype0","Canvas",800,800)
roadtot_type0.Draw()
roadfrac_type0.Draw("Same")
leg_type0 = ROOT.TLegend(.3,.7,.5,.9,"Number of Roads")
leg_type0.AddEntry(roadtot_type0,"Total")
leg_type0.AddEntry(roadfrac_type0,">0.5 Barcode Frac")
leg_type0.Draw("Same")

c0_tot2 = ROOT.TCanvas("c0_tot2","c0_tot2",800,800)
roadtot_type0.Draw()

c0_roadfrac = ROOT.TCanvas("Canvas","Canvas",800,800)
roadfrac_type0.Draw()

c1 = ROOT.TCanvas("Canvas","Canvas",800,800)
barcodefrac_type1.Draw()

c1_tot = ROOT.TCanvas("Canvastype1","Canvas",800,800)
roadtot_type1.Draw()
roadfrac_type1.Draw("Same")
leg_type1 = ROOT.TLegend(.3,.7,.5,.9,"Number of Roads")
leg_type1.AddEntry(roadtot_type1,"Total")
leg_type1.AddEntry(roadfrac_type1,">0.5 Barcode Frac")
leg_type1.Draw("Same")

c1_roadfrac = ROOT.TCanvas("Canvas","Canvas",800,800)
roadfrac_type1.Draw()

c2 = ROOT.TCanvas("Canvas","Canvas",800,800)
barcodefrac_type2.Draw()

c2_tot = ROOT.TCanvas("Canvastype2","Canvas",800,800)
roadtot_type2.Draw()
roadfrac_type2.Draw("Same")
leg_type2 = ROOT.TLegend(.3,.7,.5,.9,"Number of Roads")
leg_type2.AddEntry(roadtot_type2,"Total")
leg_type2.AddEntry(roadfrac_type2,">0.5 Barcode Frac")
leg_type2.Draw("Same")

c2_roadfrac = ROOT.TCanvas("Canvas","Canvas",800,800)
roadfrac_type2.Draw()

c3 = ROOT.TCanvas("Canvas","Canvas",800,800)
barcodefrac_type3.Draw()

c3_tot = ROOT.TCanvas("Canvastype3","Canvas",800,800)
roadtot_type3.Draw()
roadfrac_type3.Draw("Same")
leg_type3 = ROOT.TLegend(.3,.7,.5,.9,"Number of Roads")
leg_type3.AddEntry(roadtot_type3,"Total")
leg_type3.AddEntry(roadfrac_type3,">0.5 Barcode Frac")
leg_type3.Draw("Same")

c3_roadfrac = ROOT.TCanvas("Canvas","Canvas",800,800)
roadfrac_type3.Draw()

c5 = ROOT.TCanvas("Canvas5","Canvas5",800,800)
#eff_d0_5.SetMinimum(0)
#barcodefrac_type0.SetMaximum(50*barcodefrac_type0.GetMaximum())
barcodefrac_type1.Draw()
barcodefrac_type0.SetLineColor(6)
barcodefrac_type2.SetLineColor(2)
barcodefrac_type3.SetLineColor(4)
barcodefrac_type0.Draw("Same")
barcodefrac_type2.Draw("Same")
barcodefrac_type3.Draw("Same")
leg = ROOT.TLegend(.3,.7,.5,.9,"Frac of 0.5 Frac Roads")
leg.AddEntry(barcodefrac_type0,"1 Hit/Layer")
leg.AddEntry(barcodefrac_type1,"Missing Layers Only")
leg.AddEntry(barcodefrac_type2,"Extra Hits Only")
leg.AddEntry(barcodefrac_type3,"Both Missing + Extra")
leg.Draw("Same")

out_file = ROOT.TFile("road_barcodefrac_nopileup_2around_avg.root","RECREATE")
ctot.Write()
c0_tot.Write()
c1_tot.Write()
c2_tot.Write()
c3_tot.Write()
roadtot_all.Write()
roadtot_type0.Write()
roadtot_type1.Write()
roadtot_type2.Write()
roadtot_type3.Write()
roadfrac_all.Write()
roadfrac_type0.Write()
roadfrac_type1.Write()
roadfrac_type2.Write()
roadfrac_type3.Write()
barcodefrac_all.Write()
barcodefrac_type0.Write()
barcodefrac_type1.Write()
barcodefrac_type2.Write()
barcodefrac_type3.Write()
c5.Write()

