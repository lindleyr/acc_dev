from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH1F, TGraph, TLegend
from array import array
import math as m
#plt.rcParams["figure.figsize"] = [7.00, 3.50]
#plt.rcParams["figure.autolayout"] = True
#columns = ["Name", "Marks"]

layer1 = 6
layer2 = 7

xlayer1 = array('d')
ylayer1 = array('d')
zlayer1 = array('d')

xlayer2 = array('d')
ylayer2 = array('d')
zlayer2 = array('d')

xlayer3 = array('d')
ylayer3 = array('d')
zlayer3 = array('d')

xlayer4 = array('d')
ylayer4 = array('d')
zlayer4 = array('d')

xlayer5 = array('d')
ylayer5 = array('d')
zlayer5 = array('d')

xlayer6 = array('d')
ylayer6 = array('d')
zlayer6 = array('d')

xlayer7 = array('d')
ylayer7 = array('d')
zlayer7 = array('d')

xlayer8 = array('d')
ylayer8 = array('d')
zlayer8 = array('d')


data = pd.read_csv("txtfiles/merge_new.txt",sep=' ',header=None)
data.apply(pd.to_numeric)

print(data[3])
print(type(data[3]))

nevents = len(data[3])
print(nevents)

#get the correct datatype for the values
x = pd.Series(data[3]).values
y = pd.Series(data[4]).values
z = pd.Series(data[5]).values
layer = pd.Series(data[1]).values

#select only events from certain layers. the 'continue' ensures we only get the first instance of a layer per event
for i in range(nevents):
  if layer[i] == 0:
    xlayer1.append(x[i])
    ylayer1.append(y[i])
    zlayer1.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 1:
    xlayer2.append(x[i])
    ylayer2.append(y[i])
    zlayer2.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 2:
    xlayer3.append(x[i])
    ylayer3.append(y[i])
    zlayer3.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 3:
    xlayer4.append(x[i])
    ylayer4.append(y[i])
    zlayer4.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 4:
    xlayer5.append(x[i])
    ylayer5.append(y[i])
    zlayer5.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 5:
    xlayer6.append(x[i])
    ylayer6.append(y[i])
    zlayer6.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 6:
    xlayer7.append(x[i])
    ylayer7.append(y[i])
    zlayer7.append(z[i])
    continue
for i in range(nevents):
  if layer[i] == 7:
    xlayer8.append(x[i])
    ylayer8.append(y[i])
    zlayer8.append(z[i])
    continue

nbins = 50
xmin = 100.
xmax = 500.

print(m.sqrt((xlayer1[0]-xlayer2[0])*(xlayer1[0]-xlayer2[0])+(ylayer1[0]-ylayer2[0])*(ylayer1[0]-ylayer2[0])+(zlayer1[0]-zlayer2[0])*(zlayer1[0]-zlayer2[0])))
h1v2 = ROOT.TH1F("h1v2","Distance;distance;events",nbins,xmin,xmax)
h2v3 = ROOT.TH1F("h2v3","Distance;distance;events",nbins,xmin,xmax)
h3v4 = ROOT.TH1F("h3v4","Distance;distance;events",nbins,xmin,xmax)
h4v5 = ROOT.TH1F("h4v5","Distance;distance;events",nbins,xmin,xmax)
h5v6 = ROOT.TH1F("h5v6","Distance;distance;events",nbins,xmin,xmax)
h6v7 = ROOT.TH1F("h6v7","Distance;distance;events",nbins,xmin,xmax)
h7v8 = ROOT.TH1F("h7v8","Distance;distance;events",nbins,xmin,xmax)


for i in range(len(min(xlayer1,xlayer2))):
  h1v2.Fill(m.sqrt((xlayer1[i]-xlayer2[i])*(xlayer1[i]-xlayer2[i])+(ylayer1[i]-ylayer2[i])*(ylayer1[i]-ylayer2[i])+(zlayer1[i]-zlayer2[i])*(zlayer1[i]-zlayer2[i])))
for i in range(len(min(xlayer2,xlayer3))):
  h2v3.Fill(m.sqrt((xlayer2[i]-xlayer3[i])*(xlayer2[i]-xlayer3[i])+(ylayer2[i]-ylayer3[i])*(ylayer2[i]-ylayer3[i])+(zlayer2[i]-zlayer3[i])*(zlayer2[i]-zlayer3[i])))
for i in range(len(min(xlayer3,xlayer4))):
  h3v4.Fill(m.sqrt((xlayer3[i]-xlayer4[i])*(xlayer3[i]-xlayer4[i])+(ylayer3[i]-ylayer4[i])*(ylayer3[i]-ylayer4[i])+(zlayer3[i]-zlayer4[i])*(zlayer3[i]-zlayer4[i])))
#for i in range(len(min(xlayer4,xlayer5))):
for i in range(len(xlayer5)):
  h4v5.Fill(m.sqrt((xlayer4[i]-xlayer5[i])*(xlayer4[i]-xlayer5[i])+(ylayer4[i]-ylayer5[i])*(ylayer4[i]-ylayer5[i])+(zlayer4[i]-zlayer5[i])*(zlayer4[i]-zlayer5[i])))
for i in range(len(min(xlayer5,xlayer6))):
  h5v6.Fill(m.sqrt((xlayer5[i]-xlayer6[i])*(xlayer5[i]-xlayer6[i])+(ylayer5[i]-ylayer6[i])*(ylayer5[i]-ylayer6[i])+(zlayer5[i]-zlayer6[i])*(zlayer5[i]-zlayer6[i])))
for i in range(len(xlayer7)):
  h6v7.Fill(m.sqrt((xlayer6[i]-xlayer7[i])*(xlayer6[i]-xlayer7[i])+(ylayer6[i]-ylayer7[i])*(ylayer6[i]-ylayer7[i])+(zlayer6[i]-zlayer7[i])*(zlayer6[i]-zlayer7[i])))
for i in range(len(min(xlayer7,xlayer8))):
  h7v8.Fill(m.sqrt((xlayer7[i]-xlayer8[i])*(xlayer7[i]-xlayer8[i])+(ylayer7[i]-ylayer8[i])*(ylayer7[i]-ylayer8[i])+(zlayer7[i]-zlayer8[i])*(zlayer7[i]-zlayer8[i])))

h1v2.SetLineColor(2)
h2v3.SetLineColor(4)
h4v5.SetLineColor(6)
h6v7.SetLineColor(8)
c1 = TCanvas( 'c1', 'xy', 200, 10, 700, 500 )
h1v2.SetMaximum(100.)
h1v2.Draw()
h2v3.Draw("Same")
h4v5.Draw("Same")
h6v7.Draw("Same")
leg = ROOT.TLegend(.1,.7,.3,.9,"DeltaR")
leg.AddEntry(h1v2,"Layer 1 vs 2")
leg.AddEntry(h2v3,"Layer 2 vs 3")
leg.AddEntry(h4v5,"Layer 4 vs 5")
leg.AddEntry(h6v7,"Layer 6 vs 7")
leg.Draw("Same")
out_file = ROOT.TFile("deltaij.root","RECREATE")
c1.Write("deltaij")

