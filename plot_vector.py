from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH2I
#plt.rcParams["figure.figsize"] = [7.00, 3.50]
#plt.rcParams["figure.autolayout"] = True
#columns = ["Name", "Marks"]

data = pd.read_csv("1_outputvec_nopileup_withOuterLayers.txt",sep=',',header=None)
data.apply(pd.to_numeric)

print(data[2])
print(type(data[2]))

m_d0_range = 120
m_qOverPt_range = 0.002 

image_size_x = 432
image_size_y = 432

z = pd.Series(data[2]).values
x = pd.Series(data[0]).values
y = pd.Series(data[1]).values
hit1 = pd.Series(data[3]).values
hit2 = pd.Series(data[4]).values
print(z)
print(len(z))
h1 = ROOT.TH2I("h1","Hough Transform;d0;qOverPt",image_size_x,-m_d0_range,m_d0_range,image_size_y,-m_qOverPt_range,m_qOverPt_range)

h_road1 = ROOT.TH2I("h_road1","Hough Transform;d0;qOverPt",image_size_x,-m_d0_range,m_d0_range,image_size_y,-m_qOverPt_range,m_qOverPt_range)

h_road2 = ROOT.TH2I("h_road2","Hough Transform;d0;qOverPt",image_size_x,-m_d0_range,m_d0_range,image_size_y,-m_qOverPt_range,m_qOverPt_range)

h_road3 = ROOT.TH2I("h_road3","Hough Transform;d0;qOverPt",image_size_x,-m_d0_range,m_d0_range,image_size_y,-m_qOverPt_range,m_qOverPt_range)

for i in range(len(z)):
    h1.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 1):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 2):
        h_road1.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 3):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 4):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 5):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 0 and hit2[i] == 6):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 1 and hit2[i] == 3):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 1 and hit2[i] == 4):
        h_road1.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 1 and hit2[i] == 5):
        h_road2.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))
    if(hit1[i] == 3 and hit2[i] == 5):
        h_road1.SetBinContent(int(x[i])+1,int(y[i])+1,float(z[i]))

c= ROOT.TCanvas("Canvas","Canvas",800,800)
h1.Draw("COLZ")

c_road1 = ROOT.TCanvas("c_road1","c_road1",800,800)

c_road2 = ROOT.TCanvas("c_road2","c_road2",800,800)

#c_road3 = ROOT.TCanvas("c_road3","c_road3",800,800)
#h_road3.Draw("COLZ")

out_file = ROOT.TFile("1_image_nopileup_thresh15_hits_432.root","RECREATE")
h1.Write()
h_road1.Write()
h_road2.Write()
#h_road3.Write()
#with open('outputvec.txt','r') as datafile:
#     plotting = csv.reader(datafile,delimiter=',')
#     for ROWS in plotting:
#         print(ROWS)
#         x_0 = ROWS[0]
#         y_1 = ROWS[1]
#         z = ROWS[2]
         #x.append(float(ROWS[0]))
         #y.append(float(ROWS[1]))
         #np.append(z,ROWS[])
#print(y_0)
fig = plt.figure
ax = plt.axes(projection='3d')

#ax.plot_surface(data[0],data[1],data[2],edgecolor='green')
#plt.show
#df = pd.read_csv("outputvec.txt", usecols=columns)
#plt.plot(df.Name, df.Marks)
#plt.show()
