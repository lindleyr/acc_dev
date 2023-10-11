from mpl_toolkits import mplot3d
import pandas as pd
from matplotlib import pyplot as plt
import csv
import numpy as np
import ROOT
from ROOT import TCanvas, TFile, TH2I, TGraph, TGraph2D
from array import array
#plt.rcParams["figure.figsize"] = [7.00, 3.50]
#plt.rcParams["figure.autolayout"] = True
#columns = ["Name", "Marks"]

#x1 = [0.,0.,0.,0.,0.,0.,0.,0.]
#y1 = [0.,0.,0.,0.,0.,0.,0.,0.]
#z1 = [0.,0.,0.,0.,0.,0.,0.,0.]

x1 = array('d')
y1 = array('d')
z1 = array('d')

data = pd.read_csv("txtfiles/merge_new.txt",sep=' ',header=None)
data.apply(pd.to_numeric)

print(data[3])
print(type(data[3]))

#get the correct datatype for the values
x = pd.Series(data[3]).values
y = pd.Series(data[4]).values
z = pd.Series(data[5]).values

print(x[0])
#get the first event only
for i in range(0,7,1):
  x1.append(x[i]) 
  y1.append(y[i])
  z1.append(z[i])

print(x1[0])
g1 = ROOT.TGraph(9,x1,y1)
g1.SetTitle("Hits in XY; x; y")

g1.SetLineColor( 2 )
g1.SetLineWidth( 4 )
g1.SetMarkerColor( 4 )
g1.SetMarkerStyle( 7 )

g2 = ROOT.TGraph(9,y1,z1)
g2.SetTitle("Hits in YZ; y; z")
g2.SetMarkerColor( 4 )
g2.SetMarkerStyle( 7 )


g3 = ROOT.TGraph(9,x1,z1)
g3.SetTitle("Hits in XZ; x; z")
g3.SetMarkerColor( 4 )
g3.SetMarkerStyle( 7 )

g3d = ROOT.TGraph2D(9,x1,y1,z1)
g3d.SetTitle("Hits in XYZ; x; y; z")
g3d.SetMarkerColor( 4 )
g3d.SetMarkerStyle( 7 )

c1 = TCanvas( 'c1', 'xy', 200, 10, 700, 500 )
g1.Draw( 'ACP' )
out_filexy = ROOT.TFile("event_image_1.root","RECREATE")
g1.Write("xy")

c2 = TCanvas( 'c2', 'yz', 200, 10, 700, 500 )
g2.Draw( 'AP' )
#out_fileyz = ROOT.TFile("event_image_yz.root", "RECREATE")
g2.Write("yz")

c3 = TCanvas( 'c3', 'xz', 200, 10, 700, 500 )
g3.Draw( 'AP' )
#out_filexz = ROOT.TFile("event_image_xz.root", "RECREATE")
g3.Write("xz")

c4 = TCanvas( 'c4', 'xyz', 200, 10, 700, 500 )
g3d.Draw()
g3d.Write("xyz")
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
#fig = plt.figure
#ax = plt.axes(projection='3d')

#ax.plot_surface(data[0],data[1],data[2],edgecolor='green')
#plt.show
#df = pd.read_csv("outputvec.txt", usecols=columns)
#plt.plot(df.Name, df.Marks)
#plt.show()
