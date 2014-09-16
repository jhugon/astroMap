#!/usr/bin/env python

from mpl_toolkits.basemap import Basemap
import numpy as numpy
import matplotlib.pyplot as plt

class HipEntry(object):
  def __init__(self,row):
    assert(row[0]=="H")
    self.HIP = int(row[8:14])  #hipparcos number
    self.Vmag = row[41:46]  # Johnson V Magnitude
    self.RA = row[51:63]  # RA in degrees (ICRS, J1991.25)
    self.DE = row[64:76]  # DE in degrees (ICRS, J1991.25)
    self.Plx = row[79:86]  # Trigonometric parallax in millarcseconds
    self.BmV = row[245:251]  # Johnson B-V color
    for attr in ["Vmag","RA","DE","Plx","BmV"]:
      tmpMember = getattr(self,attr)
      try:
        setattr(self,attr,float(tmpMember))
      except:
        setattr(self,attr,float("NaN"))

  def getHIP(self):
    return self.HIP

  def getRA(self):
    return self.RA

  def getDE(self):
    return self.DE

  def getVmag(self):
    return self.Vmag

  def getBminusV(self):
    return self.VmV

  def getPlx(self):
    return self.Plx

  def getDistance(self): # in parsecs
    if self.Plx == 0.:
      return float("NaN")
    return 1.0/(self.Plx*1000.)

  def __str__(self):
    result = "{0:7} {1:10.4f} {2:10.4f} {3:10.4f} {4:10.4f} {5:10.4f}".format(self.HIP,self.RA,self.DE,self.Plx,self.Vmag,self.BmV)
    return result

def readHip(infileName="hip_main.dat"):
  result = []
  infile = open(infileName)
  for row in infile:
    result.append(HipEntry(row))
  return result

class NGCEntry(object):
  def __init__(self,row):
    self.NGC = row[0:5]  # NGC or IC number
    self.Type = row[6:9]  # Type
    self.RA = float(row[10:12])+float(row[13:17])/60.  # RA in hours
    self.DEsign = row[19]  # DE in degrees
    self.DE = float(row[20:22])+float(row[23:25])/60.  # DE in degrees
    if self.DEsign == "-": # take into account sign of DE
      self.DE *= -1.
    self.RA *= 15. # convert RA from hours to degrees

  def getRA(self):
    return self.RA

  def getDE(self):
    return self.DE

  def getType(self):
    return self.Type

  def getNGC(self):
    return self.NGC

def readNGC(infileName="ngc2000.dat"):
  result = []
  infile = open(infileName)
  for row in infile:
    result.append(NGCEntry(row))
  return result

dataObjs = readHip()
dataArray = numpy.zeros((len(dataObjs),3))
for i,hd in enumerate(dataObjs):
  if hd.getVmag() < 9.:
    dataArray[i][0] = hd.getRA()
    dataArray[i][1] = hd.getDE()
    dataArray[i][2] = hd.getVmag()

ngcObjs = readNGC()
ngcGx = []
ngcOC = []
ngcGb = []
ngcNb = []
for ngc in ngcObjs:
  if "Gx" in ngc.getType():
    ngcGx.append([ngc.getRA(),ngc.getDE()])
  if "OC" in ngc.getType() or "C+N" in ngc.getType():
    ngcOC.append([ngc.getRA(),ngc.getDE()])
  if "Gb" in ngc.getType():
    ngcGb.append([ngc.getRA(),ngc.getDE()])
  if "Nb" in ngc.getType() or "Pl" in ngc.getType():
    ngcNb.append([ngc.getRA(),ngc.getDE()])
ngcGx = numpy.array(ngcGx)
ngcOC = numpy.array(ngcOC)
ngcGb = numpy.array(ngcGb)
ngcNb = numpy.array(ngcNb)

## llcrnrlat,llcrnrlon,urcrnrlat,urcrnrlon
## are the lat/lon values of the lower left and upper right corners
## of the map.
## resolution = 'c' means use crude resolution coastlines.
m = Basemap(projection='kav7',lat_0=0,lon_0=0,celestial=True)
#m = Basemap(projection='robin',lat_0=0,lon_0=90,celestial=True)
#m = Basemap(projection='nplaea',boundinglat=30,lon_0=0,celestial=True)
#m = Basemap(projection='splaea',boundinglat=-30,lon_0=0,celestial=True)
#m = Basemap(projection='cyl',celestial=True)
# draw parallels and meridians.
m.drawparallels(numpy.arange(-90.,91.,30.),labels=[True,True,True])
m.drawmeridians(numpy.arange(-180.,181.,60.),labels=[True,True,True])
mapx,mapy = m(dataArray[:,0],dataArray[:,1])
m.scatter(mapx,mapy,s=10./numpy.sqrt(dataArray[:,2]),marker=".",c='k',linewidths=0)
mapx,mapy = m(ngcGx[:,0],ngcGx[:,1])
m.scatter(mapx,mapy,marker=".",c='r',linewidths=0)
mapx,mapy = m(ngcOC[:,0],ngcOC[:,1])
m.scatter(mapx,mapy,marker=".",c='g',linewidths=0)
mapx,mapy = m(ngcGb[:,0],ngcGb[:,1])
m.scatter(mapx,mapy,marker=".",c='b',linewidths=0)
mapx,mapy = m(ngcNb[:,0],ngcNb[:,1])
m.scatter(mapx,mapy,marker=".",c='y',linewidths=0)
plt.show()


