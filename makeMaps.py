#!/usr/bin/env python

import sys
import re
import urllib2
import gzip
from mpl_toolkits.basemap import Basemap
import numpy as numpy
import matplotlib.pyplot as plt
import catalogCrossRef

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
    result = self.RA
    if result > 180.:
      result -= 360.
    return result

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

def readHip(localFn="data/hip_main.dat",url="ftp://cdsarc.u-strasbg.fr/pub/cats/I/239/hip_main.dat.gz"):
  result = []
  infile = None
  try:
    infile = gzip.GzipFile(localFn)
  except:
    print "readHip: Couldn't find local data file, attempting to download..."
    infile = None
    try:
      u = urllib2.urlopen(url)
      infile = open(localFn, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (url, file_size)
      
      file_size_dl = 0
      block_sz = 8192
      while True:
          buffer = u.read(block_sz)
          if not buffer:
              break
          file_size_dl += len(buffer)
          infile.write(buffer)
          status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
          status = status + chr(8)*(len(status)+1)
          print status,
      infile.close()
      print "Done downloading %s to %s" % (url, localFn)
      infile = gzip.GzipFile(localFn)
    except:
      print "Error in readHip: Couldn't find local data file and couldn't download"
      sys.exit(1)
  for row in infile:
    result.append(HipEntry(row))
  infile.close()
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
    result = self.RA
    if result > 180.:
      result -= 360.
    return result

  def getDE(self):
    return self.DE

  def getType(self):
    return self.Type

  def getNGC(self):
    return self.NGC

def readNGC(localFn="data/ngc2000.dat",url="ftp://cdsarc.u-strasbg.fr/pub/cats/VII/118/ngc2000.dat"):
  result = []
  infile = None
  try:
    infile = open(localFn)
  except:
    print "readNGC: Couldn't find local data file, attempting to download..."
    infile = None
    try:
      u = urllib2.urlopen(url)
      infile = open(localFn, 'wb')
      meta = u.info()
      file_size = int(meta.getheaders("Content-Length")[0])
      print "Downloading: %s Bytes: %s" % (url, file_size)
      
      file_size_dl = 0
      block_sz = 8192
      while True:
          buffer = u.read(block_sz)
          if not buffer:
              break
          file_size_dl += len(buffer)
          infile.write(buffer)
          status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
          status = status + chr(8)*(len(status)+1)
          print status,
      infile.close()
      print "Done downloading %s to %s" % (url, localFn)
      infile = open(localFn)
    except:
      print "Error in readNGC: Couldn't find local data file and couldn't download"
      sys.exit(1)
  for row in infile:
    result.append(NGCEntry(row))
  infile.close()
  return result

def drawConstLines(baseMap):
  ccr = catalogCrossRef.CatalogCrossRef()
  lineFile = open("data/pp3Lines.dat")
  goodLines = []
  for l in lineFile:
    if l[0]=="#" or re.match(r"^\s*$",l):
        continue
    l = re.sub(r"#.*","",l)
    l = l.replace('\n','')
    goodLines.append(l)
  lines = " ".join(goodLines).split(";")
  for line in lines:
    stars = re.findall(r"([a-zA-Z]+)\s+([0-9]+)",line)
    #print "line: "
    starXs = []
    starYs = []
    for star in stars:
      const= star[0]
      num = int(star[1])
      if const == "HD":
        star = ccr.findByHD(num)
      else:
        star = ccr.findByFl(const,num)
      #print " ", const,num,star.getRAd(),star.getDEd()
      ra = star.getRAd()
      if ra > 180.:
        ra -= 360.
      de = star.getDEd()
      mapx,mapy = m(ra,de)
      starXs.append(mapx)
      starYs.append(mapy)
    m.plot(starXs,starYs,"-m",alpha=0.7)

class ConstBoundaries(object):
  def __init__(self,localFn="data/bound_20.dat",url="ftp://cdsarc.u-strasbg.fr/cats/VI/49/bound_20.dat"):
    infile = None
    try:
      infile = open(localFn)
    except:
      print "ConstBoundaries: Couldn't find local data file, attempting to download..."
      infile = None
      try:
        u = urllib2.urlopen(url)
        infile = open(localFn, 'wb')
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        print "Downloading: %s Bytes: %s" % (url, file_size)
        
        file_size_dl = 0
        block_sz = 8192
        while True:
            buffer = u.read(block_sz)
            if not buffer:
                break
            file_size_dl += len(buffer)
            infile.write(buffer)
            status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
            status = status + chr(8)*(len(status)+1)
            print status,
        infile.close()
        print "Done downloading %s to %s" % (url, localFn)
        infile = open(localFn)
      except:
        print "Error in ConstBoundaries: Couldn't find local data file and couldn't download"
        sys.exit(1)
    constBoundRaw = {}
    for row in infile:
        RAh = float(row[0:10])
        if RAh > 180.:
            RAh -= 360.
        DE = float(row[11:22])
        cst = row[23:27].strip()
        #pointType = float(row[28])
        RA = RAh * 15.
        if constBoundRaw.has_key(cst):
            constBoundRaw[cst].append([RA,DE])
        else:
            constBoundRaw[cst] = [[RA,DE]]
    infile.close()
    self.constBoundRaw = constBoundRaw
  def draw(self,baseMap):
    constBoundRaw = self.constBoundRaw
    for cst in constBoundRaw:
      mapXs = []
      mapYs = []
      for point in constBoundRaw[cst]:
        mapx,mapy = m(point[0],point[1])
        mapXs.append(mapx)
        mapYs.append(mapy)
      m.plot(mapXs,mapYs,"--c",alpha=0.7)

dataObjs = readHip()
dataArray = numpy.zeros((len(dataObjs),3))
for i,hd in enumerate(dataObjs):
  if hd.getVmag() < 7.:
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
#m = Basemap(projection='kav7',lat_0=0,lon_0=0,celestial=True)
m = Basemap(projection='robin',lat_0=0,lon_0=0,celestial=True)
#m = Basemap(projection='nplaea',boundinglat=30,lon_0=0,celestial=True)
#m = Basemap(projection='splaea',boundinglat=-30,lon_0=0,celestial=True)
## Cylindrical
projection = 'merc'
#m = Basemap(projection=projection,llcrnrlat=-70,urcrnrlat=70,
#            llcrnrlon=-180,urcrnrlon=180,lat_ts=20,celestial=True)
#m = Basemap(projection='cyl',lat_0=0,lon_0=0,celestial=True)
#m = Basemap(projection='merc',lat_0=0,lon_0=0,celestial=True)
#m = Basemap(projection='mill',lat_0=0,lon_0=0,celestial=True)
##
#m = Basemap(celestial=True,
#        #projection="laea",
#        #projection="lcc",
#        projection="aeqd",
#        #projection='stere',
#        rsphere=1./(2.*numpy.pi),
#        ## LMC-ish
#        #width=0.10,
#        #height=0.1,
#        #lat_0=-70,
#        #lon_0=-90,
#        ## Orion-ish
#        #width=0.1,
#        #height=0.1,
#        #lat_0=0,
#        #lon_0=-85,
#        width=0.3,
#        height=0.3,
#        lat_0=21.6,
#        lon_0=114,
#    )
# draw parallels and meridians.
m.drawparallels(numpy.arange(-90.,91.,30.),labels=[True,True,True])
m.drawmeridians(numpy.arange(-180.,181.,60.),labels=[True,True,True])
mapx,mapy = m(dataArray[:,0],dataArray[:,1])
m.scatter(mapx,mapy,s=10./(numpy.sqrt(dataArray[:,2])),marker=".",c='k',linewidths=0)
#mapx,mapy = m(ngcGx[:,0],ngcGx[:,1])
#m.scatter(mapx,mapy,marker=".",c='r',linewidths=0)
#mapx,mapy = m(ngcOC[:,0],ngcOC[:,1])
#m.scatter(mapx,mapy,marker=".",c='g',linewidths=0)
#mapx,mapy = m(ngcGb[:,0],ngcGb[:,1])
#m.scatter(mapx,mapy,marker=".",c='b',linewidths=0)
#mapx,mapy = m(ngcNb[:,0],ngcNb[:,1])
#m.scatter(mapx,mapy,marker=".",c='c',linewidths=0)
drawConstLines(m)
cbs = ConstBoundaries()
cbs.draw(m)
plt.show()


