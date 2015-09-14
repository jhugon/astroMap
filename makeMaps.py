#!/usr/bin/env python

import sys
import re
import urllib2
import gzip
import types
from mpl_toolkits.basemap import Basemap
import numpy as numpy
import matplotlib.pyplot as mpl
from matplotlib import rcParams
import catalogCrossRef
from constNames import ConstNames

def polarAxisWrapper(axis,projection,zeroAt=0.):
  """
  args:
    axis: The axis or Basemap instance
    projection: a string describing the map projection
    zeroAt: At what position you want zero to be at in degrees
  """
  def project(self,xList,yList):
    if isinstance(self,Basemap):
      return self(xList, yList)
    mapRho = None
    mapTheta = None
    if self.projection == "npaeqd":
      mapTheta = (-xList+self.zeroAt)*numpy.pi/180.
      mapRho = 90 - yList
    elif self.projection == "spaeqd":
      mapTheta = (xList+self.zeroAt)*numpy.pi/180.
      mapRho = 90 + yList 
    return mapTheta, mapRho

  axis.project = types.MethodType(project,axis)
  axis.rho_lim = None
  if not isinstance(axis,Basemap):
    axis.projection = projection
    axis.zeroAt = zeroAt
    # y-axis
    ylim = axis.get_ylim()
    tickPos = []
    tickLabel = []
    for pos in numpy.arange(ylim[0],ylim[1],10):
      tickPos.append(pos)
      if axis.projection[:2] == "np":
        tickLabel.append(u"{0:.0f}\xb0".format(90.-pos))
      elif axis.projection[:2] == "sp":
        tickLabel.append(u"{0:.0f}\xb0".format(-90.+pos))
      else:
        raise Exception("Not np or sp")
    tickPos.pop(0)
    tickLabel.pop(0)
    axis.set_rgrids(tickPos,labels=tickLabel,angle=0.)
    if axis.projection[:2] == "np":
      axis.rho_lim = [90.-ylim[1],90.-ylim[0]]
    elif axis.projection[:2] == "sp":
      axis.rho_lim = [-90.+ylim[0],-90.+ylim[1]]
    # x-axis
    tickPos = []
    tickLabel = []
    nTicks = 24
    for iTick in numpy.arange(nTicks):
      pos = iTick*360./nTicks
      label = iTick*360./nTicks-axis.zeroAt
      if label < 0:
        label += 360.
      elif label > 360.:
        label -= 360.
      if axis.projection[:2] == "np":
        label = 360.-label
      label /= 15.
      tickPos.append(pos)
      tickLabel.append(u"{0:.0f}h".format(label))
    axis.set_thetagrids(tickPos,labels=tickLabel,frac=1.05)
  return axis

def drawLinesAroundBounderies(ax,xs,ys,color='k',linestyle="-",marker=None,alpha=1.0,tooFar = 180.):
  def getSlopeIntercept(x1,y1,x2,y2):
    slope = (y2-y1)/(x2-x1)
    intercept = y1 - slope*x1
    return slope, intercept

  assert(len(xs)==len(ys))

  badList = [
    "nplaea",
    "splaea",
    "npstere",
    "spstere",
    "npaeqd",
    "spaeqd",
  ]

  if not isinstance(ax,Basemap):
    ax.plot(xs,ys,alpha=alpha,color=color,linestyle=linestyle,marker=marker)
    return 

  for badProj in badList:
    if ax.projection == badProj:
      ax.plot(xs,ys,alpha=alpha,color=color,linestyle=linestyle,marker=marker)
      return 

  nPoints = len(xs)
  if nPoints > 1:
    #print min(xs),min(ys),max(xs),max(ys)
    #print m(min(xs),min(ys),inverse=True),m(max(xs),max(ys),inverse=True)
    #print m.llcrnrlon,m.llcrnrlat,m.urcrnrlon,m.urcrnrlat
    #print m(m.llcrnrlon,m.llcrnrlat),m(m.urcrnrlon,m.urcrnrlat)
    #print "###############################################################"
    xsLists = [[]]
    ysLists = [[]]
    for i in range(nPoints):
      x = xs[i]
      y = ys[i]
      if i >= nPoints-1:
        xsLists[-1].append(x)
        ysLists[-1].append(y)
      else:
        xNext = xs[i+1]
        yNext = ys[i+1]
        xDeg, yDeg = ax(x,y,inverse=True)
        xNextDeg, yNextDeg = ax(xNext,yNext,inverse=True)
        xdiff = xNextDeg - xDeg
        #xLeftDeg = None
        #yLeftDeg = None
        #xRightDeg = None
        #yRightDeg = None
        if abs(xdiff) < 180:
          pass
          xsLists[-1].append(x)
          ysLists[-1].append(y)
        else:
          if xDeg > xNextDeg:
            #print xDeg,yDeg,xNextDeg,yNextDeg
            #####Left
            xDegM360 = xDeg - 360.
            slope, yintercept = getSlopeIntercept(xDegM360,yDeg,xNextDeg,yNextDeg)
            xLeftDeg = ax.llcrnrlon+1e-6
            yLeftDeg = slope*xLeftDeg+yintercept
            #print "leftDeg: ",xLeftDeg,yLeftDeg
            ####Right
            xNextDegP360 = xNextDeg + 360.
            #print xNextDegP360
            slope, yintercept = getSlopeIntercept(xDeg,yDeg,xNextDegP360,yNextDeg)
            xRightDeg = ax.urcrnrlon-1e-6
            yRightDeg = slope*xRightDeg+yintercept
            #print "rightDeg: ",xRightDeg,yRightDeg
            xLeft,yLeft = ax(xLeftDeg,yLeftDeg)
            xRight,yRight = ax(xRightDeg,yRightDeg)
            #### Append Points
            xsLists[-1].append(x)
            ysLists[-1].append(y)
            xsLists[-1].append(xRight)
            ysLists[-1].append(yRight)
            xsLists.append([xLeft])
            ysLists.append([yLeft])
          else:
            #print xDeg,yDeg,xNextDeg,yNextDeg
            #####Left
            xNextDegM360 = xNextDeg - 360.
            #print xNextDegM360
            slope, yintercept = getSlopeIntercept(xDeg,yDeg,xNextDegM360,yNextDeg)
            xLeftDeg = ax.llcrnrlon+1e-6
            yLeftDeg = slope*xLeftDeg+yintercept
            #print "leftDeg: ",xLeftDeg,yLeftDeg
            ####Right
            xDegP360 = xDeg + 360.
            #print xDegP360
            slope, yintercept = getSlopeIntercept(xDegP360,yDeg,xNext,yNextDeg)
            xRightDeg = ax.urcrnrlon-1e-6
            yRightDeg = slope*xRightDeg+yintercept
            #print "rightDeg: ",xRightDeg,yRightDeg
            xLeft,yLeft = ax(xLeftDeg,yLeftDeg)
            xRight,yRight = ax(xRightDeg,yRightDeg)
            #### Append Points
            xsLists[-1].append(x)
            ysLists[-1].append(y)
            xsLists[-1].append(xLeft)
            ysLists[-1].append(yLeft)
            xsLists.append([xRight])
            ysLists.append([yRight])
    for xsToPlot, ysToPlot in zip(xsLists,ysLists):
      ax.plot(xsToPlot,ysToPlot,alpha=alpha,color=color,linestyle=linestyle,marker=marker)

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

def drawConstLines(baseMap,color="k"):
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
      mapx,mapy = baseMap.project(ra,de)
      starXs.append(mapx)
      starYs.append(mapy)
    drawLinesAroundBounderies(baseMap,starXs,starYs,color=color,linestyle='-',marker=None)

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
        DE = float(row[11:22])
        cst = row[23:27].strip()
        #pointType = float(row[28])
        RA = RAh * 15.
        if RA > 180.:
            RA -= 360.
        if constBoundRaw.has_key(cst):
            constBoundRaw[cst].append([RA,DE])
        else:
            constBoundRaw[cst] = [[RA,DE]]
    infile.close()
    self.constBoundRaw = constBoundRaw
  def draw(self,baseMap,color="k"):
    constBoundRaw = self.constBoundRaw
    for cst in constBoundRaw:
      mapXs = []
      mapYs = []
      for point in constBoundRaw[cst]:
        mapx,mapy = baseMap.project(point[0],point[1])
        mapXs.append(mapx)
        mapYs.append(mapy)
      #m.plot(mapXs,mapYs,"--c",alpha=0.7)
      drawLinesAroundBounderies(baseMap,mapXs,mapYs,color=color,linestyle='-',marker=None)

class StarMapper(object):

  def __init__(self):
    self.setupData()

  def setupData(self):
    dataObjs = readHip()
    dataArray = numpy.zeros((len(dataObjs),3))
    for i,hd in enumerate(dataObjs):
      if hd.getVmag() < 7.:
        dataArray[i][0] = hd.getRA()
        dataArray[i][1] = hd.getDE()
        dataArray[i][2] = hd.getVmag()
    self.dataArray = dataArray
    
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
    self.ngcGx = numpy.array(ngcGx)
    self.ngcOC = numpy.array(ngcOC)
    self.ngcGb = numpy.array(ngcGb)
    self.ngcNb = numpy.array(ngcNb)

    #print len(self.ngcOC)
    #self.ngcOC = self.ngcOC[100:150]

  def createMap(self,basemapArgs):
    if type(basemapArgs) != dict:
      raise Exception("StarMapper.createMap requires dict argument of Basemap args")
    basemapArgs['celestial'] = True
    m = Basemap(**basemapArgs)
    return m

  def drawStars(self,basemap):
    mapx,mapy = basemap.project(self.dataArray[:,0],self.dataArray[:,1])
    basemap.scatter(mapx,mapy,s=10./(numpy.sqrt(self.dataArray[:,2])),marker=".",c='k',linewidths=0)

  def drawConsts(self,basemap):
    ConstBoundaries().draw(basemap,color="0.7")
    drawConstLines(basemap,color="0.3")

  def drawGx(self,basemap,c='r'):
    mapx,mapy = basemap.project(self.ngcGx[:,0],self.ngcGx[:,1])
    basemap.scatter(mapx,mapy,marker=".",c=c,linewidths=0)
  def drawOC(self,basemap,c='g'):
    mapx,mapy = basemap.project(self.ngcOC[:,0],self.ngcOC[:,1])
    basemap.scatter(mapx,mapy,marker=".",c=c,linewidths=0)
  def drawGb(self,basemap,c='b'):
    mapx,mapy = basemap.project(self.ngcGb[:,0],self.ngcGb[:,1])
    basemap.scatter(mapx,mapy,marker=".",c=c,linewidths=0)
  def drawNb(self,basemap,c='m'):
    mapx,mapy = basemap.project(self.ngcNb[:,0],self.ngcNb[:,1])
    basemap.scatter(mapx,mapy,marker=".",c=c,linewidths=0)

  def drawGrid(self,basemap,label=True):
    if isinstance(basemap,Basemap):
      labels=[False,False,False,False]
      if label:
        labels=[True,True,False,False]
      basemap.drawparallels(numpy.arange(-90.,91.,10.),
                            labels=labels,
                            dashes=[],
                            color = '0.8'#,
                            #linewidth = 1.
        )
      labels=[False,False,False,False]
      if label:
        labels=[False,False,True,True],
      basemap.drawmeridians(numpy.arange(-180.,181.,360/24.),
                            labels=labels,
                            fmt=lambda x: "{0:.0f}h".format(x/15.),
                            dashes=[],
                            color = '0.8'#,
                            #linewidth = 1.
        )

  def drawEcliptic(self,basemap,color="0.8",linestyle="-",marker=None):
    xData = numpy.linspace(-180,180,1000)
    yData = 23.44*numpy.sin(xData*numpy.pi/180.)
    xsToPlot, ysToPlot = basemap.project(xData,yData)
    basemap.plot(xsToPlot,ysToPlot,color=color,linestyle=linestyle,marker=marker)

if __name__ == "__main__":

  rcParams["font.size"] = 20.0
  rcParams["font.family"] = "Optima Nova LT Pro"
  rcParams["grid.linestyle"] = "-"
  rcParams["grid.linewidth"] = 1
  rcParams["grid.color"] = "0.8"

  sm = StarMapper()
  cns = ConstNames()

  basemapArgs = {
  
    'projection':'robin',
    #'projection':'kav7',
    'lat_0':0,
    'lon_0':0,
  
#    'projection':'nplaea',
#    'lon_0':0,
#    'boundinglat':30,
#  
#    'projection':'splaea',
#    'lon_0':0,
#    'boundinglat':-30,
#  
#    #'projection':'cyl',
#    #'projection':'merc',
#    'projection':'mill',
#    'llcrnrlat':-70,
#    'urcrnrlat':70,
#    'llcrnrlon':-180,
#    'urcrnrlon':180,
#    'lat_ts':20,
  
  }


  #fig, ax = mpl.subplots()
  #basemapArgs['ax']=ax

  #m = sm.createMap(basemapArgs)
  ##sm.drawStars(m)
  ##sm.drawGb(m)
  ##sm.drawGx(m)
  ##sm.drawNb(m)
  ##sm.drawOC(m)
  #sm.drawConsts(m)
  #
  #fig.savefig("map.pdf")
  #fig.savefig("map.png")

  ######################################################
  ######################################################
  ######################################################

  fig = mpl.figure(figsize=(36,48),dpi=300)
  axMain = fig.add_axes([0.07,0.23,0.86,0.54]) # left, bottom, width, height in fraction of fig
  axNP = fig.add_axes([0.07,0.69,0.86,0.28],projection="polar") # left, bottom, width, height in fraction of fig
  axSP = fig.add_axes([0.07,0.03,0.86,0.28],projection="polar") # left, bottom, width, height in fraction of fig
  axSP.set_ylim(0,50)
  axSP.projection = "spaeqd"
  axNP.set_ylim(0,50)
  axNP.projection = "npaeqd"

  mMain = sm.createMap({
    'projection':'cyl',
    #'projection':'merc',
    #'projection':'mill',
    'llcrnrlat':-80,
    'urcrnrlat':80,
    'llcrnrlon':-180,
    'urcrnrlon':180,
    'lat_ts':20,
    'ax':axMain,
  })

  polarAxisWrapper(mMain,"")
  polarAxisWrapper(axNP,"npaeqd",zeroAt=-90.)
  polarAxisWrapper(axSP,"spaeqd",zeroAt=90.)

  maps = [mMain,axNP,axSP]
  for m in maps:
    sm.drawGrid(m)
    sm.drawEcliptic(m)
    sm.drawConsts(m)
    cns.drawConstNames(m,fontsize="small",color="0.3")
    #sm.drawGb(m)
    #sm.drawGx(m)
    #sm.drawNb(m)
    #sm.drawOC(m)
    #sm.drawStars(m)

  fig.text(0.93,0.335,"Midnight\nZenith Month",size="small",ha="center",va="center")
  dataToDisplay = mMain.ax.transData
  displayToFigure = fig.transFigure.inverted()
  for iMonth in range(0,12):
    iDegrees = iMonth*30
    if iDegrees > 180:
      iDegrees -= 360
    xyAxis = mMain(iDegrees,0)
    xyDisplay = dataToDisplay.transform(xyAxis)
    xyFig = displayToFigure.transform(xyDisplay)
    iMonth += 10
    iMonth %= 12
    months = ["Dec","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov"]
    fig.text(xyFig[0],0.335,"{0}".format(months[iMonth]),ha="center",va="center")

  for location,y in [("NM Skies",34.5),("GNV",29.5),("Spain",38.),("Australia",-31.25)]:
    xyAxis = mMain(0,y)
    xyDisplay = dataToDisplay.transform(xyAxis)
    xyFig = displayToFigure.transform(xyDisplay)
    fig.text(0.04,xyFig[1],"{0}".format(location),ha="right",va="center",size="small")

  #################################################

  axGb = fig.add_axes([0.75,0.03,0.2,0.08]) # left, bottom, width, height in fraction of fig
  axGx = fig.add_axes([0.75,0.13,0.2,0.08]) # left, bottom, width, height in fraction of fig
  axOC = fig.add_axes([0.05,0.03,0.2,0.08]) # left, bottom, width, height in fraction of fig
  axNb = fig.add_axes([0.05,0.13,0.2,0.08]) # left, bottom, width, height in fraction of fig

  fig.text(0.85,0.11,"Globular Clusters",ha="center",va="bottom",size="large")
  fig.text(0.85,0.21,"Galaxies",ha="center",va="bottom",size="large")
  fig.text(0.15,0.11,"Open Clusters",ha="center",va="bottom",size="large")
  fig.text(0.15,0.21,"Nebulae",ha="center",va="bottom",size="large")

  mNb = sm.createMap({
    'projection':'robin',
    'lat_0':0,
    'lon_0':0,
    'ax':axNb,
  })
  mOC = sm.createMap({
    'projection':'robin',
    'lat_0':0,
    'lon_0':0,
    'ax':axOC,
  })
  mGb = sm.createMap({
    'projection':'robin',
    'lat_0':0,
    'lon_0':0,
    'ax':axGb,
  })
  mGx = sm.createMap({
    'projection':'robin',
    'lat_0':0,
    'lon_0':0,
    'ax':axGx,
  })

  polarAxisWrapper(mNb,"")
  polarAxisWrapper(mOC,"")
  polarAxisWrapper(mGb,"")
  polarAxisWrapper(mGx,"")

  sm.drawGrid(mNb,False)
  sm.drawGrid(mOC,False)
  sm.drawGrid(mGb,False)
  sm.drawGrid(mGx,False)

  sm.drawNb(mNb,c='k')
  sm.drawOC(mOC,c='k')
  sm.drawGb(mGb,c='k')
  sm.drawGx(mGx,c='k')

  #################################################

  fig.savefig('map.png')
  fig.savefig('map.pdf')
