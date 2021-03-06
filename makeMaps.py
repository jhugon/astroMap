#!/usr/bin/env python3

import sys
import re
import requests
import gzip
import types
from mpl_toolkits.basemap import Basemap
import numpy as numpy
import matplotlib.pyplot as mpl
from matplotlib import rcParams
import matplotlib.lines as mlines
import catalogCrossRef
from constNames import ConstNames

ngcTypeSymbols = [('Gx',"$⬬$","Galaxy"), ('Pl',"$+$","Planetary Nebula"), ('Nb',"$□$","Nebula"), ('C+N',"$⊡$","Cluster+Nebulosity"), ('OC',"$○$","Open Cluster"), ('Gb',"$⊕$","Globular Cluster"), ("","$∅$","Other")]

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
        tickLabel.append("{0:.0f}\xb0".format(90.-pos))
      elif axis.projection[:2] == "sp":
        tickLabel.append("{0:.0f}\xb0".format(-90.+pos))
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
      tickLabel.append("{0:.0f}h".format(label))
    axis.set_thetagrids(tickPos,labels=tickLabel)#,frac=1.1)
  return axis

def drawLinesAroundBounderies(ax,xs,ys,color='k',linestyle="-",marker=None,alpha=1.0,tooFar = 180.,linewidth=1):
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
    ax.plot(xs,ys,alpha=alpha,color=color,linestyle=linestyle,marker=marker,linewidth=linewidth)
    return 

  for badProj in badList:
    if ax.projection == badProj:
      ax.plot(xs,ys,alpha=alpha,color=color,linestyle=linestyle,marker=marker,linewidth=linewidth)
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
      ax.plot(xsToPlot,ysToPlot,alpha=alpha,color=color,linestyle=linestyle,marker=marker,linewidth=linewidth)

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

class NGCEntry(object):
  def __init__(self,row):
    self.NGC = row[0:5].replace(" ","")  # NGC or IC number
    self.Type = row[6:9].strip(" ")  # Type
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

class NGCNameEntry(object):
  def __init__(self,row):
    self.name = row[0:35].strip(" ")  # name mapping to NGC or IC number
    self.NGC = row[36:41].replace(" ","")  # NGC or IC number

  def getName(self):
    return self.name

  def getNGC(self):
    return self.NGC

class HCGEntry(object):
  """
  Hickson's Compact groups of Galaxies
  N = 100
  """
  def __init__(self,row):
    self.HCG = int(row[0:3])  # HCG number
    self.RA = float(row[4:6])+float(row[6:8])/60.+float(row[8:10])/60./60.  # RA in hours (1950)
    self.DEsign = row[10]  # DE in degrees
    self.DE = float(row[11:13])+float(row[13:15])/60.+float(row[15:17])/60./60.  # DE in degrees (1950)
    if self.DEsign == "-": # take into account sign of DE
      self.DE *= -1.
    self.RA *= 15. # convert RA from hours to degrees
    self.Type = row[19:21]  # Type
    self.MCount = int(row[21:23])  # Number of member galaxies
    self.AngSize = float(row[23:28])  # AngSize in arcmin

  def getRA(self):
    result = self.RA
    if result > 180.:
      result -= 360.
    return result

  def getDE(self):
    return self.DE

  def getType(self):
    return self.Type

  def getHCG(self):
    return self.HCG

  def getAngSize(self):
    return self.AngSize

class RADEObj:
  def __init__(self,ra,de):
    self.RA = float(ra)
    self.DE = float(de)
  def getRA(self):
    result = self.RA
    if result > 180.:
      result -= 360.
    return result
  def getDE(self):
    return self.DE

def makeMessierDict(ngcList,ngcNameList):
  """
  Dict of Messier numbers mapped to NGC entries

  M40 is missing due to it being a double star
  """
  messiers  = {}
  caldwells = {}
  ngcMap = {}
  for ngcObj in ngcList:
    ngcMap[ngcObj.getNGC()] = ngcObj
  for entry in ngcNameList:
    if entry.getName()[:2] == "M " and len(entry.getNGC()) != 0:
        messiers[int(entry.getName()[2:].strip(" "))] = ngcMap[entry.getNGC()]
  with open("data/caldwell.txt") as cfile:
    for line in cfile:
        caldwell = int(line[:3].strip(" "))
        if line[4] == "!":
            lineLst = line.split()
            ra = lineLst[2]
            de = lineLst[3]
            caldwells[caldwell] = RADEObj(ra,de)
        else:
            ngc = line[6:-1].replace(" ","")
            caldwell = int(line[:3].strip(" "))
            caldwells[caldwell] = ngcMap[ngc]
  
  return messiers, caldwells

def drawConstLines(baseMap,color="k",linewidth=1):
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
    drawLinesAroundBounderies(baseMap,starXs,starYs,color=color,linestyle='-',marker=None,linewidth=linewidth)

class ConstBoundaries(object):
  def __init__(self,localFn="data/bound_20.dat",url="https://cdsarc.unistra.fr/ftp/VI/49/bound_20.dat.gz"):
    infile = None
    try:
      infile = gzip.GzipFile(localFn)
    except:
      print(f"ConstBoundaries: Couldn't find local data file '{localFn}', attempting to download '{url}'...")
      infile = None
      try:
        r = requests.get(url)
        if r.status_code != 200:
          raise Exception(f"Couldn't download '{url}', status_code: {r.status_code}")
        infile = open(localFn, 'wb')
        infile.write(r.content)
        print("Done downloading %s to %s" % (url, localFn))
        infile.close()
        infile = gzip.GzipFile(localFn)
      except Exception as e:
        print(f"Error in ConstBoundaries: Couldn't find local data file and couldn't download. {e}")
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
        if cst in constBoundRaw:
            constBoundRaw[cst].append([RA,DE])
        else:
            constBoundRaw[cst] = [[RA,DE]]
    infile.close()
    self.constBoundRaw = constBoundRaw

  def draw(self,baseMap,color="k",linewidth=1):
    constBoundRaw = self.constBoundRaw
    for cst in constBoundRaw:
      mapXs = []
      mapYs = []
      for point in constBoundRaw[cst]:
        mapx,mapy = baseMap.project(point[0],point[1])
        mapXs.append(mapx)
        mapYs.append(mapy)
      #m.plot(mapXs,mapYs,"--c",alpha=0.7)
      drawLinesAroundBounderies(baseMap,mapXs,mapYs,color=color,linestyle='-',marker=None,linewidth=linewidth)

def readCatalog(localFn,url,processLineFunc):
  result = []
  infile = None
  try:
    infile = open(localFn)
  except:
    print(f"readCatalog: Couldn't find local data file '{localFn}', attempting to download '{url}'...")
    infile = None
    try:
      r = requests.get(url)
      if r.status_code != 200:
        raise Exception(f"Couldn't download '{url}', status_code: {r.status_code}, text: '{r.text}'")
      infile = open(localFn, 'wb')
      infile.write(r.content)
      print("Done downloading %s to %s" % (url, localFn))
      infile.close()
      infile = open(localFn)
    except Exception as e:
      print(f"Error in readCatalog: Couldn't find local data file and couldn't download b/c: {e}")
      print("Exiting.")
      sys.exit(1)
  for row in infile:
    result.append(processLineFunc(row))
  infile.close()
  return result

def readHip():
  return readCatalog("data/hip_main.dat","https://cdsarc.unistra.fr/ftp/I/239/hip_main.dat",HipEntry)

def readNGC():
  return readCatalog("data/ngc2000.dat","https://cdsarc.unistra.fr/ftp/VII/118/ngc2000.dat",NGCEntry)

def readNGCNames():
  return readCatalog("data/ngc2000_names.dat","https://cdsarc.unistra.fr/ftp/VII/118/names.dat",NGCNameEntry)

def readHCG():
  return readCatalog("data/hcg_groups.dat","https://cdsarc.unistra.fr/ftp/VII/213/groups.dat",HCGEntry)

def readCaldwell():
  return readCatalog("data/caldwell.txt","",CaldwellEntry)

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

    ngcNames = readNGCNames()
    self.messiers, self.caldwells = makeMessierDict(ngcObjs,ngcNames)

    self.hcgObjs = readHCG()

  def createMap(self,basemapArgs):
    if type(basemapArgs) != dict:
      raise Exception("StarMapper.createMap requires dict argument of Basemap args")
    basemapArgs['celestial'] = True
    m = Basemap(**basemapArgs)
    return m

  def drawStars(self,basemap):
    #self.dataArray=self.dataArray[:15000]
    mapx,mapy = basemap.project(self.dataArray[:,0],self.dataArray[:,1])
    mags = self.dataArray[:,2]
    maxSize = 300.0
    minSize = 5.0
    maxMag = max(mags)
    minMag = min(mags)
    sizes = (maxMag-mags)*(maxSize-minSize)/(maxMag-minMag)+minSize
    basemap.scatter(mapx,mapy,s=sizes,marker=".",facecolor="0.9",edgecolor="0.9",lw=0)

  def drawConsts(self,basemap):
    ConstBoundaries().draw(basemap,color="0.6")
    drawConstLines(basemap,color="k",linewidth=1.5)

  def drawGx(self,basemap,c='r'):
    mapx,mapy = basemap.project(self.ngcGx[:,0],self.ngcGx[:,1])
    basemap.plot(mapx,mapy,color=c,marker=".",linestyle="None",markeredgewidth=0.,markersize=9)
  def drawOC(self,basemap,c='g'):
    mapx,mapy = basemap.project(self.ngcOC[:,0],self.ngcOC[:,1])
    basemap.plot(mapx,mapy,color=c,marker=".",linestyle="None",markeredgewidth=0.,markersize=9)
  def drawGb(self,basemap,c='b'):
    mapx,mapy = basemap.project(self.ngcGb[:,0],self.ngcGb[:,1])
    basemap.plot(mapx,mapy,color=c,marker=".",linestyle="None",markeredgewidth=0.,markersize=9)
  def drawNb(self,basemap,c='m'):
    mapx,mapy = basemap.project(self.ngcNb[:,0],self.ngcNb[:,1])
    basemap.plot(mapx,mapy,color=c,marker=".",linestyle="None",markeredgewidth=0.,markersize=9)

  def drawMessiers(self,basemap,ax,c='b'):
    messier_keys = sorted(self.messiers)
    messiers = [self.messiers[x] for x in messier_keys]
    types = [x.getType() for x in messiers]
    types = []
    for x in messiers:
        try:
            types.append(x.getType())
        except AttributeError:
            types.append("")
    types = numpy.array(types)
    rade = [[x.getRA(),x.getDE()] for x in messiers]
    rade = numpy.array(rade)
    mapx,mapy = basemap.project(rade[:,0],rade[:,1])
    for t, sym, _ in ngcTypeSymbols:
        basemap.plot(mapx[types==t],mapy[types==t],color=c,marker=sym,linestyle="None",markeredgewidth=0.,markersize=9)
    for s, x, y in zip(messier_keys, mapx, mapy):
        ax.annotate("M"+str(s),(x,y),color=c,textcoords="offset points",xytext=(0,2),ha="center",va="bottom",fontsize=10)

  def drawCaldwells(self,basemap,ax,c='r'):
    caldwell_keys = sorted(self.caldwells)
    caldwells = [self.caldwells[x] for x in caldwell_keys]
    types = []
    for x in caldwells:
        try:
            types.append(x.getType())
        except AttributeError:
            types.append("")
    types = numpy.array(types)
    rade = [[x.getRA(),x.getDE()] for x in caldwells]
    rade = numpy.array(rade)
    mapx,mapy = basemap.project(rade[:,0],rade[:,1])
    for t, sym, _ in ngcTypeSymbols:
        basemap.plot(mapx[types==t],mapy[types==t],color=c,marker=sym,linestyle="None",markeredgewidth=0.,markersize=9)
    for s, x, y in zip(caldwell_keys, mapx, mapy):
        ax.annotate("C"+str(s),(x,y),color=c,textcoords="offset points",xytext=(0,2),ha="center",va="bottom",fontsize=10)

  def drawHCG(self,basemap,ax,c='g'):
    hcg_nums = [x.getHCG() for x in self.hcgObjs]
    rade = [[x.getRA(),x.getDE()] for x in self.hcgObjs]
    rade = numpy.array(rade)
    mapx,mapy = basemap.project(rade[:,0],rade[:,1])
    basemap.plot(mapx,mapy,color=c,marker='.',linestyle="None",markeredgewidth=0.,markersize=9)
    for s, x, y in zip(hcg_nums, mapx, mapy):
        ax.annotate("HCG"+str(s),(x,y),color=c,textcoords="offset points",xytext=(0,2),ha="center",va="bottom",fontsize=10)

  def drawGrid(self,basemap,label=True,color="0.85"):
    if isinstance(basemap,Basemap):
      labels=[False,False,False,False]
      if label:
        labels=[True,True,False,False]
      basemap.drawparallels(numpy.arange(-90.,91.,10.),
                            labels=labels,
                            dashes=[],
                            color = color#,
                            #linewidth = 1.
        )
      labels=[False,False,False,False]
      if label:
        labels=[False,False,True,True]
      basemap.drawmeridians(numpy.arange(-180.,181.,360/24.),
                            labels=labels,
                            fmt=lambda x: "{0:.0f}h".format(x/15.),
                            dashes=[],
                            color = color#,
                            #linewidth = 1.
        )

  def drawEcliptic(self,basemap,color="0.85",linestyle="-",marker=None):
    xData = numpy.linspace(-180,180,1000)
    yData = 23.44*numpy.sin(xData*numpy.pi/180.)
    xsToPlot, ysToPlot = basemap.project(xData,yData)
    basemap.plot(xsToPlot,ysToPlot,color=color,linestyle=linestyle,marker=marker)

class DeepSkyLegend:

  def __init__(self,fig):
    self.symbols = []
    for t, sym, title in ngcTypeSymbols:
        self.symbols.append(mlines.Line2D([],[],color="k",marker=sym,linestyle="None",markeredgewidth=0.,markersize=10,label=title))
    self.colors = []
    for c, title in [("b","Messier Objects"), ("r","Caldwell Objects"), ("g","Hickson Compact Groups (of Galaxies)")]:
        self.colors.append(mlines.Line2D([],[],color=c,marker="$●$",linestyle="None",markeredgewidth=0.,markersize=10,label=title))

    axloc_pos = [0.75,0.75,0.20,0.1]
    loc = "upper right"
    self.legend = fig.legend(handles=self.symbols,title="Deep Sky Object Types",loc=loc,bbox_to_anchor=axloc_pos,mode="expand")
    axloc_pos = [0.75,0.85,0.20,0.1]
    loc = "upper right"
    self.legend2 = fig.legend(handles=self.colors,title="Deep Sky Catalogs",loc=loc,bbox_to_anchor=axloc_pos,mode="expand")

if __name__ == "__main__":

  rcParams["font.size"] = 20.0
  rcParams["font.family"] = "Optima Nova LT Pro"
  rcParams["grid.linestyle"] = "-"
  rcParams["grid.linewidth"] = 1
  rcParams["grid.color"] = "0.85"

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
  axes = [axMain,axNP,axSP]
  for m,ax in zip(maps,axes):
    sm.drawGrid(m)
    sm.drawEcliptic(m)
    sm.drawConsts(m)
    cns.drawConstNames(m,fontsize="small",color="k",weight="bold")
    #sm.drawGb(m)
    #sm.drawGx(m)
    #sm.drawNb(m)
    #sm.drawOC(m)
    sm.drawMessiers(m,ax)
    sm.drawCaldwells(m,ax)
    sm.drawHCG(m,ax)
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

  for location,y in [("NM Skies",34.5),("Spain",38.),("Australia",-31.25)]:
    xyAxis = mMain(0,y)
    xyDisplay = dataToDisplay.transform(xyAxis)
    xyFig = displayToFigure.transform(xyDisplay)
    fig.text(0.04,xyFig[1],"{0}".format(location),ha="right",va="center",size="small")

  drawNGCPlots = True
  useRobinson = True
  if drawNGCPlots:
    axGb = fig.add_axes([0.75,0.03,0.2,0.08]) # left, bottom, width, height in fraction of fig
    axGx = fig.add_axes([0.75,0.13,0.2,0.08]) # left, bottom, width, height in fraction of fig
    axOC = fig.add_axes([0.05,0.03,0.2,0.08]) # left, bottom, width, height in fraction of fig
    axNb = fig.add_axes([0.05,0.13,0.2,0.08]) # left, bottom, width, height in fraction of fig

    fig.text(0.85,0.11,"NGC Globular Clusters",ha="center",va="bottom",size="large",fontweight="bold")
    fig.text(0.85,0.21,"NGC Galaxies",ha="center",va="bottom",size="large",fontweight="bold")
    fig.text(0.15,0.11,"NGC Open Clusters",ha="center",va="bottom",size="large",fontweight="bold")
    fig.text(0.15,0.21,"NGC Nebulae",ha="center",va="bottom",size="large",fontweight="bold")

    if useRobinson:
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
    else:
        mNb = sm.createMap({
          'projection':'cyl',
          'llcrnrlat':-80,
          'urcrnrlat':80,
          'llcrnrlon':-180,
          'urcrnrlon':180,
          'lat_ts':20,
          'ax':axNb,
        })
        mOC = sm.createMap({
          'projection':'cyl',
          'llcrnrlat':-80,
          'urcrnrlat':80,
          'llcrnrlon':-180,
          'urcrnrlon':180,
          'lat_ts':20,
          'ax':axOC,
        })
        mGb = sm.createMap({
          'projection':'cyl',
          'llcrnrlat':-80,
          'urcrnrlat':80,
          'llcrnrlon':-180,
          'urcrnrlon':180,
          'lat_ts':20,
          'ax':axGb,
        })
        mGx = sm.createMap({
          'projection':'cyl',
          'llcrnrlat':-80,
          'urcrnrlat':80,
          'llcrnrlon':-180,
          'urcrnrlon':180,
          'lat_ts':20,
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

  DeepSkyLegend(fig)
  fig.text(0.05,0.95,"Map of Constellations \n& Deep Sky Objects",ha="left",va="top",size=60,fontweight="bold")
  fig.text(0.05,0.90,"Justin Hugon",ha="left",va="top",size="large",fontweight="bold")

  fig.savefig('map.png')
  fig.savefig('map.pdf')
