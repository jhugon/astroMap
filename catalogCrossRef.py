#!/usr/bin/env python2

import urllib2
import sys
import sqlalchemy
from sqlalchemy import Column, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

Base = declarative_base()

class CrossRefEntry(Base):

  __tablename__ = "crossref"

  HD = Column(Integer,primary_key=True)
  DM = Column(String(12))
  GC = Column(Integer)
  HR = Column(Integer)
  HIP = Column(Integer)
  RAh = Column(Integer)
  RAm = Column(Integer)
  RAs = Column(Float())
  DEsign = Column(String(1))
  DEd = Column(Integer)
  DEm = Column(Integer)
  DEs = Column(Float())
  Vmag = Column(Float())
  Fl = Column(Integer)
  Bayer = Column(String(5))
  Cst = Column(String(3))

  def __repr__(self):
    result = "<CrossRefEntry("
    for i in ["HD","DM","GC","HR","HIP","Fl","Bayer","Cst"]:
        result += i+"="+str(getattr(self,i))+", "
    result += "RA={0:.1f}, ".format(self.getRAd())
    result += "DE={0:.1f}".format(self.getDEd())
    result += ")>"
    return result

  def getRAhms(self):
    return self.RAh,self.RAm,self.RAs
  def getDEdms(self):
    return self.DEsign,self.DEd,self.DEm,self.DEs
  def getRAd(self):
    inHours = float(self.RAh)+float(self.RAm)/60.+self.RAs/3600.
    return inHours*15.
  def getDEd(self):
    result = float(self.DEd)+float(self.DEm)/60.+float(self.DEs)/3600.
    if self.DEsign == "-":
      result *= -1.
    return result
  def getHD(self): return self.HD
  def getDM(self): return self.DM
  def getGC(self): return self.GC
  def getHR(self): return self.HR
  def getHIP(self): return self.HIP
  def getVmag(self): return self.Vmag
  def getFl(self): return self.Fl
  def getBayer(self): return self.Bayer
  def getCst(self): return self.Cst

class CatalogCrossRef(object):
  def __init__(self,localFn="data/crossRefCat.dat",url="ftp://cdsarc.u-strasbg.fr/pub/cats/IV/27A/catalog.dat"):
    self.localFile = None
    try:
      self.localFile = open(localFn)
    except:
      print "CatalogCrossRef: Couldn't find local data file, attempting to download..."
      self.localFile = None
      try:
        self.download(localFn,url)
        self.localFile = open(localFn)
      except:
        print "Error in CatalogCrossRef: Couldn't find local data file and couldn't download"
        sys.exit(1)

    ## Now DB stuff
    #engine = sqlalchemy.create_engine('sqlite:///:memory:',echo=True)
    engine = sqlalchemy.create_engine('sqlite:///:memory:')
    self.engine = engine
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    self.Session = Session
    session = Session()

    for line in self.localFile:
      argDict = {
        "HD" : line[0:6],
        "DM" : line[7:19],
        "GC" : line[20:25],
        "HR" : line[26:30],
        "HIP" : line[31:37],
        "RAh" : line[38:40],
        "RAm" : line[40:42],
        "RAs" : line[42:47],
        "DEsign" : line[48],
        "DEd" : line[49:51],
        "DEm" : line[51:53],
        "DEs" : line[53:57],
        "Vmag" : line[58:63],
        "Fl" : line[64:67],
        "Bayer" : line[68:73],
        "Cst" : line[74:77],
      }
      for i in ["HD","GC","HR","HIP","RAh","RAm","DEd","DEm","Fl"]:
        try:
          argDict[i] = int(argDict[i])
        except:
          argDict[i] = None
      for i in ["RAs","DEs","Vmag"]:
        try:
          argDict[i] = float(argDict[i])
        except:
          argDict[i] = None
      for i in ["Cst","Bayer","DM"]:
        argDict[i] = argDict[i].upper().strip()
        if argDict[i] == "":
          argDict[i] = None
      tmpCre = CrossRefEntry(
          **argDict
        )

      session.add(tmpCre)
    self.localFile.close()
    session.flush()

  def download(self,localFn,url):
    u = urllib2.urlopen(url)
    f = open(localFn, 'wb')
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
        f.write(buffer)
        status = r"%10d  [%3.2f%%]" % (file_size_dl, file_size_dl * 100. / file_size)
        status = status + chr(8)*(len(status)+1)
        print status,
    f.close()
    print "Done downloading %s to %s" % (url, localFn)

  def findByHD(self,HD):
    session = self.Session()
    result = session.query(CrossRefEntry).filter_by(HD=HD).first()
    return result
  def findByHIP(self,HIP):
    session = self.Session()
    result = session.query(CrossRefEntry).filter_by(HIP=HIP).first()
    return result
  def findByFl(self,Cst,Fl):
    session = self.Session()
    result = session.query(CrossRefEntry).filter_by(Cst=Cst.upper(),Fl=Fl).first()
    return result
  def findByBayer(self,Cst,Bayer):
    session = self.Session()
    result = session.query(CrossRefEntry).filter_by(Cst=Cst.upper(),Bayer=Bayer.upper()).first()
    return result

if __name__ == "__main__":
  ccr = CatalogCrossRef()
  print "CAS 45:",ccr.findByFl("CAS",45)
  print "VIR 3:",ccr.findByFl("VIR",3)
  print "HD 20010:",ccr.findByHD(20010)
  print "HD 175813:",ccr.findByHD(175813)
  print "Should be same:",ccr.findByHD(20010),ccr.findByFl("FOR",1),ccr.findByBayer("FOR","ALF")
  print "Should be same:",ccr.findByHD(36597),ccr.findByFl("Col",1),ccr.findByBayer("Col","ALF")
  print "Should be same:",ccr.findByHD(43834),ccr.findByFl("Men",1),ccr.findByBayer("Men","Alf")
  print
  print "For Vir 3:"
  vir3 = ccr.findByFl("VIR",3)
  print " RA:",vir3.getRAhms()," DE:",vir3.getDEdms()
  print " RA:",vir3.getRAd()," DE:",vir3.getDEd()
  print " HD: ",vir3.getHD(), "Vmag: ",vir3.getVmag()
  print "For Vir 3:"
  hd352 = ccr.findByHD(352)
  print " RA:",hd352.getRAhms()," DE:",hd352.getDEdms()
  print " RA:",hd352.getRAd()," DE:",hd352.getDEd()
  print "  HIP:",hd352.getHIP(), "Vmag: ",hd352.getVmag()
