#!/usr/bin/env python2

import json
import csv
import numpy

class ConstNames(object):

  def __init__(self):
    self.boxes = None
    with open("data/constBoxes.json") as rawfile:
      self.boxes = json.load(rawfile)
    #print self.boxes

    self.abbrevs = {}
    with open("data/constAbbrevs.csv") as abbsCSV:
      abbsCSVReader = csv.reader(abbsCSV)
      for line in abbsCSVReader:
        self.abbrevs[line[0]] = line[1]
    #print self.abbrevs

    self.makeBoxCenters()
    self.makeCentroids()
      
  def makeBoxCenters(self):
    self.centers = {}
    for const in self.boxes:
      centers = []
      for box in self.boxes[const]:
        #print const,box
        raCenter = (box[0][0] + box[1][0])*0.5
        deCenter = (box[0][1] + box[1][1])*0.5
        #print raCenter, deCenter
        centers.append([raCenter,deCenter])
      self.centers[const] = centers
      #print const,self.centers[const]

  def makeCentroids(self):
    self.centroids = {}
    for const in self.centers:
      ras = []
      des = []
      ra0 = self.centers[const][0][0]
      for center in self.centers[const]:
        ra = center[0]
        de = center[1]
        ra -= ra0
        if ra > 0:
          if ra > 180:
            ra -= 360.
        else:
          if ra < -180:
            ra += 360.
        ras.append(ra)
        des.append(de)
      #print ras
      raCentroid = numpy.mean(ras)
      raCentroid += ra0
      if raCentroid < 0:
        raCentroid += 360.*2
      raCentroid = raCentroid % 360.
      #print ra0, raCentroid
      deCentroid = numpy.mean(des)
      self.centroids[const] = (raCentroid,deCentroid)

  def drawConstNames(self,ax):
    for const in self.centroids:
      ra, de = self.centroids[const]
      if ra > 180.:
        ra -= 360.
      if ax.rho_lim:
        if de > ax.rho_lim[1] or de < ax.rho_lim[0]:
          continue
      else:
        if de > ax.urcrnrlat or de < ax.llcrnrlat:
          continue
      #print const, ra, de
      ra, de = ax.project(ra,de)
      try:
        ax.text(ra,de,self.abbrevs[const],fontsize=16,va="center",ha="center")
      except:
        ax.ax.text(ra,de,self.abbrevs[const],fontsize=15,va="center",ha="center")
    

if __name__ == "__main__":
  cn = ConstNames()
  
