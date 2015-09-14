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
      for center in self.centers[const]:
        ras.append(center[0])
        des.append(center[1])
      raCentroid = numpy.mean(ras)
      deCentroid = numpy.mean(des)
      self.centroids[const] = [raCentroid,deCentroid]

  def drawConstNames(self):
    pass

if __name__ == "__main__":
  cn = ConstNames()
  
