#!/usr/local/bin/python
from ROOT import*
import sys,os
from array import array
import numpy

class WindowTuner:
    """
    The window tuning in prdonwstream has unit problems
    plot all the hits in an event, then we can plot the window
    """
    def __init__(self, x1hits,uhits,vhits,x2hits,TSeed):
        self.x1hits = x1hits
        self.uhits = uhits
        self.vhits = vhits
        self.x2hits = x2hits
        self.TSeed = TSeed
        #for x1, draw all the hits in x.
        self.x1HitPositions = []
        self.x1HitPositions_MC = []
        self.x1Zpositions = []
        self.x1Zpositions_MC = []

        self.x2HitPositions = []
        self.x2HitPositions_MC = []
        self.x2Zpositions = []
        self.x2Zpositions_MC = []

        for hit in x1hits:
            if hit.MatchedToSeed==1:
                self.x1HitPositions_MC.append(hit.x)
                self.x1Zpositions_MC.append(hit.zmid)
            else:
                self.x1HitPositions.append(hit.x)
                self.x1Zpositions.append(hit.zmid)
        for hit in x2hits:
            if hit.MatchedToSeed==1:
                self.x2HitPositions_MC.append(hit.x)
                self.x2Zpositions_MC.append(hit.zmid)
            else:
                self.x2HitPositions.append(hit.x)
                self.x2Zpositions.append(hit.zmid)
#now we have all the hits
        if len(self.x1Zpositions_MC)>0:
            self.x1hitsMCgraph = TGraph(len(self.x1HitPositions_MC),
                                        array('f',self.x1HitPositions_MC),
                                        array('f',self.x1Zpositions_MC))
            self.x1hitsMCgraph.SetMarkerStyle(22)
            self.x1hitsMCgraph.SetMarkerColor(kGreen+2)
            self.x1hitsMCgraph.SetMarkerSize(2);
            
        if len(self.x2Zpositions_MC)>0:
            self.x2hitsMCgraph = TGraph(len(self.x2HitPositions_MC),
                                        array('f',self.x2HitPositions_MC),
                                        array('f',self.x2Zpositions_MC))
            self.x2hitsMCgraph.SetMarkerStyle(22)
            self.x2hitsMCgraph.SetMarkerColor(kGreen+2)
            self.x2hitsMCgraph.SetMarkerSize(2);

        self.x1hitsOther = TGraph(len(self.x1HitPositions),
                                  array('f',self.x1HitPositions),
                                  array('f',self.x1Zpositions))

        self.x2hitsOther = TGraph(len(self.x2HitPositions),
                                  array('f',self.x2HitPositions),
                                  array('f',self.x2Zpositions))
        for othergr in [self.x1hitsOther,self.x2hitsOther]:
            othergr.SetMarkerStyle(20)
            othergr.SetMarkerColor(kRed+4)
            othergr.SetMarkerSize(2)
        can = TCanvas('can','can',1200,600)
        can.Divide(2,1)
        can.cd(1)
        #self.x1hitsOther.Draw('ap')
        self.mg = TMultiGraph()
        self.mg.SetTitle("hits in x;x[mm];z[mm]")
        self.mg.Add(self.x1hitsOther)
        self.mg.Add(self.x2hitsOther)

        if len(self.x1Zpositions_MC)>0:
            self.mg.Add(self.x1hitsMCgraph)
        if len(self.x2Zpositions_MC)>0:
            self.mg.Add(self.x2hitsMCgraph)
        self.mg.Draw('ap')
        #hits are now drawn, so now draw the selection window.
        can.Update()
        raw_input('continue?')
        
