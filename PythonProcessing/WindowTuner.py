#!/usr/local/bin/python
from ROOT import*
import sys,os
from os import listdir
from os.path import isfile, join
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
                print 'setting nonmatched x position to %f'%(hit.x)
                self.x1HitPositions.append(hit.x)
                self.x1Zpositions.append(hit.zmid)
        for hit in x2hits:
            if hit.MatchedToSeed==1:
                self.x2HitPositions_MC.append(hit.x)
                self.x2Zpositions_MC.append(hit.zmid)
            else:
                self.x2HitPositions.append(hit.x)
                self.x2Zpositions.append(hit.zmid)
    def makeGraphs(self):
#now we have all the hits
        if len(self.x1Zpositions_MC)>0:
            print'setting x1 mc z positions to ',self.x1Zpositions_MC
            print 'setting x1 mc x positions to ',self.x1HitPositions_MC
            self.x1hitsMCgraph = TGraph(len(self.x1HitPositions_MC),
                                        array('f',self.x1Zpositions_MC),
                                        array('f',self.x1HitPositions_MC))
            self.x1hitsMCgraph.SetMarkerStyle(22)
            self.x1hitsMCgraph.SetMarkerColor(kGreen+2)
            self.x1hitsMCgraph.SetMarkerSize(2);
            
        if len(self.x2Zpositions_MC)>0:
            print'setting x2 mc z positions to ',self.x2Zpositions_MC
            print 'setting x2 mc x positions to ',self.x2HitPositions_MC
            self.x2hitsMCgraph = TGraph(len(self.x2HitPositions_MC),
                                        array('f',self.x2Zpositions_MC),
                                        array('f',self.x2HitPositions_MC))
            self.x2hitsMCgraph.SetMarkerStyle(22)
            self.x2hitsMCgraph.SetMarkerColor(kGreen+2)
            self.x2hitsMCgraph.SetMarkerSize(2);
        if len(self.x1HitPositions)>0:
            self.x1hitsOther = TGraph(len(self.x1HitPositions),
                                      array('f',self.x1Zpositions),
                                      array('f',self.x1HitPositions))
        if len(self.x2HitPositions)>0:
            self.x2hitsOther = TGraph(len(self.x2HitPositions),
                                      array('f',self.x2Zpositions),
                                      array('f',self.x2HitPositions))
        if len(self.x1HitPositions)>0 and len(self.x2HitPositions)>0:
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
        if len(self.x1Zpositions)>0:
            self.mg.Add(self.x1hitsOther)
        if len(self.x2Zpositions)>0:
            self.mg.Add(self.x2hitsOther)

        if len(self.x1Zpositions_MC)>0:
            self.mg.Add(self.x1hitsMCgraph)
        if len(self.x2Zpositions_MC)>0:
            self.mg.Add(self.x2hitsMCgraph)
        self.mg.Draw('ap')
        #hits are now drawn, so now draw the selection window.
        #x line looks like x(z) = xmag + (z-zmag)tx + curv (z-ztt)^2
        xfun = TF1('xfun','[0]+[1]*(x-[2])+[3]*(x-[4])*(x-[4])',2320,2700)
        xfunplus = TF1('xfunplus','[0]+[1]*(x-[2])+[3]*(x-[4])*(x-[4])+[5]',2320,2700)
        xfunminus = TF1('xfunminus','[0]+[1]*(x-[2])+[3]*(x-[4])*(x-[4])-[5]',2320,2700)
        xfunplus2 = TF1('xfunplus2','[0]+[1]*(x-[2])+[3]*(x-[4])*(x-[4])+[5]',2320,2700)
        xfunminus2 = TF1('xfunminus2','[0]+[1]*(x-[2])+[3]*(x-[4])*(x-[4])-[5]',2320,2700)

        for fun in [xfun,xfunplus,xfunminus,xfunplus2,xfunminus2]:
            fun.SetParameter(0,self.TSeed.preselXmag)
            fun.SetParameter(1,self.TSeed.preselTx)
            fun.SetParameter(2,self.TSeed.preselZmag)
            fun.SetParameter(3,1.7e-5*(self.TSeed.tx - self.TSeed.preselTx))#curvature from prdowntrack
            fun.SetParameter(4,2485.)
        tolerance = 1.5+ 20000./self.TSeed.dsP
        if self.TSeed.dsP/20000. < 1./(9.5):
            print 'using max tolerance instead'
            tolerance = 10
        xfunplus.SetParameter(5,tolerance)
        xfunminus.SetParameter(5,tolerance)
        xfunplus2.SetParameter(5,2*tolerance)
        xfunminus2.SetParameter(5,2*tolerance)
        xfun.Draw("same")
        for fun in [xfunplus,xfunminus]:
            fun.SetLineStyle(kDashed)
            fun.Draw("same")
        for fun in [xfunplus2,xfunminus2]:
            fun.SetLineStyle(kDashed)
            fun.SetLineColor(kBlue)
            fun.Draw("same")
        can.Update()
        self.pt = TPaveText(0.1,0.5,0.9,0.9)
        self.pt.AddText(('Downstream momentum = %f'%(self.TSeed.dsP)))
        self.pt.AddText(('total number of hits: x1(%d),x2(%d)'%(len(self.x1hits),len(self.x2hits))))
        self.pt.AddText(('MC matched: x1(%d),x2(%d)'%(len(self.x1HitPositions_MC),len(self.x2HitPositions_MC))))
        can.cd(2)
        self.pt.Draw()
        can.Update()
        sav = raw_input('save? [y]/n')
        #automagically generate the save name
        if not 'n' in sav:
            mypath='./search_window_in_x'
            onlyfiles = [ f for f in listdir(mypath) if isfile(join(mypath,f)) ]
            #print onlyfiles
            if len(onlyfiles)<1:
                can.SaveAs("./search_window_in_x/window_tune_ev1.pdf")
            else:
                thelist = []
                for file in onlyfiles:
                    fname=(((file.split('.'))[0]).split('ev'))[1]
                    thelist.append(int(fname))
                #print thelist
                thelist.sort()
                can.SaveAs("./search_window_in_x/window_tune_ev%d.pdf"%(thelist[-1]+1))
            #can.SaveAs
        bs=raw_input('continue? [y]/n')
        if 'n' in bs:
            exit()
