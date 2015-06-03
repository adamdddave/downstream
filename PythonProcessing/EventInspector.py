#!/usr/local/bin/python
from ROOT import*
import sys,os,re
from array import array
from HitInformation import *
from graphs2display import*
from FitComparer import *

class EventInspector:

    def __init__(self,g2d,fc,p1,p2,p3,p4):
        g2d.set_tseed_projections()
        g2d.make_x_y_z_projections()
        curr_yplot = TF1("curr_xplot","[0]+[1]*x",g2d.mg.GetXaxis().GetXmin(),g2d.mg.GetXaxis().GetXmax())
        curr_xplot = TF1("curr_yplot","[0]+[1]*x",g2d.mg2.GetXaxis().GetXmin(),g2d.mg2.GetXaxis().GetXmax()) 
        curr_xplot.SetParameters(p1,p2)
        curr_yplot.SetParameters(p3,p4)
        for plot in [curr_xplot,curr_yplot]:
            plot.SetLineColor(kGreen+2)
            plot.SetLineStyle(kDashed)
            
        g2d.c1.cd(1)
        curr_yplot.Draw("lsame")
        g2d.c1.Update()
        g2d.c1.cd(2)
        curr_xplot.Draw("lsame")
        g2d.c1.cd(3)
        g2d.leg.AddEntry(curr_xplot,"Fit from New #chi^{2}","l")
        gPad.cd(1)
        g2d.leg.Draw()
        #fc = FitComparer(chi2)
        print 'MC Particle P = %f MeV, Downstream P = %f MeV '%(g2d.mcParticleP,g2d.Tseed.dsP)
        g2d.c1.Update()
        fc.print_table()
