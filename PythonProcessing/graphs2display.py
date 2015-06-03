#!/usr/local/bin/python
from ROOT import *
from array import array
from HitInformation import *
from math import *
import numpy, scipy
class graphs2display:
    """
    this is the class which displays all of the graphs we would need
    the classes needed are Hit information, which has all the info from the hit
    and that's it.
    Remeber, the first thing that the v plane will consider is xMax, then xmin, since it is a line like
    .
     .
      .
       .
        .
    Otherwise things don't fall on the right line.
    """
    def __init__(self,hitx1,hitu,hitv,hitx2):
        self.x1=hitx1
        self.x2=hitx2
        self.u=hitu
        self.v=hitv
        #self.mcParticleP = particleP
    def set_counter_pave(self,counter):
        self.pt=TPaveText(0.1,0.1,0.8,0.9)
        self.pt.SetFillColor(0)
        self.pt.SetBorderSize(0)
        self.pt.AddText("Event %i\n"%counter)
        self.evCounter = counter
#         self.pt.AddText(('True P from x1 %f')%( (self.x1.trueHitP)))
#         self.pt.AddText(('True P from u  %f')%( (self.u.trueHitP)))
#         self.pt.AddText(('True P from v  %f')%( (self.v.trueHitP)))
#         self.pt.AddText(('True P from x2 %f')%( (self.x2.trueHitP)))
#         self.pt.AddText(('True P from MC particle %f')%( (self.mcParticleP) ))
#         self.pt.AddText(('True x1 tx %f')% (self.x1.trueTx))
#         self.pt.AddText(('True u tx %f')% (self.u.trueTx))
#         self.pt.AddText(('True v tx %f')% (self.v.trueTx))
#         self.pt.AddText(('True x2 tx %f')% (self.x2.trueTx))
#         self.pt2 =TPaveText(0.1,0.0,0.8,0.1)

#         self.pt2.SetFillColor(0)
#         self.pt2.SetBorderSize(0)
#         #self.pt2.AddText((''))
#         self.pt2.AddText(('Equation for tseed yproj:y0=y_{Tseed}-(ty_{Tseed}z_{Tseed})'))
#         self.pt2.AddText(('y_{FromSeed,At UT plane}=y0+(ty_{Tseed}z_{Plane})'))
        

    def setTseedInfo(self,tinfo):
        self.Tseed = tinfo
        self.TseedX=tinfo.x
        self.TseedY=tinfo.y
        self.TseedZ=tinfo.z
        self.TseedTx=tinfo.tx
        self.TseedTy=tinfo.ty
        self.TseedY0 = tinfo.y0
        self.TseedAtX1Y=tinfo.y0 + (tinfo.ty* self.x1.zHit)
        self.TseedAtUY=tinfo.y0 + (tinfo.ty* self.u.zHit)
        self.TseedAtVY=tinfo.y0 + (tinfo.ty* self.v.zHit)
        self.TseedAtX2Y=tinfo.y0 + (tinfo.ty* self.x2.zHit)
        #now some fanciness to try to get the U position
        #first, get it from the truth
        angle=(5.*pi/180.)
        upos= (self.u.trueXmid * cos(angle))+(self.u.trueYmid * sin(angle))
        #print 'uPos from truth = %f'%upos
        self.TseedAtUX=(upos - (self.TseedAtUY * sin(angle)))/cos(angle)
        #print 'tseed u,x projection = %f'%self.TseedAtUX
        #now try v
        vpos=(self.v.trueXmid* cos(-angle))+(self.v.trueYmid * sin(-angle))
        #print 'vpos from truth = %f'%vpos
        self.TseedAtVX=(vpos-(self.TseedAtVY*sin(-angle)))/cos(-angle)
        #print 'tseed v, xprojection = %f'%self.TseedAtVX
        #since we have the truth, now we can try to do it from the actual params
        #we don't expect these to really be on the right line
        #but let's project for xmin and xmax to see if we can possibly
        #correct the position for it
        
        #simply a copy from above
        #central value for the prediction
        upos2= (self.u.xMin * cos(angle))+(self.u.yMin * sin(angle))
        #print 'uPos2 from hit = %f'%upos
        self.TseedAtUX2=(upos2 - (self.TseedAtUY * sin(angle)))/cos(angle)
        #print 'tseed u from hit,x projection = %f'%self.TseedAtUX2
        #now try v
        vpos2=(self.v.xMin* cos(-angle))+(self.v.yMin * sin(-angle))
        #print 'vpos2 from hit = %f'%vpos2
        self.TseedAtVX2=(vpos2-(self.TseedAtVY*sin(-angle)))/cos(-angle)
        #print 'tseed v from hit, xprojection = %f'%self.TseedAtVX2
        #now using x at global y as in the line hit method
        self.UxAtyLineHit = self.u.xAtYEq0+(self.u.dxDy*self.TseedAtUY)
        self.VxAtyLineHit = self.v.xAtYEq0+(self.v.dxDy*self.TseedAtVY)
        #print 'u xposition from linehit method = %f'%self.UxAtyLineHit
        #print 'v xposition from linehit method = %f'%self.VxAtyLineHit
    
    def set_tseed_projections(self):
        y_pos_at_zpreut = self.TseedY0 + (self.Tseed.ty * 2320) #in case we want to add it to the graph.
        alphaY = 11.87
        betaY = (-6.4035e4)
        gammaY = (4.8315e8)
        
        alphaTY = 0.0043405
        betaTY = -28.455
        gammaTY = 1.7005e5
        rho = -0.766
        #rho =0 #check to see if this is the problem in the chi2
        #errorY = alphaY+betaY/self.hitx1.trueHitP + gammaY/(self.hitx1.trueHitP*self.hitx1.trueHitP)
        #errorTY = alphaTY+betaTY/self.hitx1.trueHitP + gammaTY/(self.hitx1.trueHitP*self.hitx1.trueHitP)
        #errorTYY = rho*(errorY*errorY*errorTY*errorTY)**0.5
        
        errorY = alphaY+betaY/abs(self.Tseed.dsP) + gammaY/(self.Tseed.dsP*self.Tseed.dsP)
        errorTY = alphaTY+betaTY/abs(self.Tseed.dsP) + gammaTY/(self.Tseed.dsP*self.Tseed.dsP)
        errorTYY = rho*(errorY*errorY*errorTY*errorTY)**0.5
        errorMatrix = numpy.array([[errorY*errorY,errorTYY],[errorTYY,errorTY*errorTY]])
        self.ypreUTgraph = TGraphAsymmErrors(1,
                                             array('f',[y_pos_at_zpreut]),
                                             array('f',[2320.]),
                                             array('f',[errorY]),
                                             array('f',[errorY]),
                                             array('f',[0]),
                                             array('f',[0]))

        self.ypreUTtruegraph = TGraphAsymmErrors(1,
                                                 array('f',[self.x1.trueEntryY + self.x1.trueTy*(2320. - self.x1.trueEntryZ)]),
                                                 array('f',[2320.]),
                                                 array('f',[0]),
                                                 array('f',[0]),
                                                 array('f',[0]),
                                                 array('f',[0]))
        ty_beta = 1.869
        y_beta = -3678
        tyTmeas = self.Tseed.ty
        tyTmeas+=self.Tseed.dsQ/self.Tseed.dsP*ty_beta #+ self.TSeed.dsQ*(-0.003)
        self.typreUTmomCor = TGraphErrors(1,
                                          array('f',[tyTmeas]),
                                          array('f',[2320]),
                                          array('f',[errorTY]),
                                          array('f',[0]))
        self.tyPreUTtruth = TGraph(1,
                                   array('f',[self.x1.trueTy]),
                                   array('f',[2320]))
        
        yat0 = self.Tseed.y-(self.Tseed.ty*self.Tseed.z)
        yat0_take2 = self.Tseed.y-(tyTmeas*self.Tseed.z)
        ymeas = yat0+self.Tseed.ty*2320.
        ymeas_take2 = yat0_take2 + tyTmeas*2320.
        ymeas+= y_beta*self.Tseed.dsQ/self.Tseed.dsP
        self.ypreUTMomAdjusted = TGraphAsymmErrors(1,
                                                   array('f',[ymeas]),
                                                   array('f',[2320]),
                                                   array('f',[errorY]),
                                                   array('f',[errorY]),
                                                   array('f',[0]),
                                                   array('f',[0]))
        
        self.ypreUTMomAdjustedWithSeedAdj = TGraphAsymmErrors(1,
                                                              array('f',[ymeas_take2]),
                                                              array('f',[2320]),
                                                              array('f',[errorY]),
                                                              array('f',[errorY]),
                                                              array('f',[0]),
                                                              array('f',[0]))
                                                   


    def make_x_y_z_projections(self):
        
        self.zvsygraph=TGraphAsymmErrors(4,#n
                                    array('f',[self.x1.y,self.u.y,self.v.y,self.x2.y]),#xaxis
                                    array('f',[self.x1.zHit,self.u.zHit,self.v.zHit,self.x2.zHit]),#yaxis
                                    array('f',[abs(self.x1.y-self.x1.yMin),abs(self.u.y-self.u.yMin),abs(self.v.y-self.v.yMin),abs(self.x2.y-self.x2.yMin)]),#exlow
                                    array('f',[abs(self.x1.y-self.x1.yMax),abs(self.u.y-self.u.yMax),abs(self.v.y-self.v.yMax),abs(self.x2.y-self.x2.yMax)]),#exhigh
                                    array('f',[0,0,0,0]),#eylow
                                    array('f',[0,0,0,0])#eyhigh
                                    )
        self.zvsygraph.SetTitle("hit z vs hit y; y[mm]; z[mm]")
        self.zvsxgraph=TGraphAsymmErrors(4,#n
                            array('f',[self.x1.x,self.u.x,self.v.x,self.x2.x]),#xaxis
                            array('f',[self.x1.zHit,self.u.zHit,self.v.zHit,self.x2.zHit]),#yaxis
                            array('f',[abs(self.x1.x-self.x1.xMin),abs(self.u.x-self.u.xMin),abs(self.v.x-self.v.xMin),abs(self.x2.x-self.x2.xMin)]),#exlow
                            array('f',[abs(self.x1.x-self.x1.xMax),abs(self.u.x-self.u.xMax),abs(self.v.x-self.v.xMax),abs(self.x2.x-self.x2.xMax)]),#exhigh
                            array('f',[0,0,0,0]),#eylow
                            array('f',[0,0,0,0])#eyhigh
                            )
        self.zvsxgraph.SetTitle("hit z vs hit x; x[mm]; z[mm]")

        
        #true shiz
        self.truezvsygraph=TGraphAsymmErrors(4,#n
                                array('f',[self.x1.trueYmid,self.u.trueYmid,self.v.trueYmid,self.x2.trueYmid]),#xaxis
                                array('f',[self.x1.trueZmid,self.u.trueZmid,self.v.trueZmid,self.x2.trueZmid]),#yaxis
                                array('f',[0,0,0,0]),#exlow
                                array('f',[0,0,0,0]),#exhigh
                                array('f',[0,0,0,0]),#eylow
                                array('f',[0,0,0,0])#eyhigh
                                )
        self.truezvsygraph.SetTitle("hit z vs hit y; y[mm]; z[mm]")

        self.truezvsxgraph=TGraphAsymmErrors(4,#n
                                             array('f',[self.x1.trueXmid,self.u.trueXmid,self.v.trueXmid,self.x2.trueXmid]),#xaxis
                                             array('f',[self.x1.trueZmid,self.u.trueZmid,self.v.trueZmid,self.x2.trueZmid]),#yaxis
                                             array('f',[0,0,0,0]),#exlow
                                             array('f',[0,0,0,0]),#exhigh
                                             array('f',[0,0,0,0]),#eylow
                                             array('f',[0,0,0,0])#eyhigh
                            )
        self.truezvsxgraph.SetTitle("hit z vs hit x; x[mm]; z[mm]")
        #now draw mike's projected y position
        self.tseedarrayY=TGraphAsymmErrors(4,#n
                                           array('f',[self.TseedAtX1Y,self.TseedAtUY,self.TseedAtVY,self.TseedAtX2Y]),#xaxis
                                           array('f',[self.x1.trueZmid,self.u.trueZmid,self.v.trueZmid,self.x2.trueZmid]),#yaxis
                                           array('f',[0,0,0,0]),#exlow
                                           array('f',[0,0,0,0]),#exhigh
                                           array('f',[0,0,0,0]),#eylow
                                           array('f',[0,0,0,0]))#eyhigh

        #x position from y from mike
        self.tseedarrayX=TGraphAsymmErrors(2,#n
                                           array('f',[self.TseedAtUX,self.TseedAtVX]),#xaxis
                                           array('f',[self.u.trueZmid,self.v.trueZmid]),#yaxis
                                           array('f',[0,0]),#exlow
                                           array('f',[0,0]),#exhigh
                                           array('f',[0,0]),#eylow
                                           array('f',[0,0]))#eyhigh
        #now the ones not from the truth
        #x position from y from mike
        self.tseedarrayX2=TGraphAsymmErrors(2,#n
                                           array('f',[self.TseedAtUX2,self.TseedAtVX2]),#xaxis
                                           array('f',[self.u.trueZmid,self.v.trueZmid]),#yaxis
                                           array('f',[0,0]),#exlow
                                           array('f',[0,0]),#exhigh
                                           array('f',[0,0]),#eylow
                                           array('f',[0,0]))#eyhigh
        #finally, the positions from x at y as in the line hit method
        self.tseedarrayUVLineHit = TGraphAsymmErrors(2,#n
                                                    array('f',[self.UxAtyLineHit,self.VxAtyLineHit]),#xaxis
                                                    array('f',[self.u.trueZmid,self.v.trueZmid]),#yaxis
                                                    array('f',[0,0]),#exlow
                                                    array('f',[0,0]),#eylow
                                                    array('f',[0,0]))#eyhigh
        #add ypreUT with momentum corrections, etc.
        
        self.tseedarrayX2.SetMarkerColor(kMagenta+2)
        self.tseedarrayX2.SetMarkerStyle(29)#filled star
        self.tseedarrayX2.SetMarkerSize(2.5)
        self.tseedarrayUVLineHit.SetMarkerColor(kCyan+2)
        self.tseedarrayUVLineHit.SetMarkerStyle(24)#open circle
        self.tseedarrayUVLineHit.SetMarkerSize(2.5)
        self.c1=TCanvas('c1','c1',1200,600)
        self.c1.Divide(3,1)
        print 'x1 z position = %f'%(self.x1.trueZmid)
        thegraphs=[self.zvsygraph,self.zvsxgraph,self.ypreUTgraph]
        trueGraphs=[self.truezvsygraph,self.truezvsxgraph,self.ypreUTtruegraph,self.tyPreUTtruth]
        projGraphs=[self.tseedarrayY,self.tseedarrayX]
        for graph in thegraphs:
            graph.SetMarkerSize(2)
            graph.SetMarkerStyle(3)#was 34
            graph.SetMarkerColor(kBlue)
            graph.GetYaxis().SetTitleOffset(1.5)
        for graph in trueGraphs:
            graph.SetMarkerSize(2.5)
            graph.SetMarkerStyle(kOpenDiamond)
            graph.SetMarkerColor(kGreen+2)
        for graph in projGraphs:
            graph.SetMarkerSize(2.5)
            graph.SetMarkerColor(kRed)
            graph.SetMarkerStyle(kMultiply)
            
        #now all the graphs are displayed, we really want to set a minimum and maximum for x and y
        #so that everything is shown correctly
        self.ypreUTgraph.SetLineColor(kBlue)
        self.ypreUTtruegraph.SetFillColor(kGreen+2)        
        self.ypreUTtruegraph.SetFillStyle(3144)
        self.ypreUTtruegraph.SetMarkerStyle(32)
        self.ypreUTMomAdjusted.SetLineColor(kMagenta+2)
        self.ypreUTMomAdjusted.SetMarkerColor(kMagenta+2)
        self.ypreUTMomAdjusted.SetMarkerSize(2.5)
        self.ypreUTMomAdjusted.SetMarkerStyle(30)
        self.ypreUTMomAdjustedWithSeedAdj.SetFillColor(kMagenta+2)
        self.ypreUTMomAdjustedWithSeedAdj.SetLineColor(kMagenta+2)
        self.ypreUTMomAdjustedWithSeedAdj.SetMarkerSize(2.5)
        self.ypreUTMomAdjustedWithSeedAdj.SetMarkerStyle(10)
        self.c1.cd(1)
        self.mg = TMultiGraph()
        self.mg.SetTitle("hit z vs hit y; y[mm]; z[mm]")
        # self.zvsygraph.Draw('ap')
        #self.truezvsygraph.Draw('p same')
        #self.tseedarrayY.Draw('p same')
        self.mg.Add(self.zvsygraph)
        self.mg.Add(self.truezvsygraph)
        self.mg.Add(self.ypreUTgraph)
        self.mg.Add(self.ypreUTtruegraph)
        self.mg.Add(self.ypreUTMomAdjusted)
        self.mg.Add(self.tseedarrayY)
        #self.mg.Add(self.ypreUTMomAdjustedWithSeedAdj)
        self.mg.Draw('ap')
        self.mg.GetYaxis().SetTitleOffset(1.5)
        #print "fit z vs y"
        self.truezvsygraph.Fit("pol1",'Q')
        fitted_vals = self.truezvsygraph.GetFunction("pol1")
        #print "fitted values are"
        #print fitted_vals.GetParameter(0),fitted_vals.GetParameter(1),
        #print "**"
        self.trueYZfittedYintercept = fitted_vals.GetParameter(0)
        self.trueYZfittedYslope = fitted_vals.GetParameter(1)
        gStyle.SetOptFit(111)
        
        self.c1.Update()
        self.c1.cd(2)
        #self.zvsxgraph.GetXaxis().SetRangeUser(displayXmin,displayXmax)
        #self.zvsxgraph.Draw('ap')
        #self.truezvsxgraph.Draw('p same')
        #self.tseedarrayX.Draw('p same')
        #self.tseedarrayX2.Draw('p same')
        self.mg2 = TMultiGraph()
        self.mg2.SetTitle("hit z vs hit x; x[mm]; z[mm]")
        self.mg2.Add(self.zvsxgraph)
        self.mg2.Add(self.truezvsxgraph)
        #self.mg2.Add(self.tseedarrayX)
        #self.mg2.Add(self.tseedarrayX2)
        #self.mg2.Add(self.tseedarrayUVLineHit)
        self.mg2.Draw('ap')
        self.mg2.GetYaxis().SetTitleOffset(1.5)
        #print "fit z vs x"
        self.truezvsxgraph.Fit("pol1",'Q')
        fitted_vals2 = self.truezvsxgraph.GetFunction("pol1")
        #print "fitted values are"
        #print fitted_vals2.GetParameter(0),fitted_vals2.GetParameter(1),
        #print "**"
        self.trueXZfittedYintercept = fitted_vals2.GetParameter(0)
        self.trueXZfittedYslope = fitted_vals2.GetParameter(1)
        gStyle.SetOptFit(111)
        self.c1.Update()
        self.c1.cd(3)
        gPad.Divide(2,2)
        gPad.cd(1)
        self.leg=TLegend(0.1,0.1,0.7,0.9)
        self.leg.SetFillColor(0)
        self.leg.SetFillStyle(0)
        self.leg.SetBorderSize(0)
        self.leg.AddEntry(self.zvsygraph,'Recorded Hit Position','p')
        self.leg.AddEntry(self.truezvsygraph,'True Hit Position','p')
        #self.leg.AddEntry(self.tseedarrayY,'#splitline{Projected Pos. from T Seed,}{using truth}','p')
        #self.leg.AddEntry(self.tseedarrayX2,'#splitline{Projected Pos. from T Seed,}{using hit pos}','p')
        #self.leg.AddEntry(self.tseedarrayUVLineHit,'#splitline{X position from Y}{using LineHit Method}','p')
        self.leg.AddEntry(self.tseedarrayY,'Projected Pos from T Seed','p')
        self.leg.AddEntry(self.ypreUTMomAdjusted,'y_{PreUT} w/ 1/p Correction','p')
        
        self.leg.SetTextSize(0.05)
        self.leg.Draw()
        self.c1.cd(3)
        gPad.cd(2)
        #print (-1* self.trueXZfittedYintercept / self.trueXZfittedYslope)
        self.val1 = (-1* self.trueXZfittedYintercept / self.trueXZfittedYslope)
        self.val2 = 1./self.trueXZfittedYslope
        self.val3 = (-1* self.trueYZfittedYintercept / self.trueYZfittedYslope)
        self.val4 = 1. / self.trueYZfittedYslope
        self.pt.AddText("#splitline{MC P = %f}{DS P est =%f}"%(self.mcParticleP,self.Tseed.dsP))
        self.pt.AddText('True x_intercept %f'%self.val1)
        self.pt.AddText('True x_slope %f'%self.val2)
        self.pt.AddText('True y_intercept %f'%self.val3)
        self.pt.AddText('True y_slope %f'%self.val4)
        # self.pt.AddText('x1 = %f, x1_guess=%f, chi2=%f'%(self.x1.trueXmid, self.val1+self.val2*self.x1.zmid,(self.x1.trueXmid-(self.val1+self.val2*self.x1.zmid))**2*332.41))
#         angle=(-5.*pi/180.)

#         utrue=self.u.trueXmid*cos(angle)+self.u.trueYmid*sin(angle)
#         uguess=(self.val1+self.val2*self.u.zmid)*cos(angle)+(self.val3+self.val4*self.u.zmid)*sin(angle)

#         self.pt.AddText('u = %f, u_guess=%f, chi2=%f'%(utrue,uguess,
#                                                        (utrue-uguess)**2 * 332.41))

#         vtrue=self.v.trueXmid*cos(angle)-self.v.trueYmid*sin(angle)
#         vguess=(self.val1+self.val2*self.v.zmid)*cos(angle)-(self.val3+self.val4*self.v.zmid)*sin(angle)
#         self.pt.AddText('v = %f, v_guess=%f, chi2=%f'%(vtrue,vguess,
#                                                        (vtrue-vguess)**2 * 332.41))
#         self.pt.AddText('x2 = %f, x2_guess=%f, chi2=%f'%(self.x2.trueXmid, self.val1+self.val2*self.x2.zmid,(self.x2.trueXmid-(self.val1+self.val2*self.x2.zmid))**2*332.41))

#         y0true = self.x1.trueEntryY - self.x1.trueZmid*self.x1.trueTy
#         yTtrue = y0true + self.x1.trueTy*2320.
#         yTproj = self.TseedY0+self.TseedTy*2320.
#         self.pt.AddText('yTtrue = %f, yTproj = %f, '%(yTtrue,yTproj))
#         self.pt.AddText('tyTseed = %f, guess=%f '%(self.x1.trueTy,self.TseedTy))
        
        self.pt.Draw()
        #self.pt2.Draw()
        self.c1.cd(3)
        #gPad.Update()
        gPad.cd(4)
        self.mgTseed = TMultiGraph()
        self.mgTseed.SetTitle("z PreUT vs ty;ty;z[mm]")
        
        self.typreUTmomCor.SetMarkerColor(kMagenta+2)
        self.typreUTmomCor.SetMarkerStyle(30)
        
        self.mgTseed.Add(self.tyPreUTtruth)
        self.mgTseed.Add(self.typreUTmomCor)
        
        self.mgTseed.Draw('ap')
        self.mgTseed.GetYaxis().SetRangeUser(2315,2330)
        self.c1.Update()
        #fit the lines and extract the params for the truth.
        
    def save_curr_graph(self):
        answer=raw_input('save? y or [n]:')

        if 'y' in answer:
            self.c1.SaveAs('x_and_y_projections_ev%d.pdf'%self.evCounter)
        
    
