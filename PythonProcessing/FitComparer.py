#!/usr/local/bin/python
from ROOT import*
from array import array
from math import *
from tabulate import tabulate
class FitComparer:
    def theline(self,par0,par1,z):
        return par0 + par1*z
    def __init__(self,chi2):
        self.chi2in = chi2
        #predicted parameters from chi2 fit
        self.predicted_x1 = chi2.finalParams[0]+chi2.finalParams[1]*(chi2.hitx1.zHit)
        self.predicted_u = ((chi2.finalParams[0]+chi2.finalParams[1]*(chi2.hitu.zHit))*chi2.hitu.cosTheta +
                       (chi2.finalParams[2]+chi2.finalParams[3]*(chi2.hitu.zHit))*chi2.hitu.sinTheta)
        self.predicted_v = ((chi2.finalParams[0]+chi2.finalParams[1]*(chi2.hitv.zHit))*chi2.hitv.cosTheta +
                       (chi2.finalParams[2]+chi2.finalParams[3]*(chi2.hitv.zHit))*chi2.hitv.sinTheta)
        self.predicted_x2 = chi2.finalParams[0]+chi2.finalParams[1]*(chi2.hitx2.zHit)
       #  print 'predicted vals = '
#         for pred in [self.predicted_x1,self.predicted_u,self.predicted_v,self.predicted_x2]:
#             print pred
        #true parameters
        self.true_x1 = chi2.hitx1.trueXmid
        self.true_u = chi2.hitu.trueXmid*chi2.hitu.cosTheta + chi2.hitu.trueYmid*chi2.hitu.sinTheta
        self.true_v = chi2.hitv.trueXmid*chi2.hitv.cosTheta + chi2.hitv.trueYmid*chi2.hitv.sinTheta
        self.true_x2 = chi2.hitx2.trueXmid
        #observed positions, as mike says
        self.x1_obs = chi2.hitx1.x
        self.u_obs = chi2.hitu.xMin*chi2.hitu.cosTheta + chi2.hitu.yMin*chi2.hitu.sinTheta
        self.v_obs = chi2.hitv.xMin*chi2.hitv.cosTheta + chi2.hitv.yMax*chi2.hitv.sinTheta
        self.x2_obs = chi2.hitx2.x
        #measured positions, using momentum correction to hits
        self.x1_meas = self.x1_obs
        self.u_meas = self.u_obs
        self.v_meas_correction = -1*(-0.00418*chi2.TSeed.dsQ+1309*chi2.TSeed.dsQ/chi2.TSeed.dsP)
        self.v_meas = self.v_obs + self.v_meas_correction
        self.x2_meas_correction = -1*(-0.0062195*chi2.TSeed.dsQ+1843.*chi2.TSeed.dsQ/chi2.TSeed.dsP)
        self.x2_meas = self.x2_obs + self.x2_meas_correction
        
        #truth after straight line from state before UT

        self.true_x1_from_line = self.theline((chi2.hitx1.trueXmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTx) , chi2.hitx1.trueTx, chi2.hitx1.zHit)
        
        uxline = self.theline((chi2.hitx1.trueXmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTx) , chi2.hitx1.trueTx, chi2.hitu.zHit)
        uyline = self.theline((chi2.hitx1.trueYmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTy) , chi2.hitx1.trueTy,chi2.hitu.zHit)
        
        self.true_u_from_line = uxline*chi2.hitu.cosTheta + uyline*chi2.hitu.sinTheta
        
        vxline = self.theline((chi2.hitx1.trueXmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTx) , chi2.hitx1.trueTx, chi2.hitv.zHit)
        vyline = self.theline((chi2.hitx1.trueYmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTy) , chi2.hitx1.trueTy,chi2.hitv.zHit)
        
        self.true_v_from_line = vxline*chi2.hitv.cosTheta + vyline*chi2.hitv.sinTheta
        
        self.true_x2_from_line = self.theline((chi2.hitx1.trueXmid - chi2.hitx1.trueZmid*chi2.hitx1.trueTx) , chi2.hitx1.trueTx, chi2.hitx2.zHit)
        
        

    def print_table(self):
        #make the table headers
        
        heads = ['Xi','Observed','Measured','Truth From Line','Predicted','Meas - Pred (um)','True - Pred (um)','(Xi_meas - Xi_pred)^2','sigma','Chi2']
        x1_row = ['x1',self.x1_obs,self.x1_meas,self.true_x1_from_line,self.predicted_x1,(self.x1_meas-self.predicted_x1)*1e3,(self.true_x1_from_line-self.predicted_x1)*1e3, 
                  (self.x1_meas-self.predicted_x1)**2 ,(1./self.chi2in.WeightMatrix[0][0])**0.5,(self.x1_meas-self.predicted_x1)**2 * self.chi2in.WeightMatrix[0][0] ]
        
        u_row = ['u',self.u_obs,self.u_meas,self.true_u_from_line,self.predicted_u,(self.u_meas-self.predicted_u)*1e3,(self.true_u_from_line-self.predicted_u)*1e3, 
                 (self.u_meas-self.predicted_u)**2, (1./self.chi2in.WeightMatrix[1][1])**0.5,(self.u_meas-self.predicted_u)**2 * self.chi2in.WeightMatrix[1][1] ]
        v_row = ['v',self.v_obs,self.v_meas,self.true_v_from_line,self.predicted_v,(self.v_meas-self.predicted_v)*1e3,(self.true_v_from_line-self.predicted_v)*1e3,
                 (self.v_meas-self.predicted_v)**2, (1./self.chi2in.WeightMatrix[2][2])**0.5,(self.v_meas-self.predicted_v)**2 * self.chi2in.WeightMatrix[2][2]]
        x2_row = ['x2',self.x2_obs,self.x2_meas,self.true_x2_from_line,self.predicted_x2,(self.x2_meas-self.predicted_x2)*1e3,(self.true_x2_from_line-self.predicted_x2)*1e3, 
                  (self.x2_meas-self.predicted_x2)**2 ,(1./self.chi2in.WeightMatrix[3][3])**0.5,(self.x2_meas-self.predicted_x2)**2 * self.chi2in.WeightMatrix[3][3] ]
        print tabulate([x1_row,u_row,v_row,x2_row],headers=heads,tablefmt="grid")
        yat0 = self.chi2in.TSeed.y-(self.chi2in.tyT*self.chi2in.TSeed.z)
        yTextrap = yat0 + self.chi2in.tyT*self.chi2in.zpUT
        heads_tseed = ['xi', 'Observed', 'Propagated','Measured','Predicted','Truth?', 'Meas-Pred','sigma','Chi2, No Correlation', 'Chi2']
        tyt_row= ['tyT', self.chi2in.TSeed.ty, self.chi2in.TSeed.ty, self.chi2in.tyT,self.chi2in.finalParams[3],self.chi2in.hitx1.trueTy, 
                  self.chi2in.tyT - self.chi2in.finalParams[3],self.chi2in.errorMatrix[1][1]**0.5,(self.chi2in.tyT - self.chi2in.finalParams[3])**2/(self.chi2in.errorMatrix[1][1]), '*']
        
        yt_row = ['yT', self.chi2in.TSeed.y, yTextrap, self.chi2in.yT, self.chi2in.finalParams[2]+ self.chi2in.finalParams[3]*self.chi2in.zpUT, self.chi2in.hitx1.trueEntryY + self.chi2in.hitx1.trueTy*(self.chi2in.hitx1.trueEntryZ - 2320), 
                  self.chi2in.yT-(self.chi2in.finalParams[2]+ self.chi2in.finalParams[3]*self.chi2in.zpUT),self.chi2in.errorMatrix[0][0]**0.5,(self.chi2in.yT-(self.chi2in.finalParams[2]+ self.chi2in.finalParams[3]*self.chi2in.zpUT))**2/(self.chi2in.errorMatrix[0][0]),'*']
        syf = self.chi2in.finalParams[3][0]
        yf = self.chi2in.finalParams[2][0]
        ty = self.chi2in.tyT
        yt = self.chi2in.yT
        zput = self.chi2in.zpUT
        syy = self.chi2in.errorMatrix[0][0]
        styty = self.chi2in.errorMatrix[1][1]
        styy = self.chi2in.errorMatrix[0][1]
        chi2_prop_tseed = (syf - ty)*(-yt*styy + yf*styy+syf*zput*styy-syf*syy + ty*syy)/(styy*styy-styty*syy) +(yt - yf -syf*zput)*(yt*styty - (yf + syf*zput)*styty +(syf-ty)*styy)/(-styy*styy+styty*syy)
        #print chi2_prop_tseed
        yt_row[-1] = chi2_prop_tseed
        #print len(yt_row)
        #final_chi2_Tseed = ['','','','','','','', chi2_prop_tseed]
        #print len(final_chi2_Tseed)
        print tabulate([tyt_row,yt_row],headers=heads_tseed,tablefmt="grid")
        run_chi2 = 0;
        for ding in [x1_row,u_row,v_row,x2_row]:
            run_chi2+=ding[-1][0]
        run_chi2+=yt_row[-1]
        print'Sum of chi2 Contributions = %f'%run_chi2
