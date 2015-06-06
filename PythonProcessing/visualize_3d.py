#!/usr/local/bin/python
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

class visualize_3d:
    def __init__(self,hits,chi2):
        for hit in hits:
            if hit.plane==0:
                self.hitx1 = hit
            if hit.plane==1:
                self.hitu = hit
            if hit.plane==2:
                self.hitv = hit
            if hit.plane==3:
                self.hitx2 = hit
        self.chi2 = chi2
            
    def plot(self):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.plot([self.hitx1.xMin,self.hitx1.xMid,self.hitx1.xMax],
                [self.hitx1.yMin,self.hitx1.yMid,self.hitx1.yMax],
                [self.hitx1.zmid,self.hitx1.zmid,self.hitx1.zmid],label='x1')

        ax.plot([self.hitx2.xMin,self.hitx2.xMid,self.hitx2.xMax],
                [self.hitx2.yMin,self.hitx2.yMid,self.hitx2.yMax],
                [self.hitx2.zmid,self.hitx2.zmid,self.hitx2.zmid],label='x2')
        #draw boxes for the U and V
        ax.plot([self.hitu.xMin,self.hitu.xMax,self.hitu.xMax,self.hitu.xMin,self.hitu.xMin],
                [self.hitu.yMin,self.hitu.yMin,self.hitu.yMax,self.hitu.yMax,self.hitu.yMin],
                [self.hitu.zmid,self.hitu.zmid,self.hitu.zmid,self.hitu.zmid,self.hitu.zmid],label='u')
        ax.plot([self.hitv.xMin,self.hitv.xMax,self.hitv.xMax,self.hitv.xMin,self.hitv.xMin],
                [self.hitv.yMin,self.hitv.yMin,self.hitv.yMax,self.hitv.yMax,self.hitv.yMin],
                [self.hitv.zmid,self.hitv.zmid,self.hitv.zmid,self.hitv.zmid,self.hitv.zmid],label='v')
        ax.plot([self.hitu.xMin,self.hitu.xMax],
                [self.hitu.yMin,self.hitu.yMax],
                [self.hitu.zmid,self.hitu.zmid],label='uline')
        #add the true u and v positions to be sure we are correct.
        ax.plot([self.hitv.xMin,self.hitv.xMax],
                [self.hitv.yMax,self.hitv.yMin],
                [self.hitv.zmid,self.hitv.zmid],label='vline')
        ax.scatter([self.hitu.trueXmid,self.hitv.trueXmid],
                   [self.hitu.trueYmid,self.hitv.trueYmid],
                   [self.hitu.trueZmid,self.hitv.trueZmid],c='r', marker='o')
        #now draw the fit.
        zplots = np.linspace(2320,self.hitx2.zmid, 100)#z range to draw the plots
        xplots = self.chi2.finalParams[0]+self.chi2.finalParams[1]*zplots
        yplots = self.chi2.finalParams[2]+self.chi2.finalParams[3]*zplots
        ax.plot(xplots,yplots,zplots,label='fit')
        plt.show()
