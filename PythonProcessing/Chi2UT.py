#!/usr/local/bin/python
from ROOT import *
import numpy#,scipy
from math import *
#from array import *
class Chi2UT:
    def __init__(self,hits,smear = False):
        """
        eliminated the hard coding of cos and sin theta, as well as the weights
        instead we will take this staight from the hits themselves.
        There are now two methods for guessing:
        1. calculate a guess. Good for first time
        2. Use the last parameters as a guess. May help convergence or number of iterations
        finally, using momentum corrections must be added to each step. Sets fine grain control
        5-5-15. Change constructor to list of hits.
        6-7-15. Bug in the modifications of the position. We are oly really assuming x-bending,
                but we adjust the entire v position. Corrected this to only account for v_x
        """
        for hit in hits:
            if hit.plane==0:
                self.hitx1 = hit
            if hit.plane==1:
                self.hitu = hit
            if hit.plane==2:
                self.hitv = hit
            if hit.plane==3:
                self.hitx2 = hit
        #add gaussian smearing
        varx1 = 1./sqrt(self.hitx1.weight)
        varu = 1./sqrt(self.hitu.weight)
        varv = 1./sqrt(self.hitv.weight)
        varx2 = 1./sqrt(self.hitx2.weight)
        myran = numpy.random.normal(0.,varx1)
        
        #here, init all the other stuff.
        self.zx1 = self.hitx1.zmid
        self.zu = self.hitu.zmid
        self.zv = self.hitv.zmid
        self.zx2 = self.hitx2.zmid
        self.zpUT = 2320.
        
        #print"val     true val"
        #print x1.x, x1.trueXmid
        #print u.xMin*u.cosTheta + u.yMin*u.sinTheta, u.trueXmid*u.cosTheta + u.trueYmid*u.sinTheta
        #print v.xMin*v.cosTheta + v.yMax*v.sinTheta, v.trueXmid*v.cosTheta + v.trueYmid*v.sinTheta
        #print x2.x , x2.trueXmid
        #print"***********************"
        
        self.hitx1.distPosition = self.hitx1.x
        self.hitu.distPosition = self.hitu.xMin*self.hitu.cosTheta + self.hitu.yMin*self.hitu.sinTheta
        self.hitv.distPosition = self.hitv.xMin*self.hitv.cosTheta + self.hitv.yMax*self.hitv.sinTheta#min and max to keep the orientation correct
        self.hitx2.distPosition = self.hitx2.x
        if True==smear:
            self.hitx1.distPosition+=numpy.random.normal(0.,varx1)
            self.hitu.distPosition+=numpy.random.normal(0.,varu)
            self.hitv.distPosition+=numpy.random.normal(0.,varv)
            self.hitx2.distPosition+=numpy.random.normal(0.,varx2)
        #true parameters for test
#         self.hitx1.distPosition = x1.trueXmid
#         self.hitu.distPosition = u.trueXmid*u.cosTheta + u.trueYmid*u.sinTheta
#         self.hitv.distPosition = v.trueXmid*v.cosTheta + v.trueYmid*v.sinTheta#min and max to keep the orientation correct
#         self.hitx2.distPosition = x2.trueXmid
        
        self.WeightMatrix = numpy.zeros(shape=(6,6))
        self.aMatrix = numpy.array([[1,       (self.zx1),         0.,       0.],
                                    [self.hitu.cosTheta, (self.zu)*self.hitu.cosTheta, self.hitu.sinTheta, (self.zu)*self.hitu.sinTheta],
                                    [self.hitv.cosTheta, (self.zv)*self.hitv.cosTheta, self.hitv.sinTheta, (self.zv)*self.hitv.sinTheta],
                                    [1.,       (self.zx2),         0.,       0.,],
                                    [0.,       0.,          1.,       (self.zpUT)],
                                    [0.,       0.,          0.,       1.]])
        self.LHSmatrix = numpy.zeros(shape=(4,4))
        self.hits = [self.hitx1,self.hitu,self.hitv,self.hitx2]
        self.hitsRemoved = False
        #print insay

    def removeHits(self,h2r):
        self.hitsRemoved = True
        if 'x1'in h2r:
            #print 'removing x1'
            self.hits = [self.hitu,self.hitv,self.hitx2]
        if 'u'in h2r:
            #print 'removing u'
            self.hits = [self.hitx1,self.hitv,self.hitx2]
        if 'v'in h2r:
            #print 'removing v'
            self.hits = [self.hitx1,self.hitu,self.hitx2]
        if 'x2'in h2r:
            #print 'removing x2'
            self.hits = [self.hitx1,self.hitu,self.hitv]
            
    def setTSeedParams(self,TSeedInfo,useMomentumInfo):
        #self.yT=yT
        self.TSeed = TSeedInfo
        self.tyT = TSeedInfo.ty
        #self.tyT = self.hitx1.trueTy
        #and momentum dependence
        ty_beta = 1.869
        y_beta = -3678;
        if useMomentumInfo==False:
            ty_beta = 0;
            y_beta = 0;
        #self.tyT-=self.hitx1.trueQoP*ty_beta
        self.tyT-=self.TSeed.dsQ/self.TSeed.dsP*ty_beta
        #self.tyT-=self.TSeed.dsQ/self.TSeed.dsP*ty_beta + self.TSeed.dsQ*(-0.003)
        #yat0 = TSeedInfo.y-(self.tyT*TSeedInfo.z)
        yat0 = TSeedInfo.y-(self.TSeed.ty*TSeedInfo.z)
        #yat0 = self.hitx1.trueEntryY - self.hitx1.trueEntryZ*self.hitx1.trueTy
        #self.yT = yat0 + self.hitx1.trueTy*self.zpUT
        #self.yT = TSeedInfo.y - TSeedInfo.ty*(TSeedInfo.z - self.zpUT)
        #self.yT = yat0+self.tyT*self.zpUT
        self.yT = yat0+self.TSeed.ty*self.zpUT
        #self.yT+= y_beta*self.hitx1.trueQoP
        self.yT-= y_beta*self.TSeed.dsQ/self.TSeed.dsP
        

    def modifyV_X2_positions(self):    
        #self.hitv.distPosition-= -0.00418*self.TSeed.dsQ+1309*self.TSeed.dsQ/self.TSeed.dsP
        self.hitv.distPosition = (self.hitv.xMin -(0.004214*self.TSeed.dsQ + 1312.5*self.TSeed.dsQ/self.TSeed.dsP))*self.hitv.cosTheta + self.hitv.yMax*self.hitv.sinTheta
        self.hitx2.distPosition-= -0.0062195*self.TSeed.dsQ+1843.*self.TSeed.dsQ/self.TSeed.dsP
        
    def set_WeightMatrixHitErrors(self,correctVX2errors=False):
        self.WeightMatrix[0][0]=self.hitx1.weight
        #print "(0,0) of weight matrix is set to %f"%self.hitx1.weight
        self.WeightMatrix[1][1]=self.hitu.weight
        #print "(1,1) of weight matrix is set to %f"%self.hitu.weight
        self.WeightMatrix[2][2]=self.hitv.weight
        #print "(2,2) of weight matrix is set to %f"%self.hitv.weight
        self.WeightMatrix[3][3]=self.hitx2.weight
        if correctVX2errors==True:
            #print'correcting v and x2 errors for momentum'
            sigma_corr_v = 0.0070945 + 213.7/self.TSeed.dsP + 1.0178e6/self.TSeed.dsP/self.TSeed.dsP
            sigma_corr_x2 = 0.009321 + 249.95/self.TSeed.dsP + 1.34e6/self.TSeed.dsP/self.TSeed.dsP
            #self.WeightMatrix[2][2]= 1./(1./self.hitv.weight + sigma_corr_v*sigma_corr_v)
            self.WeightMatrix[2][2] = 1./(sigma_corr_v*sigma_corr_v)
            #self.WeightMatrix[3][3] = 1./(1./self.hitx2.weight + sigma_corr_x2*sigma_corr_x2)
            self.WeightMatrix[3][3] = 1./(sigma_corr_x2*sigma_corr_x2)
            #print 'sigma_v**2 = %f, sigma_x2**2 = %f'%(1./self.hitv.weight, 1./self.hitx2.weight)
            #print 'corrections to be added in quadrature are %f for v, %f for x2'%(sigma_corr_v*sigma_corr_v,sigma_corr_x2*sigma_corr_x2)
        #print "(3,3) of weight matrix is set to %f"%self.hitx2.weight
        #print "set Weight matrix to"
        #print self.WeightMatrix

#     def correct_X2_error_with_momentum(self):
#         x2_err = sqrt(1./self.WeightMatrix[3][3])

    def set_WeightMatrixTseedMatrix(self, error_matrix):
        #print 'setting Tseed Weights to ',error_matrix
        self.errorMatrix = error_matrix
        invErrMatrix = numpy.linalg.inv(error_matrix)
        rows,columns = invErrMatrix.shape
        for row in range(rows):
            for col in range(columns):
                self.WeightMatrix[4+row][4+col] = invErrMatrix[row,col]
                #print "setting WeightMatrix[%d,%d] to %f"%(4+row,4+col,invErrMatrix[row,col])
        #print "Weight matrix now"
        #print self.WeightMatrix

    def determineTseedWeightMatrix(self,useMomentumInfo):
        alphaY = 11.87
        betaY = (-6.4035e4)
        gammaY = (4.8315e8)

        alphaTY = 0.0043405
        betaTY = -28.455
        gammaTY = 1.7005e5

        if useMomentumInfo==False:
            print 'not using momentum information'
            betaY = 0.
            gammaY = 0.
            betaTY = 0.
            gammaTY = 0.
        rho = -0.766
        #rho =0 #check to see if this is the problem in the chi2
        #errorY = alphaY+betaY/self.hitx1.trueHitP + gammaY/(self.hitx1.trueHitP*self.hitx1.trueHitP)
        #errorTY = alphaTY+betaTY/self.hitx1.trueHitP + gammaTY/(self.hitx1.trueHitP*self.hitx1.trueHitP)
        #errorTYY = rho*(errorY*errorY*errorTY*errorTY)**0.5
        errorY = alphaY+betaY/abs(self.TSeed.dsP) + gammaY/(self.TSeed.dsP*self.TSeed.dsP)
        errorTY = alphaTY+betaTY/abs(self.TSeed.dsP) + gammaTY/(self.TSeed.dsP*self.TSeed.dsP)
        errorTYY = rho*(errorY*errorY*errorTY*errorTY)**0.5
        errorMatrix = numpy.array([[errorY*errorY,errorTYY],[errorTYY,errorTY*errorTY]])
        self.set_WeightMatrixTseedMatrix(errorMatrix)

    def setLHSmatrix(self):
        sigma_yy = self.errorMatrix[0][0]
        err_yy = sigma_yy
        sigma_tyty = self.WeightMatrix[1][1]
        err_tyty=sigma_tyty
        sigma_tyy = self.errorMatrix[1][0]
        cov_tyy = sigma_tyy
        
        for hit in self.hits:
            self.LHSmatrix[0][0]+= hit.cosTheta*hit.cosTheta*hit.weight
            self.LHSmatrix[0][1]+= hit.zmid*hit.cosTheta*hit.cosTheta*hit.weight
            self.LHSmatrix[0][2]+= hit.cosTheta*hit.sinTheta*hit.weight
            self.LHSmatrix[0][3]+= hit.zmid*hit.cosTheta*hit.sinTheta*hit.weight
            self.LHSmatrix[1][1]+= hit.cosTheta*hit.cosTheta*hit.zmid*hit.zmid*hit.weight
            self.LHSmatrix[1][2]+= hit.cosTheta*hit.sinTheta*hit.zmid*hit.weight
            self.LHSmatrix[1][3]+= hit.cosTheta*hit.sinTheta*hit.zmid*hit.zmid*hit.weight
            self.LHSmatrix[2][2]+= hit.sinTheta*hit.sinTheta*hit.weight
            self.LHSmatrix[2][3]+= hit.sinTheta*hit.sinTheta*hit.zmid*hit.weight
            self.LHSmatrix[3][3]+= hit.sinTheta*hit.sinTheta*hit.zmid*hit.zmid*hit.weight
        #add tseed contributions
        #self.LHSmatrix[2][2]+= sigma_tyty/(sigma_tyty*sigma_yy - sigma_tyy*sigma_tyy)
        #self.LHSmatrix[2][3]+= (sigma_tyy - self.zpUT*sigma_tyty)/(sigma_tyy*sigma_tyy - sigma_tyty*sigma_yy)
        #self.LHSmatrix[3][3]+= (sigma_tyty*self.zpUT*self.zpUT - 2*sigma_tyy*self.zpUT + sigma_yy)/(sigma_tyty*sigma_yy - sigma_tyy*sigma_tyy)
        #try rewriting everything here in gory detail, using different names
        self.LHSmatrix[2][2] += (-1)*err_tyty /(cov_tyy*cov_tyy - err_tyty*err_yy)
        self.LHSmatrix[2][3] += (cov_tyy - self.zpUT*err_tyty)/(cov_tyy*cov_tyy - err_tyty*err_yy)
        self.LHSmatrix[3][3] += ((-2)*cov_tyy*self.zpUT + self.zpUT*self.zpUT*err_tyty + err_yy )/((-1)*cov_tyy*cov_tyy + err_tyty*err_yy)
        #print 'before symmetrizing'
        #print self.LHSmatrix
        #derivative matrix is symmetric
        rows,cols = self.LHSmatrix.shape
        for row in range(rows):
            for col in range(cols):
                if col<=row: continue
                self.LHSmatrix[col][row] = self.LHSmatrix[row][col]
        #print 'after symmetrizing'
        #print self.LHSmatrix
        #print 'determinant of matrix  = %f'%(numpy.linalg.det(self.LHSmatrix))
        #print 'matrix condition number = %f '%(numpy.linalg.cond(self.LHSmatrix))

    def CalculateGuess(self):
        sigma_yy = self.errorMatrix[0][0]
        sigma_tyy = self.errorMatrix[1][0]
        sigma_tyty = self.WeightMatrix[1][1]
        
        self.guessY = self.yT - self.tyT*self.zpUT + 0.001 #add an offset so that the cross terms aren't zero by construction
        self.guessSY = self.tyT + 0.0001
        #guess x and sx by fitting the hit positions.
        
        #self.guessSX = -(self.hitx2.distPosition - self.hitx1.distPosition)/(self.zx2-self.zx1)
        #self.guessX =0.5*( (self.hitx2.distPosition + self.hitx1.distPosition) -self.guessSX*(self.zx2+self.zx1))
        #print "self.guessSX = %f, self.guessX = %f"%(self.guessSX,self.guessX)
        #print "self.zx1 = %f, self.zx2 = %f"%(self.zx1,self.zx2)
        zmin = 2700
        zmax = 2000
        xmin,xmax = -9999,-9999
        for hit in self.hits:
            if hit.plane<2:
                if hit.zmid < zmin:
                    zmin = hit.zmid
                    #xmin = hit.distPosition*hit.cosTheta
                    xmin = hit.x
            else:
                if hit.zmid > zmax:
                    zmax = hit.zmid
                    #xmax = hit.distPosition*hit.cosTheta
                    xmax = hit.x
        #print "zmin = %f zmax = %f"%(zmin,zmax)
        
        guessSXnew = -(xmax - xmin)/(zmax-zmin)
        guessXnew = 0.5*( (xmax + xmin) -guessSXnew*(zmax+zmin))
        #print "guessSXnew = %f, guessXnew = %f"%(guessSXnew,guessXnew)
        self.guessSX = guessSXnew
        self.guessX = guessXnew
        
    def setGuess(self,xguess,sxguess,yguess,syguess):
        self.guessX = xguess
        self.guessSX = sxguess
        self.guessY = yguess
        self.guessSY = syguess

    def setRHSvector(self):
        sigma_yy = self.errorMatrix[0][0]
        sigma_tyy = self.errorMatrix[1][0]
        sigma_tyty = self.WeightMatrix[1][1]                
        #print "guessY = %f, guess SY = %f, guess X = %f, guess SX = %f"%(self.guessY,self.guessSY,self.guessX, self.guessSX)
        self.RHSvector = numpy.zeros(shape=(4,1))
        #guessPosition = 0;
        
        for hit in self.hits:
            guessPosition = (self.guessX + self.guessSX*hit.zmid)*hit.cosTheta + (self.guessY + self.guessSY*hit.zmid)*hit.sinTheta 
            
#             if hit.plane==2:#v plane correction
#                 #print'v plane test'
#                 #print 'original position  = %f, guess_position = %f'%(hit.distPosition,guessPosition)
#                      #alpha = 1341.5*self.hitx1.trueQoP
#                 #alpha = 1341.5*self.TSeed.dsQ/self.TSeed.dsP
#                 alpha = -0.00418*self.TSeed.dsQ+1309*self.TSeed.dsQ/self.TSeed.dsP
#                 #print 'alpha = %f '%alpha
#                 guessPosition-=alpha
#             # if hit.plane==3:
# #                 #alpha = 1890.*self.hitx1.trueQoP
# #                 print'x2 plane test'
# #                 print 'original position  = %f, guess_position = %f'%(hit.distPosition,guessPosition)
# #                 alpha = -0.0062195*self.TSeed.dsQ+1843.*self.TSeed.dsQ/self.TSeed.dsP
# #                 print 'alpha = %f '%alpha
# #                 guessPosition-=alpha
#             #sin theta kills off the non Y terms, and cosTheta = 1 for x terms

            self.RHSvector[0]+=(hit.distPosition-guessPosition)*hit.weight
            self.RHSvector[1]+=hit.zmid*(hit.distPosition - guessPosition)*hit.weight
            self.RHSvector[2]+=hit.sinTheta*( hit.distPosition- guessPosition)*hit.weight
            self.RHSvector[3]+=hit.zmid*hit.sinTheta*( hit.distPosition -guessPosition )*hit.weight

        self.RHSvector[2]+=(self.yT - self.guessY - self.guessSY*self.zpUT)*sigma_tyty/(sigma_tyty*sigma_yy - sigma_tyy*sigma_tyy)
        self.RHSvector[2]+=(self.tyT-self.guessSY)*sigma_tyy/(sigma_tyty*sigma_yy-sigma_tyy*sigma_tyy)

        self.RHSvector[3]+=(self.yT - self.guessY - self.guessSY*self.zpUT)*(self.zpUT*sigma_tyty - sigma_tyy)/(sigma_tyty*sigma_yy - sigma_tyy*sigma_tyy)
        self.RHSvector[3]+=(self.guessSY - self.tyT)*(self.zpUT*sigma_tyy - sigma_yy)/(sigma_tyty*sigma_yy - sigma_tyy*sigma_tyy)
        #print "set RHS to"
        #print self.RHSvector
    def calc_best_track(self):    
        
        L = numpy.linalg.cholesky(self.LHSmatrix)
        y =numpy.linalg.solve(L,self.RHSvector)
        #print y
        #now solve the other problem
        Lh = L.conj().T
        
        self.finalDelta = numpy.linalg.solve(Lh,y)
        #print 'solving gives'
        #print self.finalDelta
        #print 'guessX = %f, correction X = %f, final param = %f'%(self.guessX, self.finalDelta[0],self.guessX+self.finalDelta[0])
        #print 'guessSX = %f, correction SX = %f, final param = %f'%(self.guessSX, self.finalDelta[1],self.guessSX+self.finalDelta[1])
        #print 'guessY = %f, correction Y = %f, final param = %f'%(self.guessY, self.finalDelta[2],self.guessY+self.finalDelta[2])
        #print 'guessSY = %f, correction SY = %f, final param = %f'%(self.guessSY, self.finalDelta[3],self.guessSY+self.finalDelta[3])
        self.finalParams = numpy.zeros(shape=(4,1))
        self.finalParams[0]=self.guessX + self.finalDelta[0]
        self.finalParams[1]=self.guessSX + self.finalDelta[1]
        self.finalParams[2]=self.guessY + self.finalDelta[2]
        self.finalParams[3]=self.guessSY + self.finalDelta[3]
#         print "from fit, final params are"
#         print self.finalParams
    def Chi2ForParams(self):
        
        xiFromBefore = numpy.zeros(shape=(6,1))
        # xiFromBefore[0][0]=self.hitx1.distPosition
#         xiFromBefore[1][0]=self.hitu.distPosition
#         xiFromBefore[2][0]=self.hitv.distPosition
#         xiFromBefore[3][0]=self.hitx2.distPosition
#         xiFromBefore[4][0]=self.yT
#         xiFromBefore[5][0]=self.tyT
        #print "chi2forparams, original xi"
        #for entry in xiFromBefore:
        #    print entry
        #newXi = numpy.zeros(shape=(6,1))
        for hit in self.hits:
            xiFromBefore[hit.plane][0]=hit.distPosition
        xiFromBefore[4][0]=self.yT
        xiFromBefore[5][0]=self.tyT
        #print'xi = '
        #print xiFromBefore
        #print 'amatrix = '
        #print self.aMatrix
        #print 'final params are'
        #print self.finalParams
        #print "new xi"
        #for entry in newXi:
        #    print entry
        guessParams = numpy.dot(self.aMatrix,self.finalParams)
        if True==self.hitsRemoved:
            #find the missing hit
            run_count = 0
            index2delete = -1;
            for hit in self.hits:
                print 'looking at plane %i'%hit.plane
                run_count+=hit.plane
            if run_count==3:
                print 'missing x2'
                index2delete = 3
            if run_count==4:
                print 'missing v'
                index2delete = 2
            if run_count==5:
                print 'missing u'
                index2delete = 1
            if run_count==6:
                print 'missing x1'
                index2delete = 0
            print 'index2delete = %i'%index2delete
            guessParams[index2delete][0]=0
        #vec2mult  = numpy.dot(self.aMatrix,self.finalParams)-xiFromBefore
        vec2mult  = guessParams-xiFromBefore
        #print 'vect2mult = '
        #print vec2mult
        #print "getting chi2"
        #print 'todot = '
        todot = numpy.dot(self.WeightMatrix,vec2mult)
        #print todot
        chi2out = numpy.dot(vec2mult.transpose(),todot)
        #print chi2out
        return chi2out[0][0]
            
        
    def check_chi2_for_params(self,xin,sxin,yin,syin):
        #print 'checking chi2 for parameters x = %f, sx = %f, y=%f, sy=%f'%(xin,sxin,yin,syin)
        xiFromBefore = numpy.zeros(shape=(6,1))
        xiFromBefore[0][0]=self.hitx1.distPosition
        xiFromBefore[1][0]=self.hitu.distPosition
        xiFromBefore[2][0]=self.hitv.distPosition
        xiFromBefore[3][0]=self.hitx2.distPosition
        xiFromBefore[4][0]=self.yT
        xiFromBefore[5][0]=self.tyT

        #        print 'measurements are',xiFromBefore
        test_params=numpy.zeros(shape=(4,1))
        test_params[0]=xin
        test_params[1]=sxin
        test_params[2]=yin
        test_params[3]=syin
        #       print 'parameters as a vector are',test_params
        vec2mult = numpy.dot(self.aMatrix,test_params)-xiFromBefore
        #      print 'a.alpha -xi = ',vec2mult
        half=numpy.dot(self.WeightMatrix,vec2mult)
        rest = numpy.dot(vec2mult.transpose(),half)
        print "chi2 = %f"%rest
    
    
    def display_weight_matrix_as_image(self):
        from matplotlib import pyplot as plt
        plt.imshow(self.WeightMatrix, interpolation='nearest')
        plt.show()
    def display_xi(self):
        from matplotlib import pyplot as plt
        xiFromBefore = numpy.zeros(shape=(6,1))
        xiFromBefore[0][0]=self.hitx1.distPosition
        xiFromBefore[1][0]=self.hitu.distPosition
        xiFromBefore[2][0]=self.hitv.distPosition
        xiFromBefore[3][0]=self.hitx2.distPosition
        xiFromBefore[4][0]=self.yT
        xiFromBefore[5][0]=self.tyT
        for entry in xiFromBefore:
            print entry
        plt.imshow(xiFromBefore, interpolation='nearest')
        plt.show()
    
    def check_u_projections(self,uhit):
        print 'u.xMin = %f, u.xMid = %f, u.xMax = %f, u.x = %f, u.trueXmid = %f '%(uhit.xMin,uhit.xMid,uhit.xMax,uhit.x,uhit.trueXmid)
        print 'u.yMin = %f, u.yMid = %f, u.yMax = %f, u.y = %f, u.trueYmid = %f'%(uhit.yMin,uhit.yMid,uhit.yMax,uhit.y,uhit.trueYmid)
        cT = cos(5.*pi/180.)
        sT = -1*sin(5.*pi/180.)
        print 'u from min = %f, u from mid = %f, u from max = %f, u from nom = %f, u from true%f'%(uhit.xMin*cT + uhit.yMin*sT,
                                                                                                   uhit.xMid*cT + uhit.yMid*sT,
                                                                                                   uhit.xMax*cT + uhit.yMax*sT,
                                                                                                   uhit.x*cT + (-1)*uhit.y*sT,
                                                                                                   uhit.trueXmid*cT+uhit.trueYmid*sT
                                                                                                   )

    def check_v_projections(self,vhit):
        print 'v.xMin = %f, v.xMid = %f, v.xMax = %f, v.x = %f'%(vhit.xMin,vhit.xMid,vhit.xMax,vhit.x)
        print 'v.yMin = %f, v.yMid = %f, v.yMax = %f, v.y = %f'%(vhit.yMin,vhit.yMid,vhit.yMax,vhit.y)
        print 'v.xMin = %f, v.xMid = %f, v.xMax = %f, v.x = %f, v.trueXmid = %f '%(vhit.xMin,vhit.xMid,vhit.xMax,vhit.x,vhit.trueXmid)
        print 'v.yMin = %f, v.yMid = %f, v.yMax = %f, v.y = %f, v.trueYmid = %f'%(vhit.yMin,vhit.yMid,vhit.yMax,vhit.y,vhit.trueYmid)
        cT = cos(5.*pi/180.)
        sT = sin(5.*pi/180.)
        print 'v from min = %f, v from mid = %f, v from max = %f, v from nom = %f, v from true %f'%(vhit.xMin*cT + vhit.yMax*sT,
                                                                                                    vhit.xMid*cT + vhit.yMid*sT,
                                                                                                    vhit.xMax*cT + vhit.yMin*sT,
                                                                                                    vhit.x*cT + vhit.y*sT,
                                                                                                    vhit.trueXmid*cT+vhit.trueYmid*sT
                                                                                    )
    def check_chi2_for_guess(self):
        self.check_chi2_for_params(self.guessX,self.guessSX,self.guessY,self.guessSY)
            
        
    def check_z_magnet_from_seed(self):
        """
        Check whether the line from the rear of the tseed intersects with the line from the UT at the magnet
        bend point. First, we'll need to tuple the intersection. Then, if it's a single point for all momenta,
        we'll just make a check of what the x positions are
        """
        x0fromBackofT = self.TSeed.Tend_X - self.TSeed.Tend_Z*self.TSeed.Tend_TX
        #print 'x0fromBackofT = %f'%x0fromBackofT
        #line_from_t  = x0fromBackofT + self.TSeed.TendTX*z
        #line_from_ut = self.finalDelta[0]+self.finalDelta[1]*z
        z=(self.finalParams[0]-x0fromBackofT)/(self.TSeed.Tend_TX-self.finalParams[1])
        #print 'should return z = ',z
        return z[0]
    def check_zmag_param(self):
        p0 = 5393.
        p1 = -1.363e5
        p2 = 1.102e9
        z = p0 + p1/abs(self.TSeed.dsP) + p2/self.TSeed.dsP/self.TSeed.dsP
        return z
    def check_zmag_param_with_tx(self):
        p0 = 5393.
        p1 = -1.363e5
        p2 = 1.102e9
        p3 = 11.36
        p4 = 335.3
        z = p0 + p1/abs(self.TSeed.dsP) + p2/self.TSeed.dsP/self.TSeed.dsP + p3*self.TSeed.Tend_TX + p4*(self.TSeed.Tend_TX*self.TSeed.Tend_TX)
        return z
    def get_delta_x_interept_for_truth(self,usezmag = False):
        #return the difference between the UT and FT x magnet positions.
        x0fromBackofT = self.TSeed.Tend_X - self.TSeed.Tend_Z*self.TSeed.Tend_TX
        #
        z=5430.#mm,, peak position only
        if usezmag==True:
            z=(self.finalParams[0]-x0fromBackofT)/(self.TSeed.Tend_TX-self.finalParams[1])   
        line_from_t  = x0fromBackofT + self.TSeed.Tend_TX*z
        line_from_ut = self.finalParams[0]+self.finalParams[1]*z
        return line_from_ut - line_from_t
    def get_delta_x_from_parameterization(self,withTX=False):
        x0fromBackofT = self.TSeed.Tend_X - self.TSeed.Tend_Z*self.TSeed.Tend_TX
        p0 = 5393.
        p1 = -1.363e5
        p2 = 1.102e9
        z = p0 + p1/abs(self.TSeed.dsP) + p2/self.TSeed.dsP/self.TSeed.dsP
        if withTX==True:
            p3 = 11.36
            p4 = 335.3
            z+=p3*self.TSeed.Tend_TX + p4*(self.TSeed.Tend_TX*self.TSeed.Tend_TX)
        line_from_t  = x0fromBackofT + self.TSeed.Tend_TX*z
        line_from_ut = self.finalParams[0]+self.finalParams[1]*z

        return line_from_ut - line_from_t

    def get_delta_x_from_parameterization_of_tx(self):
        x0fromBackofT = self.TSeed.Tend_X - self.TSeed.Tend_Z*self.TSeed.Tend_TX
        p0 = 5393.
        p1 = 11.36
        p2 = 335.3
        z = p0 + p1*self.TSeed.tx + p2*self.TSeed.tx
        
        line_from_t  = x0fromBackofT + self.TSeed.Tend_TX*z
        line_from_ut = self.finalParams[0]+self.finalParams[1]*z

        return line_from_ut - line_from_t
    def get_delta_x_from_previous_magnetpos(self):
        line_from_ut = self.finalParams[0]+ self.finalParams[1]*self.TSeed.preselZmag
        return line_from_ut - self.TSeed.preselXmag
