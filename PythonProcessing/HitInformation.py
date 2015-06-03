#!/usr/local/bin/python
class HitInformation:
    """ A class to read out the information from the printout and put
    it into  a hit, then accessible. hitinfo is then expected to be a list of lines"""
    def __init__(self,hitinfo):
        if not len(hitinfo)==33:#dummy test
            print "did not pass an array of the right length! len(hitinfo)=%i"%len(hitinfo)
            return
        #print hitinfo
        self.plane = int((hitinfo[0].split())[1])
        self.length=float((hitinfo[1].split())[1])
        self.zAtYEq0=float((hitinfo[2].split())[1])
        self.zHit=float((hitinfo[3].split())[1])
        self.zmid=float((hitinfo[4].split())[1])
        self.yBegin=float((hitinfo[5].split())[1])
        self.yEnd=float((hitinfo[6].split())[1])
        self.y=float((hitinfo[7].split())[1])
        self.yMin=float((hitinfo[8].split())[1])
        self.yMid=float((hitinfo[9].split())[1])
        self.yMax=float((hitinfo[10].split())[1])
        self.isX=bool((hitinfo[11].split())[1])
        self.xAtYEq0 = float((hitinfo[12].split())[1])
        self.x=float((hitinfo[13].split())[1])
        self.xT=float((hitinfo[14].split())[1])
        self.xMin=float((hitinfo[15].split())[1])
        self.xMid=float((hitinfo[16].split())[1])
        self.xMax=float((hitinfo[17].split())[1])
        self.dxDy = float((hitinfo[18].split())[1])#for calculation of x from y using linehit method
        self.trueXmid=float((hitinfo[19].split())[1])
        self.trueYmid=float((hitinfo[20].split())[1])
        self.trueZmid=float((hitinfo[21].split())[1])
        self.trueTx=float((hitinfo[22].split())[1])
        self.trueTy=float((hitinfo[23].split())[1])
        self.trueQoP=float((hitinfo[24].split())[1])
        self.trueEntryX=float((hitinfo[25].split())[1])
        self.trueEntryY=float((hitinfo[26].split())[1])
        self.trueEntryZ=float((hitinfo[27].split())[1])
        self.trueHitP=float((hitinfo[28].split())[1]) 
        self.cosTheta = float((hitinfo[29].split())[1])
        self.sinTheta = float((hitinfo[30].split())[1])
        self.weight = float((hitinfo[31].split())[1])
        self.MatchedToSeed = float((hitinfo[32].split())[1]) 
        #print" hit p = %f" %self.trueHitP
        #initialized, now read the hitinfo
        
class TSeedInfo:
    def __init__(self,tseedinfo):
        
        if not len(tseedinfo)==5:
            print 'did not pass an array of the correct length! len(tseedinfo)=%i'%(len(tseedinfo))
            return
        
        self.x=float((tseedinfo[0].split())[1])
        self.y=float((tseedinfo[1].split())[1])
        self.z=float((tseedinfo[2].split())[1])
        self.tx=float((tseedinfo[3].split())[1])
        self.ty=float((tseedinfo[4].split())[1])
        self.y0=self.y - (self.z * self.ty)
        #print 'y0 = %f'%self.y0
        #print self.y0
        # print self.x
        # print self.y
        # print self.z
        # print self.tx
        # print self.ty
    def set_mc_momentum(self,mom):
        self.mcP = mom
    def set_ds_candidate_momentum(self,mom):
        #print 'setting downstream momentum to %f'%mom
        self.dsP = abs(mom)
        self.dsQ = mom/abs(mom)
    def set_endT_info(self,buffer):
        #print buffer
        #print 'setting tseed info to'
        #print buffer
        self.Tend_X = float((buffer[0].split())[2])
        self.Tend_Y = float((buffer[1].split())[2])
        self.Tend_Z = float((buffer[2].split())[2])
        self.Tend_TX = float((buffer[3].split())[2])
        self.Tend_TY = float((buffer[4].split())[2])
        #print 'in set tseed stuff, x=%f, y=%f,z=%f, tx=%f, ty=%f'%(self.Tend_X,self.Tend_Y,self.Tend_Z,self.Tend_TX,self.Tend_TY)
    def setPreselInfo(self,buffer):
        buf =  buffer
        #first kill all the lines in the buffer that don't have the shit we want.
        buf[:] = [line for line in buf if float((line.split())[2]1)]
        print buf
