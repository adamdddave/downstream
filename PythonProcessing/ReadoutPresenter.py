#!/usr/local/bin/python
from ROOT import *
import sys,os, re
from array import array
from HitInformation import *
from graphs2display import*
from Chi2UT import *
#from FitComparer import *
#from EventInspector import *
import numpy
class ReadoutPresenter:
    """reads output for mike and makes a chi2 calculation
      Takes a chunk of text, imported as "chunk" and then gets the 
      values of the hit positions, slopes, momentum, etc and calculates
      a pseudo chi2 out of them.
      """
    pass
    #first,some tests
    def __init__(self,file):
        self.numberHits=0
        self.hitDict={}
        
        
    def chi2(self):
        return 'printing chi2'

#main loop on file
if not len(sys.argv)>1:
    print 'expecting an argument!'
    print 'give a file, and optionally if you want only mc hits y/[n]'
    sys.exit()
f=open(sys.argv[1],'r')
#loop over lines in the file, till we find a "good to go" line, which
#says that we have 4 mc matched hits
mc_matched_hits_only = False
if len(sys.argv) == 3:
    if sys.argv[2] == 'y':
        mc_matched_hits_only = True
found_good=False
good_seed=False
is_x2=False
is_x1=False
is_u=False
is_v=False
x1_done=False
x2_done=False
u_done=False
v_done=False
has_DS_momentum = False
buffer=[]
a=0
ax1=0
au=0
av=0
tinf=0
ev_count=0
chi2val = 0
#this is the final, completely general version. must get all hits first.
hits2consider=[]
x12consider = []
u2consider=[]
v2consider=[]
x22consider=[]
continue_hit_iterations = 'y'
#ty_residual = TH1F("ty_residual","ty_residual;(ty_{T Seed} - ty_{Ideal})",100,10,10)#ty residual
chi2_distribution = TH1D("chi2_distribution","#chi^{2} distribution;#chi^{2};Entries / 1",100,0,100)
delta_chi2_distribution = TH1D("delta_chi2_distribution","#chi^{2}_{1 It} - #chi^{2}_{0 It}; #chi^{2}_{1 It} - #chi^{2}_{0 It}; Entries / 0.5", 100, -50, 50);

# chi2_distribution_mom_range_5_10_3hit = TH1D("chi2_distribution_mom_range_5_10_3hit","#chi^{2} distribution, 3 hits, 5<p<10 GeV;#chi^{2};Entries / 1",100,0,100);
# chi2_distribution_mom_range_5_10_3hit = TH1D("chi2_distribution_mom_range_10_20_3hit","#chi^{2} distribution, 3 hits, 10<p<20 GeV;#chi^{2};Entries / 1",100,0,100);
# chi2_distribution_mom_range_5_10_3hit = TH1D("chi2_distribution_mom_range_20_3hit","#chi^{2} distribution, 3 hits, p>20 GeV;#chi^{2};Entries / 1",100,0,100);

# chi2_distribution_mom_range_5_10_4hit = TH1D("chi2_distribution_mom_range_5_10_4hit","#chi^{2} distribution, 4 hits, 5<p<10 GeV;#chi^{2};Entries / 1",100,0,100);
# chi2_distribution_mom_range_5_10_4hit = TH1D("chi2_distribution_mom_range_10_20_4hit","#chi^{2} distribution, 4 hits, 10<p<20 GeV;#chi^{2};Entries / 1",100,0,100);
# chi2_distribution_mom_range_5_10_4hit = TH1D("chi2_distribution_mom_range_20_4hit","#chi^{2} distribution, 4 hits, p>20 GeV;#chi^{2};Entries / 1",100,0,100);


yat0_good_chi2 = TH1D("yat0_good_chi2","y(z=0),#chi^{2}<15; y(z=0)[mm]; Entries / 10 mm", 200,-1000,1000)
yat0_bad_chi2 = TH1D("yat0_bad_chi2","y(z=0),#chi^{2}>50; y(z=0)[mm]; Entries / 10 mm", 200,-1000,1000)
chi2_vs_delta_chi2 = TH2D("chi2_vs_delta_chi2","#chi^{2} vs #delta #chi^{2}; #delta #chi^2; #chi^2",200,-1000.,1000.,100,0.,100.)
z_magnet = TH1D("z_magnet","z intercept of Tseed and UT; z[mm]; Entries/ 10 mm",150,5000.,6500.)
z_magnet_old = TH1D("z_magnet_old","z intercept of Tseed and UT from MC tuning; z[mm]; Entries/ 10 mm",150,5000.,6500.)
z_magnet_vs_p = TH2D("z_magnet_vs_p","z intercept of Tseed and UT vs p;p[MeV]; z[mm]",1000,0.,100000.,150,5000.,6500.)
z_magnet_vs_1overp = TH2D("z_magnet_vs_1overp","z intercept of Tseed and UT vs 1/p;1/p_{true}[1/MeV]; z[mm]",100,0,1./1000.,150,5000.,6500.)
z_magnet_vs_tx = TH2D("z_magnet_vs_tx","z intercept of Tseed and UT vs tx of tseed;tx_{T Seed}; z[mm]",100, -2., 2. ,150,5000.,6500.)
z_magnet_vs_ty = TH2D("z_magnet_vs_ty","z intercept of Tseed and UT vs ty of tseed;ty_{T Seed}; z[mm]",100, -2., 2. ,150,5000.,6500.)
z_magnet_vs_tx_vs_1overp = TH3D("z_magnet_vs_tx_vs_1overp","z intercept of Tseed and UT vs tx and 1/p;1/p_{True}[1/MeV];tx;z_magnet[mm]",100,0,1./1000,100, -2,2, 150,5000,6500)
deltax_ut_tseed = TH1D("deltax_ut_tseed","#delta x_{Magnet} Between UT and FT;#delta x_{Magnet}[mm]; Entries / 1 mm",100, -50, 50)
deltax_ut_tseed_atzmag = TH1D("deltax_ut_tseed_atzmag","#delta x_{Magnet} at z_{Magnet} Between UT and FT;#delta x_{Magnet}[mm]; Entries / 1 mm",100, -50, 50)
deltax_ut_tseed_param = TH1D("deltax_ut_tseed_param","#delta x_{Magnet} Between UT and FT,parameterized;#delta x_{Magnet}[mm]; Entries / 1 mm",100, -50, 50)
deltax_ut_tseed_param_take2 = TH1D("deltax_ut_tseed_param_take2","#delta x_{Magnet} Between UT and FT,parameterized with 1/p and tx;#delta x_{Magnet}[mm]; Entries / 1 mm",100, -50, 50)
good_seed_counter=0
while not found_good:
    
    line=f.readline()
    if not line:
        print 'eof'
        break
    if 'Start of Tseed Loop' in line:
        #print 'getting seed'
        for i in range(4):
            line=f.readline()#skip the labels
        for i in range(5):
            line=f.readline()
            buffer.append(line.strip('\n'))
        tinf=TSeedInfo(buffer)#got the tseed info
        buffer=[]
    if "Seed Reconstructible              1" in line:
        #print 'found a good seed'
        good_seed_counter +=1
        good_seed=True
        #skip correlations
        for i in range(7):
            line = f.readline()
        for i in range(5):
            line = f.readline()
            #print line
            buffer.append(line.strip('\n'))
        #this is all the information at the back of the T station
        tinf.set_endT_info(buffer)
        buffer = []
        for i in range(20):
            line = f.readline()#skip the zpost UT shiz, as well as the preselection info
            if good_seed_counter>1:
                has_DS_momentum = False
            if 'moment' in line and has_DS_momentum==False:
                #print 'moment test',line
                #print line
                #print good_seed_counter
                #print 'calling set_ds_candidate_momentum'
                tinf.set_ds_candidate_momentum(float(line.split()[1]))
                #print 'called set_ds_candidate_momentum, dsP = %f,dsQ =%f'%(tinf.dsP,tinf.dsQ)
                has_DS_momentum=True
            if 'zMagnet' in line:
                z_magnet_old.Fill(float(line.split()[1]))

    if "Printing Table of Hit Information" in line and good_seed:
        #print 'getting hit parameters'
        line=f.readline()
        #line=f.readline()
        for i in range(36):#was 35, now added matched info
            line=f.readline()
            if 'True Hit info' in line:
                continue
            if line=='\n':
                continue
            if 'DDDB' in line:
                continue
            if "angles are" in line:
                continue
            if 'Plane              3' in line:
                is_x2=True
            if 'Plane              0' in line:
                is_x1=True
            if 'Plane              1' in line:
                is_u=True
            if 'Plane              2' in line:
                is_v=True
            buffer.append(line.strip('\n'))
        #print'buffer before checking plane'
        #print buffer
        #found the hit, now make hit info out of it
        if is_x2:
            #print'in x2, good seed counter = %i'%good_seed_counter
            #print 'got x2 hit'
            #a=HitInformation(buffer)
            x22consider.append(HitInformation(buffer))
            #print buffer
            buffer=[]
            is_x2=False
            x2_done=True
        elif is_x1:
            #print 'got x1'
            #ax1=HitInformation(buffer)
            x12consider.append(HitInformation(buffer))
            is_x1=False
            x1_done=True
            buffer=[]
        elif is_v:
            #print 'found v'
            #av=HitInformation(buffer)
            v2consider.append(HitInformation(buffer))
            buffer=[]
            is_v=False
            v_done=True
        elif is_u:
            #print 'found u'
            #au=HitInformation(buffer)
            u2consider.append(HitInformation(buffer))
            buffer=[]
            is_u=False
            u_done=True
    #if x1_done and x2_done and v_done and u_done and good_seed:
    #two different cases to consider:
    #1. Seed is reconstructible and we found 4 good hits: we'll see "good to go" in the line'
    #2. seed is reconstructible, but didn't find the hits, only get "looping over x hits. There
    # will be no MC momenutm attached here.
    # This is the 5% in efficiency drop from preselection
    if 'Print Good To Go' in line:
            line=f.readline()
            #print line
            if not "m_mc_part->p()" in line:
                x1_done = False
                x2_done = False
                v_done = False
                u_done = False
                continue
            partp = float(line.split('=')[1])
            tinf.set_mc_momentum(partp)
            #print "tseed mc momentum = %f"%partp
            
    if 'Looping on x hits' in line and good_seed:#instead of just looking for good hits, loop over all possibilities
        if mc_matched_hits_only == True:
            #print 'stripping down hits'
            #get rid of all the hits in the containers that aren't matched.
            #also check that the MC momentum is within some tolerance of the seed            
            x12consider[:] = [hit for hit in x12consider if hit.MatchedToSeed==1]
            u2consider[:] = [hit for hit in u2consider if hit.MatchedToSeed==1]
            v2consider[:] = [hit for hit in v2consider if hit.MatchedToSeed==1]
            x22consider[:] = [hit for hit in x22consider if hit.MatchedToSeed==1]
        #do this here because sometimes we won't have "good to go" for a seed.
        if mc_matched_hits_only == False: #don't have mc matched hits, strip those away
            x12consider[:] = [hit for hit in x12consider if hit.MatchedToSeed==0]
            u2consider[:] = [hit for hit in u2consider if hit.MatchedToSeed==0]
            v2consider[:] = [hit for hit in v2consider if hit.MatchedToSeed==0]
            x22consider[:] = [hit for hit in x22consider if hit.MatchedToSeed==0]
        
#         for container in [x12consider,u2consider,v2consider,x22consider]:
#             for hit in container:
#                 print "Plane =%i, momentum = %f, matched  = %i"%(hit.plane,hit.trueHitP, hit.MatchedToSeed)
        #prune down the number of combinations
        for hitlist in [x12consider,u2consider,v2consider,x22consider]:
            if len(hitlist)>30:
                for i in range(len(hitlist)-30):
                    hitlist.pop()
        ev_count+=1
        
        print 'Getting choices for %i x1, %i u, %i v and %i x2 hits'%(len(x12consider),len(u2consider),len(v2consider),len(x22consider))
        ####replace the for loop with itertools. this is just much nicer. to print out the stuff, say list(blah)
        ###first, loop over containers, put into the list to consider if it's not empty
        containers_to_consider=[]
        for container in [x12consider,u2consider,v2consider,x22consider]:
            #print 'looking at container'
            #print 'length of this container = ',len(container)
            if len(container)!=0:
                containers_to_consider.append(container)
        #print 'length of containers_to_consider = %i'%len(containers_to_consider)
        combinations_of_hits=[]
        from itertools import *
        #3 hits is not implemented yet.
        # if(len(containers_to_consider)==3):
        #     #print '3 hit combo!'
        #     combinations_of_hits=list(product(containers_to_consider[0],
        #                                       containers_to_consider[1],
        #                                       containers_to_consider[2]))
        #el
        if(len(containers_to_consider)==4):
            #print '4 hit combo'
            combinations_of_hits=list(product(containers_to_consider[0],
                                              containers_to_consider[1],
                                              containers_to_consider[2],
                                              containers_to_consider[3]))
        #print 'combination of hits = ', combinations_of_hits
        
        #now we only have filled containers
        #for x2hit,x1hit,uhit,vhit in [(x2hit,x1hit,uhit,vhit) for x2hit in x22consider for x1hit in x12consider for uhit in u2consider for vhit in v2consider]:
        for hitcombo in combinations_of_hits:
            # if 'n' in continue_hit_iterations:
            #     break
            # au = uhit
            # av = vhit
            # ax1 = x1hit
            # a = x2hit
            
            ###HERE
            #g2d=graphs2display(ax1,au,av,a)#,partp)
            #g2d.set_counter_pave(ev_count)
            #g2d.setTseedInfo(tinf)
            #g2d.set_tseed_projections()
            #g2d.make_x_y_z_projections()
            #g2d.save_curr_graph()
            #print "tseed ty = %f, tseed y = %f"%(g2d.TseedTy, g2d.TseedY)
            #do we want to calculate the chi2 and put it on the graph?
            #calcChi2 = raw_input('Solve the input, calculate a chi2 and plot on graph? y/[n]')
            #print 'considering hits'
            # for hit in hitcombo:
            #     print 'hit.plane= ',hit.plane
            calcChi2 = 'y'
            if 'y' in calcChi2:
                
                chi2 = Chi2UT(hitcombo)
                #print 'setting tseed params with dsp=%f,dsq=%f'%(tinf.dsP, tinf.dsQ)
                chi2.setTSeedParams(tinf,True)
                chi2.set_WeightMatrixHitErrors(True)
                chi2.determineTseedWeightMatrix(True)
                #chi2.setTSeedParams(tinf,False)
                #chi2.determineTseedWeightMatrix(False)
                chi2.modifyV_X2_positions()
                #chi2.removeHits('x1')
                chi2.setLHSmatrix()
                chi2.CalculateGuess()
                chi2.setRHSvector()
                chi2.calc_best_track()
                chi2val = chi2.Chi2ForParams()
                #print 'chi2val = %f'%chi2val
            #print 'for extracted params'
            #chi2.check_chi2_for_params(chi2.finalParams[0],chi2.finalParams[1],chi2.finalParams[2],chi2.finalParams[3])
            #print 'for initial guess'
            #chi2.check_chi2_for_guess()
            #chi2.check_chi2_for_params(-10.296112,-0.031333,17.361532,-0.081736)
            #print 'for truth'
                #chi2.check_chi2_for_params(g2d.val1,g2d.val2,g2d.val3,g2d.val4)
                #             #chi2.check_chi2_for_true_vals()
                #             curr_yplot = TF1("curr_xplot","[0]+[1]*x",g2d.mg.GetXaxis().GetXmin(),g2d.mg.GetXaxis().GetXmax())
                #             curr_xplot = TF1("curr_yplot","[0]+[1]*x",g2d.mg2.GetXaxis().GetXmin(),g2d.mg2.GetXaxis().GetXmax())
                #             #the parameters have to be inverted again to draw on the same graph
                #             par1=-chi2.finalParams[0]/chi2.finalParams[1]
                #             par2 = 1./chi2.finalParams[1]
                #             par3=-chi2.finalParams[2]/chi2.finalParams[3]
                #             par4 = 1./chi2.finalParams[3]
                #             curr_xplot.SetParameters(par1,par2)
                #             curr_yplot.SetParameters(par3,par4)
                #             for plot in [curr_xplot,curr_yplot]:
                #                 plot.SetLineColor(kGreen+2)
                #                 plot.SetLineStyle(kDashed)
                
                #             g2d.c1.cd(1)
                #             curr_yplot.Draw("lsame")
                #             g2d.c1.Update()
                #             g2d.c1.cd(2)
                #             curr_xplot.Draw("lsame")
                #             g2d.c1.cd(3)
                #             g2d.leg.AddEntry(curr_xplot,"Fit from New #chi^{2}","l")
                #             gPad.cd(1)
                #             g2d.leg.Draw()
                #             fc = FitComparer(chi2)
                #             print 'MC Particle P = %f MeV, Downstream P = %f MeV '%(g2d.mcParticleP,g2d.Tseed.dsP)
                #             fc.print_table()
                #             g2d.c1.Update()
                #g2d.save_curr_graph()
                #chi2.check_chi2_for_params(-10.296112,-0.031333,17.361532,-0.081736)
                #doIterate = raw_input('do another iteration? y / [n]')
                doIterate = 'y'
                if'y' in doIterate:
                    counter = 0
                    while 'y' in doIterate:
                        counter+=1
                    #print 'iteration %d'%counter
                        chi2.setGuess(chi2.finalParams[0],chi2.finalParams[1],chi2.finalParams[2],chi2.finalParams[3])
                        chi2.setRHSvector()
                        chi2.calc_best_track()
                        chi2before = chi2val
                        chi2val=chi2.Chi2ForParams()
                        if(chi2val<15):
                            yat0_good_chi2.Fill(chi2.TSeed.y0)
                        if(chi2val>50):
                            yat0_bad_chi2.Fill(chi2.TSeed.y0)
                        #print 'chi2val = %f'%chi2val
                        delta_chi2_distribution.Fill(chi2val-chi2before)
                        chi2_distribution.Fill(chi2val)
                        #if chi2.TSeed.mcP=>5e3 and chi2.TSeed.mcP<10e3:
                            
                        chi2_vs_delta_chi2.Fill(chi2val-chi2before,chi2val)
                    #print 'delta chi2 = %f'%(chi2val-chi2before)
                        par1=-chi2.finalParams[0]/chi2.finalParams[1]
                        par2 = 1./chi2.finalParams[1]
                        par3=-chi2.finalParams[2]/chi2.finalParams[3]
                        par4 = 1./chi2.finalParams[3]
                        z_magnet.Fill(chi2.check_z_magnet_from_seed())
                        z_magnet_vs_tx.Fill(chi2.TSeed.Tend_TX,chi2.check_z_magnet_from_seed())
                        z_magnet_vs_ty.Fill(chi2.TSeed.Tend_TY,chi2.check_z_magnet_from_seed())
                        # z_magnet_vs_p.Fill(chi2.TSeed.mcP,chi2.check_z_magnet_from_seed())
                        # z_magnet_vs_1overp.Fill(1./chi2.TSeed.mcP,chi2.check_z_magnet_from_seed())
                        z_magnet_vs_p.Fill(abs(chi2.TSeed.dsP),chi2.check_z_magnet_from_seed())
                        z_magnet_vs_1overp.Fill(abs(1./chi2.TSeed.dsP),chi2.check_z_magnet_from_seed())
                        z_magnet_vs_tx_vs_1overp.Fill(1./abs(chi2.TSeed.dsP),chi2.TSeed.Tend_TX,chi2.check_z_magnet_from_seed());
                        #print chi2.check_z_magnet_from_seed()
                        deltax_ut_tseed.Fill(chi2.get_delta_x_interept_for_truth())
                        deltax_ut_tseed_atzmag.Fill(chi2.get_delta_x_interept_for_truth(True))
                        deltax_ut_tseed_param.Fill(chi2.get_delta_x_from_parameterization())
                        deltax_ut_tseed_param_take2.Fill(chi2.get_delta_x_from_parameterization(True))
                        #                     curr_xplot.SetParameters(par1,par2)
                        #                     curr_yplot.SetParameters(par3,par4)
                        #                     g2d.c1.cd(1)
                        #                     curr_yplot.Draw("lsame")
                        #                     g2d.c1.cd(2)
                        #                     curr_xplot.Draw("lsame")
                        #                     g2d.c1.Update()
                        #fc = FitComparer(chi2)
                        #                     fc.print_table()
                        #                     g2d.c1.cd(3)
                        #                     gPad.cd(3)
                        #                     ptnew = TPaveText(0.1,0.1,0.8,0.9)
                        #                     ptnew.SetFillColor(0)
                        #                     ptnew.SetBorderSize(0)
                        #                     ptnew.AddText ("new #chi^{2} = %f"%chi2val)
                        #                     ptnew.Draw()
                        #                     g2d.c1.Update()
                        #                     g2d.save_curr_graph()
                        #
                        #if(chi2val>50):
                        #    ein = EventInspector(g2d,fc,par1,par2,par3,par4)
                        #    doIterate = raw_input('do another iteration? y / [n]')
                        #else:
                        #    doIterate = 'n'
                        
                        doIterate = 'n'
            continue_hit_iterations = 'y'
            chi2_distribution.Fill(chi2val)
            good_seed = False
#        move2next=raw_input('Move to next event? [y] / n')
        move2next = 'y'
        if mc_matched_hits_only==False and ev_count > 300:
            #cut off at 300 combinations.
            found_good = True
#             #print 'breaking loop'
#             #continue
        if not 'n' in move2next:
            #print 'moving to next'
            good_seed = False
            good_seed_counter=0
            x1_done=False
            x2_done=False
            u_done=False
            v_done=False
            has_DS_momentum=False
            x12consider = []
            x22consider = []
            u2consider = []
            v2consider = []
            #print chi2val
            #
            if (ev_count % 100) ==0:
                print 'ev %i'%ev_count
            continue
        #print 'done for now'
        found_good=True
#else:
    #print 'Exited Loop'
f.close()
saveName = './plots/'
typeName = '4hitGhost'
if mc_matched_hits_only==True:
    typeName = '4hitMCMatched'
cfinal = TCanvas("cfinal","cfinal",1200,1200)
cfinal.cd()
gStyle.SetOptStat(111111)
chi2_distribution.Draw()
#print chi2_distribution.GetMean(),chi2_distribution.GetRMS()
cfinal.Update()
cfinal.SaveAs(saveName+"chi2_distribution_"+typeName+".C")
cfinal.SaveAs(saveName+"chi2_distribution_"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"chi2_distribution_"+typeName+"_logy.pdf")
cfinal.SetLogy(False)
#continuenow = raw_input('good?')
cfinal.Clear()
delta_chi2_distribution.Draw()
cfinal.Update()
cfinal.SaveAs(saveName+"delta_chi2_distribution_"+typeName+".C")
cfinal.SaveAs(saveName+"delta_chi2_distribution_"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"delta_chi2_distribution_"+typeName+"_logy.pdf")
cfinal.SetLogy(False)
#continuenow = raw_input('good?')

cfinal.Clear()
yat0_good_chi2.Draw();
cfinal.Update()
cfinal.SaveAs(saveName+"yat0_chi2_lt_15_"+typeName+".C")
cfinal.SaveAs(saveName+"yat0_chi2_lt_15_"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"yat0_chi2_lt_15"+typeName+"_logy.pdf")
cfinal.SetLogy(False)
#continuenow = raw_input('good?')
cfinal.Clear()
yat0_bad_chi2.Draw()
cfinal.Update()
cfinal.SaveAs(saveName+"yat0_chi2_gt_50_"+typeName+".C")
cfinal.SaveAs(saveName+"yat0_chi2_gt_50_"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"yat0_chi2_gt_50"+typeName+"_logy.pdf")
cfinal.SetLogy(False)
#continuenow = raw_input('good?')
cfinal.Clear()
chi2_vs_delta_chi2.Draw("colz")
cfinal.Update()
cfinal.SaveAs(saveName+"chi2_vs_delta_chi2_"+typeName+".C")
cfinal.SaveAs(saveName+"chi2_vs_delta_chi2_"+typeName+".pdf")
cfinal.SetLogz(True)
cfinal.SaveAs(saveName+"chi2_vs_delta_chi2_"+typeName+"_logz.pdf")
cfinal.SetLogz(False)
#continuenow = raw_input('good?')
cfinal.Clear()
z_magnet.Draw()
cfinal.Update()
cfinal.SaveAs(saveName+"zmagnet"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"zmagnet"+typeName+"_logy.pdf")
cfinal.SetLogy(False)
#continuenow = raw_input('good?')
cfinal.Clear()
z_magnet_vs_p.Draw("colz")
cfinal.Update()
cfinal.SaveAs(saveName+"zmagnet_vs_dsp_"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet_vs_dsp_"+typeName+".pdf")
cfinal.SetLogz(True)
cfinal.SaveAs(saveName+"zmagnet_vs_dsp_"+typeName+"_logz.pdv")
cfinal.SetLogz(False)

#continuenow = raw_input('good?')
cfinal.Clear()
z_magnet_vs_1overp.Draw("colz")
cfinal.Update()
cfinal.SaveAs(saveName+"zmagnet_vs_1overdsp_"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet_vs_1overdsp_"+typeName+".pdf")
cfinal.SetLogz(True)
cfinal.SaveAs(saveName+"zmagnet_vs_1overdsp_"+typeName+"_logz.pdv")
cfinal.SetLogz(False)

cfinal.Clear()
z_magnet_old.Draw()
cfinal.Update()
cfinal.SaveAs(saveName+"zmagnet_old"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet_old"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"zmagnet_old"+typeName+"_logy.pdf")
cfinal.SetLogy(False)

cfinal.Clear()
deltax_ut_tseed.Draw()
cfinal.SaveAs(saveName+"delta_xmagnet"+typeName+".C")
cfinal.SaveAs(saveName+"delta_xmagnet"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"delta_xmagnet"+typeName+"_logy.pdf")
cfinal.SetLogy(False)

cfinal.Clear()
deltax_ut_tseed_atzmag.Draw()
cfinal.SaveAs(saveName+"delta_xmagnet_atzmag"+typeName+".C")
cfinal.SaveAs(saveName+"delta_xmagnet_atzmag"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"delta_xmagnet_atzmag"+typeName+"_logy.pdf")
cfinal.SetLogy(False)

cfinal.Clear()
deltax_ut_tseed_param.Draw()
cfinal.SaveAs(saveName+"delta_xmagnet_param"+typeName+".C")
cfinal.SaveAs(saveName+"delta_xmagnet_param"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"delta_xmagnet_param"+typeName+"_logy.pdf")
cfinal.SetLogy(False)

cfinal.Clear()
z_magnet_vs_tx.Draw()
cfinal.SaveAs(saveName+"zmagnet_param_vs_tx"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet_param_vs_tx"+typeName+".pdf")
cfinal.SetLogz(True)
cfinal.SaveAs(saveName+"zmagnet_param_vs_tx"+typeName+"_logz.pdf")

cfinal.Clear()
z_magnet_vs_ty.Draw()
cfinal.SaveAs(saveName+"zmagnet_param_vs_ty"+typeName+".C")
cfinal.SaveAs(saveName+"zmagnet_param_vs_ty"+typeName+".pdf")
cfinal.SetLogz(True)
cfinal.SaveAs(saveName+"zmagnet_param_vs_ty"+typeName+"_logz.pdf")

cfinal.Clear()
deltax_ut_tseed_param_take2.Draw()
cfinal.SaveAs(saveName+"delta_xmagnet_param_with_tx"+typeName+".C")
cfinal.SaveAs(saveName+"delta_xmagnet_param_with_tx"+typeName+".pdf")
cfinal.SetLogy(True)
cfinal.SaveAs(saveName+"delta_xmagnet_param_with_tx"+typeName+"_logy.pdf")
cfinal.SetLogy(False)

cfinal.Clear()
outfile= TFile(saveName+"zmag_vs_tx_vs_1overp.root","recreate")
outfile.cd()
z_magnet_vs_tx_vs_1overp.Write()
outfile.Close()
###Print the chi2 for the last shiz.


