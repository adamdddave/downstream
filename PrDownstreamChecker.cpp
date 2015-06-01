// Include files 

#ifndef PrDownstreamChecker_CPP
#define PrDownstreamChecker_CPP
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/PhysicalConstants.h"
#include "Event/FitNode.h"
#include "Event/STCluster.h"
#include "Kernel/HitPattern.h"
#include "Linker/AllLinks.h"
#include "Linker/LinkedFrom.h"
#include "Event/StateParameters.h"

#include "PrDownstreamChecker.h"


DECLARE_ALGORITHM_FACTORY( PrDownstreamChecker )

bool dAncestor(const LHCb::MCParticle* mcPart){
  bool fromD = false;
  const LHCb::MCParticle* mother = mcPart->mother();
  while ( mother !=0 && fromD == false) {
    fromD = mother->particleID().hasCharm();
    mother = mother->mother();
  }
  return fromD;
}



//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrDownstreamChecker::PrDownstreamChecker(const std::string& name, ISvcLocator* pSvcLocator ) :
  GaudiTupleAlg( name , pSvcLocator ),
  m_associator(0,""), m_associator_seed(0,""), m_associator_updatedseed(0,""),evt(-1)
{
  //locations of the downstream and seed containers
  declareProperty("DownstreamContainer", m_DownstreamContainer = LHCb::TrackLocation::Downstream);
  declareProperty("TSeedContainer", m_TSeedContainer = LHCb::TrackLocation::Seed);
  declareProperty("UpdatedTSeed_Container",m_Test_Container = LHCb::TrackLocation::Downstream +"UpdatedSeed");
}

//=============================================================================
// Destructor
//=============================================================================
PrDownstreamChecker::~PrDownstreamChecker() {;}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrDownstreamChecker::initialize()
{
  // Mandatory initialization of GaudiAlgorithm
  StatusCode sc = GaudiTupleAlg::initialize();
  if ( sc.isFailure() ) { return sc; }
  
  m_selector = tool<IMCReconstructible>("MCReconstructible","Selector",this);

  return StatusCode::SUCCESS;
}


//=============================================================================
// Execute
//=============================================================================
StatusCode PrDownstreamChecker::execute()
{
  StatusCode sc = StatusCode::SUCCESS;
  
  evt++;
  
  // check input data
  if (!exist<LHCb::Tracks>(m_DownstreamContainer))
    return Warning( m_DownstreamContainer+" not found", StatusCode::SUCCESS, 0);
  if (!exist<LHCb::Tracks>(m_TSeedContainer))
    return Warning( m_TSeedContainer+" not found", StatusCode::SUCCESS, 0);

  //ad
  if (!exist<LHCb::Tracks>(m_Test_Container))
    return Warning( m_Test_Container+" not found", StatusCode::SUCCESS, 0);
  //ad
    
  LHCb::Tracks* downstream = get<LHCb::Tracks>(m_DownstreamContainer);
  nDownstreamTracks = downstream->size();
  LHCb::Tracks* seed = get<LHCb::Tracks>(m_TSeedContainer);
  nSeedTracks = seed->size();

  if ( msgLevel(MSG::DEBUG) ){  
    //ad. Test the downstream test container and make sure we can get the extra info
    LHCb::Tracks* downstream_test = get<LHCb::Tracks>(m_DownstreamContainer);
    //  nDownstreamTestTracks = downstream_test->size();
    for(LHCb::Tracks::const_iterator itT = downstream_test->begin();
        downstream_test->end()!=itT; ++itT){
      LHCb::Track* tester = *itT;
      debug()<<"test info 1800 = "<<tester->info(1800,-999)<<endmsg;
      //info()<<"test info 21 = "<<tester->info(LHCb::Track::FitTChi2,-999)<<endmsg;
    }
    LHCb::Tracks* updatedSeed_test = get<LHCb::Tracks>(m_Test_Container);
    //  nDownstreamTestTracks = downstream_test->size();
    for(LHCb::Tracks::const_iterator itT = updatedSeed_test->begin();
        updatedSeed_test->end()!=itT; ++itT){
      LHCb::Track* tester = *itT;
      debug()<<"test info for seed 1800 = "<<tester->info(1800,-999)<<endmsg;
      //info()<<"test info 21 = "<<tester->info(LHCb::Track::FitTChi2,-999)<<endmsg;
    }
    
  }
  
    //end test ad
  //nSeedTracks = 
  //int nDownstreamTracks;
  
  // get the association tables
  m_associator  = AsctTool(evtSvc(), m_DownstreamContainer);
  m_directTable = m_associator.direct();
  if (!m_directTable)
    return Error("Failed to find direct table for Downstream tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  m_inverseTable = m_associator.inverse();
  if (!m_inverseTable)
    return Error("Failed to find inverse table for Downstream tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  
  m_associator_seed  = AsctTool(evtSvc(), m_TSeedContainer);
  m_directTable_seed = m_associator_seed.direct();
  if (!m_directTable_seed)
    return Error("Failed to find direct table for Seed tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  m_inverseTable_seed = m_associator_seed.inverse();
  if (!m_inverseTable_seed)
    return Error("Failed to find inverse table for Seed tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);

  //AD 1-10-14
  
  m_associator_updatedseed  = AsctTool(evtSvc(), m_Test_Container);
  m_directTable_updatedseed = m_associator_updatedseed.direct();
  if (!m_directTable_updatedseed)
    return Error("Failed to find direct table for Updated Seed tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  m_inverseTable_updatedseed = m_associator_updatedseed.inverse();
  if (!m_inverseTable_updatedseed)
    return Error("Failed to find inverse table for Updated Seed tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
    
  // loop over MC particles
  if ( loopOverMCParticles().isFailure() ) return Error("Failed to loop over MC particles", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  
  // loop over tracks
  if ( loopOverTracks().isFailure() ) return Error("Failed to loop over tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  
  // loop over tracks
  if ( loopOverSeeds().isFailure() ) return Error("Failed to loop over tracks", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  
  //loop over updated tseeds. AD 1-10-14
  if ( loopOverUpdatedSeeds().isFailure() ) return Error("Failed to loop over updated tseeds", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
  return StatusCode::SUCCESS;
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrDownstreamChecker::finalize()
{
  return GaudiTupleAlg::finalize();
}


//=============================================================================
// Loop over MC particles and look for reconstructed tracks
//=============================================================================
StatusCode PrDownstreamChecker::loopOverMCParticles()
{
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> loopOverMCParticles" << endmsg;
  
  n_in_evt = -1;
  
  MCTrackInfo trackInfo( evtSvc(), msgSvc() );
  
  // book the ntuple
  Tuples::Tuple ntuple = GaudiTupleAlg::nTuple("MCParticles","MCParticles", CLID_ColumnWiseTuple );
  LHCb::Tracks* downstream_tracks = get<LHCb::Tracks>(m_DownstreamContainer);     
  // loop over MC particles
  const LHCb::MCParticles* particles = get<LHCb::MCParticles>(LHCb::MCParticleLocation::Default);
  int nTrack(0);
  for (LHCb::MCParticles::const_iterator ip = particles->begin(); ip != particles->end(); ++ip)
    {
      
      debug()   << "processing MC particle #" << nTrack << "/" << particles->size() <<  endmsg ;
      
      nTrack++;
      if ( msgLevel(MSG::DEBUG) ) debug() << "==> PID = " << fabs((*ip)->particleID().pid()) << endmsg;
      
      //reject electrons at this stage
      if ( fabs((*ip)->particleID().pid())==11 ) continue;
      
      // information for the tuple
      int nHits = -999;
      //rest moved to "fill MC particle. 1-15-14
      /*
      bool reconstructible_asLong(false);
      bool B_child(false);
      bool D_child(false);
      bool isLong(false);
      bool isDown(false);
      bool isInVelo(false);*/
      bool over5(false);
      int parentABSID = -999;//add to MC partcilefill. AD 1-12-14
      /*
      int has_all_info(-1);

      //need to make sure we have full track info.
      has_all_info = 
        trackInfo.fullInfo( *ip );

      reconstructible_asLong = 
        m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,*ip);

      
      isLong = 
        trackInfo.hasVeloAndT( *ip );
      
      isDown = 
        trackInfo.hasT( *ip ) &&  trackInfo.hasTT( *ip );
      */
      over5 = 
        5000. < fabs( (*ip)->p() );
      /*
      isInVelo = 
        trackInfo.hasVelo( *ip );
      */
      if ( 0 != (*ip)->originVertex() ) {
        const LHCb::MCParticle* mother =  (*ip)->originVertex()->mother();
        if ( 0 != mother ) {
          if ( 0 != mother->originVertex() ) {
            //double rOrigin = mother->originVertex()->position().rho();
            parentABSID = abs( mother->particleID().pid() );
          }
        }
      }

      /*
      if(bAncestor(*ip)){
        B_child = true;
      }
      if(dAncestor(*ip)){
        D_child = true;
      }*/
      //here
      nHits = getUT_hits(*ip);
      
      
      if(!over5) continue;
      //if(!over5) continue;
      if(parentABSID!=310) continue;
      
      /// good candidate
      n_in_evt++;
      
      /// fill the tuple with MC information
      /*      ntuple->column("parentABSID", parentABSID);
      ntuple->column("fullInfo", has_all_info);
      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      ntuple->column("isLong", isLong);
      ntuple->column("isDown", isDown);
      ntuple->column("over5", over5);
      ntuple->column("isInVelo", isInVelo);
      ntuple->column("B_child", B_child);
      ntuple->column("D_child", D_child);*/
      ntuple->column("nUT_hits", nHits);
      fillMCParticle(*ip,ntuple);
      
      /// now start to look for the reconstructed particle
      AsctTool::InverseType::Range range = m_inverseTable->relations(*ip);
      bool reconstructed(false);
      bool reconstructed_seed(false);
      unsigned int clone(0u);
      double purity(-1.);
      
      reconstructed = !(range.empty());
      const LHCb::Track *track(0);
      if (reconstructed)
	{
	  track = range.begin()->to();
	  clone = range.size()-1u;
	  purity = range.begin()->weight();
	}
      
      range = m_inverseTable_seed->relations(*ip);
      reconstructed_seed = !(range.empty());
      
      fillTrack(track,ntuple);
      ntuple->column("reconstructed", reconstructed);
      ntuple->column("nClones", clone);
      ntuple->column("purity", purity);
      ntuple->column("reconstructed_seed", reconstructed_seed);
      
      
      
      // Loop over the Tracks      
      LinkedFrom<LHCb::STCluster,LHCb::MCParticle> ttLink(evtSvc(),msgSvc(), LHCb::STClusterLocation::UTClusters);
      /*if (!ttLink.notFound()){
	debug() << "now loop over UT clusters for this MC particle" << endmsg;
	int nhits = 0;
	const LHCb::STCluster* TTCluster = ttLink.first((*ip));
	for ( ; 0 != TTCluster; TTCluster = ttLink.next()) {
	  //if ( !TTCluster->isUT() ) continue;
	  //LHCbID : { UT STChannelID : 20128867 : type=2 strip=99 sector=73 detRegion=1 layer=2 station=1 } }
	  nhits++;
	  debug() << "hit " << nhits 
	    //<< ": type=" << TTCluster->type() 
		  << ": ID=" << TTCluster->channelID()
		  << ": strip=" << TTCluster->strip() 
		  << ": sector=" << TTCluster->sector() 
		  << ": detRegion=" << TTCluster->detRegion() 
		  << ": layer=" << TTCluster->layer() 
		  << ": station=" << TTCluster->station() 
		  << endmsg;
	}
      }else{
	debug () << "no ttLink" << endmsg;
	}*/
      if (!ttLink.notFound()){
	debug() << "begin loop over reco downstream tracks" << endmsg;
	for( LHCb::Tracks::const_iterator it = downstream_tracks->begin(); downstream_tracks->end() != it; ++it ) {
	  //const Track* tr = *it;
	  //debug () << "looping over lhcbIDs" << endmsg;
	  int IDnumber = 0;
	  for( std::vector<LHCb::LHCbID>::const_iterator iId = (*it)->lhcbIDs().begin();
	       (*it)->lhcbIDs().end() != iId; ++iId ){
	    IDnumber++;
	    //debug() << "id " << IDnumber << " = " << *iId << endmsg;
	    int TrackSTID((*iId).lhcbID());
	    if (!ttLink.notFound()){
	      debug() << "now loop over UT clusters for this MC particle" << endmsg;
	      int nhits = 0;
	      const LHCb::STCluster* TTCluster = ttLink.first((*ip));
	      for ( ; 0 != TTCluster; TTCluster = ttLink.next()) {
		//if ( !TTCluster->isUT() ) continue;
		//LHCbID : { UT STChannelID : 20128867 : type=2 strip=99 sector=73 detRegion=1 layer=2 station=1 } }
		nhits++;
		if(TTCluster->channelID() == TrackSTID){
		  debug() << "matching hit " << nhits 
			  << ": type=" << TTCluster->channelID() 
			  << ": strip=" << TTCluster->strip() 
			  << ": sector=" << TTCluster->sector() 
			  << ": detRegion=" << TTCluster->detRegion() 
			  << ": layer=" << TTCluster->layer() 
			  << ": station=" << TTCluster->station() 
			  << endmsg;
		}else{
		  debug() << "STID = " << TrackSTID << ", MC ID = " << TTCluster->channelID()  << endmsg;
		}
	      }
	    }else{
	      debug () << "no ttLink" << endmsg;
	    }
	  }
	}
	//break;
      }else{
	debug() << "no TT link found" << endmsg;
      }
      StatusCode sc = ntuple->write();
      if (sc.isFailure()) return Error("Failed to fill the ntuple", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
      
    }//end of loop over MC particles
  
  return StatusCode::SUCCESS;
}
//=============================================================================
// Loop over tracks and look for matched MC particles
//=============================================================================
StatusCode PrDownstreamChecker::loopOverTracks()
{
  
  debug()   << "Loopng over reconstructed downstream tracks" << endmsg ;
  // book the ntuple
  Tuples::Tuple ntuple = 
    GaudiTupleAlg::nTuple("DownstreamTracks","DownstreamTracks", CLID_ColumnWiseTuple );
    
  // Loop over tracks
  LHCb::Tracks* tracks = get<LHCb::Tracks>(m_DownstreamContainer);

  int nTrack = -1;
  n_in_evt = -1;
  for (LHCb::Tracks::const_iterator it = tracks->begin(); it!=tracks->end(); ++it)
    {
      n_in_evt++;
      nTrack++;
      //debug()   << "processing reconstructed downstream track #" 
      //<< nTrack << "/" << tracks->size() <<  endmsg ;
      fillTrack(*it,ntuple);


      //add n tracks info. AD 1-24-14

      int ngoodx1(-2), ngoodx2(-2),ngoodu(-2),ngoodv(-2);
      int ntotx1(-2), ntotx2(-2),ntotu(-2),ntotv(-2);
      ngoodx1 = (*it)->info(1900,-1);//good 
      ngoodu = (*it)->info(1901,-1);//good 
      ngoodv = (*it)->info(1902,-1);//good 
      ngoodx2 = (*it)->info(1903,-1);//good 
      ntotx1 = (*it)->info(1904,-1);//tot 
      ntotu = (*it)->info(1905,-1);//tot 
      ntotv = (*it)->info(1906,-1);//tot 
      ntotx2 = (*it)->info(1907,-1);//tot 
      
      //initialized. add it to the tuple
      ntuple->column("Goodx1Hits",ngoodx1);
      ntuple->column("GooduHits",ngoodu);
      ntuple->column("GoodvHits",ngoodv);
      ntuple->column("Goodx2Hits",ngoodx2);

      ntuple->column("Totx1Hits",ntotx1);
      ntuple->column("TotuHits",ntotu);
      ntuple->column("TotvHits",ntotv);
      ntuple->column("Totx2Hits",ntotx2);

      AsctTool::DirectType::Range range = m_directTable->relations(*it);
      bool ghost = range.empty();
      
      const LHCb::MCParticle* particle(0);
      if (!ghost) {
	particle = range.begin()->to();
	if( particle && particle->particleID().threeCharge()==0 ) {
	  particle = 0;
	}
      }
      
      fillMCParticle(particle,ntuple);
      
      /*     bool reconstructible_asLong(false);
      reconstructible_asLong = 
	m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,particle);
      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      bool reconstructible_asDown(false);
      reconstructible_asDown =
        m_selector->isReconstructibleAs(IMCReconstructible::ChargedDownstream,particle);
      ntuple->column("reconstructible_asDown",reconstructible_asDown);//AD 1-10-14
      */     
      LHCb::Track* seedTr = new LHCb::Track;
      LHCb::Track* seedUTTr = *it;
      
      SmartRefVector<LHCb::Track>& ancestor = seedUTTr->ancestors();
      for( SmartRefVector<LHCb::Track>::iterator trIt = ancestor.begin();
	   ancestor.end() != trIt; trIt++) {
	seedTr = *trIt;
      }
      
      AsctTool::DirectType::Range range_seed = m_directTable_seed->relations(seedTr);
      bool ghost_seed = range_seed.empty();
      ntuple->column("ghost_seed", ghost_seed);
      

            
      std::vector<double>new_chi2s;new_chi2s.reserve(100);
      std::vector<double>new_ndfs;new_ndfs.reserve(100);
      std::vector<int>new_chi2_is_mc_cand;new_chi2_is_mc_cand.reserve(100);
      new_chi2s.clear();
      new_ndfs.clear();
      new_chi2_is_mc_cand.clear();

      
      int size_of_track_candidates = (*it)->info(36000,-999.);
      if(size_of_track_candidates >100)
      {
        size_of_track_candidates=100;}
      
      debug()<<"size_of_track_candidates = "<<size_of_track_candidates<<endmsg;
      ntuple->column("n_track_candidates_considered",size_of_track_candidates);
      
      for(int i=1; i<size_of_track_candidates;++i){
        debug()<<"i = "<<i<<endmsg;
        debug()<<"info 36000 + i = "<<(*it)->info(36000+i,-9999.)<<endmsg;
        if(i>100)continue;
        new_chi2s.push_back((*it)->info(36000+i,-9999.));
        new_ndfs.push_back((*it)->info(36500+i,-9999.));
        new_chi2_is_mc_cand.push_back((*it)->info(37000+i,-1));
        //new_chi2_is_matched_2_store.push_back((*it)->info(38000+i,-1));
        
        debug()<<"new_chi2s[i] = "<<new_chi2s[i]<<endmsg;
      }
      debug()<<"makking array"<<endmsg;
      ntuple->farray("new_chi2",new_chi2s.begin(), new_chi2s.end(),"size_of_track_candidates",100);
      ntuple->farray("new_ndf",new_ndfs.begin(), new_ndfs.end(),"size_of_track_candidates",100);
      ntuple->farray("new_chi2_is_good_mc_cand",new_chi2_is_mc_cand.begin(),new_chi2_is_mc_cand.end(),
                     "size_of_track_candidates",100);
      //ntuple->farray("new_chi2_is_final_track",new_chi2_is_matched_2_store.begin(), new_chi2_is_matched_2_store.end(),
      //             "size_of_track_candidates",100);
      debug()<<"array done"<<endmsg;
      ntuple->column("saved_track_new_chi2", (*it)->info(38000,-1.) );
      ntuple->column("saved_track_new_ndf", (*it)->info(39000,-1.) );

      ntuple->column("ft_begin_true_xmid",(*it)->info(99000,-99999.));
      ntuple->column("ft_begin_true_ymid",(*it)->info(99001,-99999.));
      ntuple->column("ft_begin_true_zmid",(*it)->info(99002,-99999.));
      ntuple->column("ft_begin_true_tx",(*it)->info(99003,-99999.));
      ntuple->column("ft_begin_true_ty",(*it)->info(99004,-99999.));
      ntuple->column("ft_begin_true_qop",(*it)->info(99005,-99999.));
      ntuple->column("ft_begin_true_entryx",(*it)->info(99006,-99999.));
      ntuple->column("ft_begin_true_entryy",(*it)->info(99007,-99999.));
      ntuple->column("ft_begin_true_entryz",(*it)->info(99008,-99999.));

      StatusCode sc = ntuple->write();
      if (sc.isFailure()) return Error("Failed to fill the ntuple", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
    
      
      
    }//end of loop over tracks
  
  return StatusCode::SUCCESS;
}
//=============================================================================
// Loop over tracks and look for matched MC particles
//=============================================================================
StatusCode PrDownstreamChecker::loopOverSeeds()
{
  // book the ntuple
  Tuples::Tuple ntuple = GaudiTupleAlg::nTuple("SeedTracks","SeedTracks", CLID_ColumnWiseTuple );
    
  // Loop over tracks
  LHCb::Tracks* tracks = get<LHCb::Tracks>(m_TSeedContainer);
  
  n_in_evt = -1;
  for (LHCb::Tracks::const_iterator it = tracks->begin(); it!=tracks->end(); ++it)
    {
      n_in_evt++;
      fillTrack(*it,ntuple);
    
      //m_directTable_seed 
      AsctTool::DirectType::Range range = m_directTable_seed->relations(*it);
      bool ghost = range.empty();
      
      const LHCb::MCParticle* particle(0);
      if (!ghost) {
	particle = range.begin()->to();
	if( particle && particle->particleID().threeCharge()==0 ) {
	  particle = 0;
	}
      }
      
      fillMCParticle(particle,ntuple);
      //hack the ideal states myself
      
      bool reconstructible_asLong(false);
      reconstructible_asLong = 
	m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,particle);
      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      
      LHCb::Track* seedTr = new LHCb::Track;
      LHCb::Track* seedUTTr = *it;
      
      SmartRefVector<LHCb::Track>& ancestor = seedUTTr->ancestors();
      for( SmartRefVector<LHCb::Track>::iterator trIt = ancestor.begin();
	   ancestor.end() != trIt; trIt++) {
	seedTr = *trIt;
      }
      
      AsctTool::DirectType::Range range_seed = m_directTable_seed->relations(seedTr);
      bool ghost_seed = range_seed.empty();
      ntuple->column("ghost_seed", ghost_seed);
      
      StatusCode sc = ntuple->write();
      if (sc.isFailure()) return Error("Failed to fill the ntuple", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
      
    }//end of loop over tracks
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Loop over tracks and look for matched MC particles
//=============================================================================
StatusCode PrDownstreamChecker::loopOverUpdatedSeeds()
{
  // book the ntuple
  Tuples::Tuple ntuple = GaudiTupleAlg::nTuple("UpdatedSeedTracks","UpdatedSeedTracks", CLID_ColumnWiseTuple );
    
  // Loop over tracks
  LHCb::Tracks* tracks = get<LHCb::Tracks>(m_Test_Container);
  
  n_in_evt = -1;

  //stuff for hit info
      std::vector<double>x1_linehit_mc_x,x1_linehit_mc_y;
      std::vector<double>presel_track_at_x1_mc_x, presel_track_at_x1_mc_y;
      std::vector<double>final_track_at_x1_mc_x, final_track_at_x1_mc_y;
      std::vector<double>x1_truth_mc_x, x1_truth_mc_y,x1_truth_mc_z;
      
      std::vector<double>u_linehit_mc_x,u_linehit_mc_y;
      std::vector<double>presel_track_at_u_mc_x, presel_track_at_u_mc_y;
      std::vector<double>final_track_at_u_mc_x, final_track_at_u_mc_y;
      std::vector<double>u_truth_mc_x, u_truth_mc_y,u_truth_mc_z;
      
      std::vector<double>v_linehit_mc_x,v_linehit_mc_y;
      std::vector<double>presel_track_at_v_mc_x, presel_track_at_v_mc_y;
      std::vector<double>final_track_at_v_mc_x, final_track_at_v_mc_y;
      std::vector<double>v_truth_mc_x, v_truth_mc_y,v_truth_mc_z;
      
      std::vector<double>x2_linehit_mc_x,x2_linehit_mc_y;
      std::vector<double>presel_track_at_x2_mc_x, presel_track_at_x2_mc_y;
      std::vector<double>final_track_at_x2_mc_x, final_track_at_x2_mc_y;
      std::vector<double>x2_truth_mc_x, x2_truth_mc_y,x2_truth_mc_z;
      //final track candid hits
      std::vector<double>x1_linehit_final_candid_x,x1_linehit_final_candid_y;
      std::vector<double>presel_track_at_x1_final_candid_x, presel_track_at_x1_final_candid_y;
      std::vector<double>final_track_at_x1_final_candid_x, final_track_at_x1_final_candid_y;
      std::vector<double>x1_truth_final_candid_x, x1_truth_final_candid_y;
      
      std::vector<double>u_linehit_final_candid_x,u_linehit_final_candid_y;
      std::vector<double>presel_track_at_u_final_candid_x, presel_track_at_u_final_candid_y;
      std::vector<double>final_track_at_u_final_candid_x, final_track_at_u_final_candid_y;
      std::vector<double>u_truth_final_candid_x, u_truth_final_candid_y;
      
      std::vector<double>v_linehit_final_candid_x,v_linehit_final_candid_y;
      std::vector<double>presel_track_at_v_final_candid_x, presel_track_at_v_final_candid_y;
      std::vector<double>final_track_at_v_final_candid_x, final_track_at_v_final_candid_y;
      std::vector<double>v_truth_final_candid_x, v_truth_final_candid_y;
      
      std::vector<double>x2_linehit_final_candid_x,x2_linehit_final_candid_y;
      std::vector<double>presel_track_at_x2_final_candid_x, presel_track_at_x2_final_candid_y;
      std::vector<double>final_track_at_x2_final_candid_x, final_track_at_x2_final_candid_y;
      std::vector<double>x2_truth_final_candid_x, x2_truth_final_candid_y;
      //hacked state params
      std::vector<double>hacked_state_x, hacked_state_y, hacked_state_tx,hacked_state_ty;//,hacked_state_qop;
      hacked_state_x.reserve(20);hacked_state_y.reserve(20);
      hacked_state_tx.reserve(20);hacked_state_ty.reserve(20);
      //hack zpostUT
      std::vector<double>hacked_state_zpUT_x, hacked_state_zpUT_y, hacked_state_zpUT_tx,hacked_state_zpUT_ty;//,hacked_state_zpUT_qop;
      hacked_state_zpUT_x.reserve(20);hacked_state_zpUT_y.reserve(20);
      hacked_state_zpUT_tx.reserve(20);hacked_state_zpUT_ty.reserve(20);
      //hack endvelo
      std::vector<double>hacked_EndVelo_state_x, hacked_EndVelo_state_y, hacked_EndVelo_state_tx,hacked_EndVelo_state_ty;//,hacked_EndVelo_state_qop;
      hacked_EndVelo_state_x.reserve(20);hacked_EndVelo_state_y.reserve(20);
      hacked_EndVelo_state_tx.reserve(20);hacked_EndVelo_state_ty.reserve(20);

      // hack midTT
      std::vector<double>hacked_midTT_state_x, hacked_midTT_state_y, hacked_midTT_state_tx,hacked_midTT_state_ty;//,hacked_midTT_state_qop;
      hacked_midTT_state_x.reserve(20);hacked_midTT_state_y.reserve(20);
      hacked_midTT_state_tx.reserve(20);hacked_midTT_state_ty.reserve(20);
      //begin T
      std::vector<double>hacked_begT_state_x, hacked_begT_state_y, hacked_begT_state_tx,hacked_begT_state_ty;//,hacked_begT_state_qop;
      hacked_begT_state_x.reserve(20);hacked_begT_state_y.reserve(20);
      hacked_begT_state_tx.reserve(20);hacked_begT_state_ty.reserve(20);

      x1_linehit_mc_x.reserve(20);x1_linehit_mc_y.reserve(20);
      presel_track_at_x1_mc_x.reserve(20); presel_track_at_x1_mc_y.reserve(20);
      final_track_at_x1_mc_x.reserve(20); final_track_at_x1_mc_y.reserve(20);
      x1_truth_mc_x.reserve(20); x1_truth_mc_y.reserve(20);x1_truth_mc_z.reserve(20);
      
      u_linehit_mc_x.reserve(20);u_linehit_mc_y.reserve(20);
      presel_track_at_u_mc_x.reserve(20); presel_track_at_u_mc_y.reserve(20);
      final_track_at_u_mc_x.reserve(20); final_track_at_u_mc_y.reserve(20);
      u_truth_mc_x.reserve(20); u_truth_mc_y.reserve(20);u_truth_mc_z.reserve(20);
      
      v_linehit_mc_x.reserve(20);v_linehit_mc_y.reserve(20);
      presel_track_at_v_mc_x.reserve(20); presel_track_at_v_mc_y.reserve(20);
      final_track_at_v_mc_x.reserve(20); final_track_at_v_mc_y.reserve(20);
      v_truth_mc_x.reserve(20); v_truth_mc_y.reserve(20);v_truth_mc_z.reserve(20);
      
      x2_linehit_mc_x.reserve(20);x2_linehit_mc_y.reserve(20);
      presel_track_at_x2_mc_x.reserve(20); presel_track_at_x2_mc_y.reserve(20);
      final_track_at_x2_mc_x.reserve(20); final_track_at_x2_mc_y.reserve(20);
      x2_truth_mc_x.reserve(20); x2_truth_mc_y.reserve(20);x2_truth_mc_z.reserve(20);
      //final track candid hits
      x1_linehit_final_candid_x.reserve(20);x1_linehit_final_candid_y.reserve(20);
      presel_track_at_x1_final_candid_x.reserve(20); presel_track_at_x1_final_candid_y.reserve(20);
      final_track_at_x1_final_candid_x.reserve(20); final_track_at_x1_final_candid_y.reserve(20);
      x1_truth_final_candid_x.reserve(20); x1_truth_final_candid_y.reserve(20);
      
      u_linehit_final_candid_x.reserve(20);u_linehit_final_candid_y.reserve(20);
      presel_track_at_u_final_candid_x.reserve(20); presel_track_at_u_final_candid_y.reserve(20);
      final_track_at_u_final_candid_x.reserve(20); final_track_at_u_final_candid_y.reserve(20);
      u_truth_final_candid_x.reserve(20); u_truth_final_candid_y.reserve(20);
      
      v_linehit_final_candid_x.reserve(20);v_linehit_final_candid_y.reserve(20);
      presel_track_at_v_final_candid_x.reserve(20); presel_track_at_v_final_candid_y.reserve(20);
      final_track_at_v_final_candid_x.reserve(20); final_track_at_v_final_candid_y.reserve(20);
      v_truth_final_candid_x.reserve(20); v_truth_final_candid_y.reserve(20);
      
      x2_linehit_final_candid_x.reserve(20);x2_linehit_final_candid_y.reserve(20);
      presel_track_at_x2_final_candid_x.reserve(20); presel_track_at_x2_final_candid_y.reserve(20);
      final_track_at_x2_final_candid_x.reserve(20); final_track_at_x2_final_candid_y.reserve(20);
      x2_truth_final_candid_x.reserve(20); x2_truth_final_candid_y.reserve(20);
      
      std::vector<double>new_chi2s;new_chi2s.reserve(100);
      std::vector<double>new_ndfs;new_ndfs.reserve(100);
      std::vector<int>new_chi2_is_mc_cand;new_chi2_is_mc_cand.reserve(100);
      //std::vector<int>new_chi2_is_matched_2_store;new_chi2_is_matched_2_store.reserve(100);
      
  for (LHCb::Tracks::const_iterator it = tracks->begin(); it!=tracks->end(); ++it)
    {
      n_in_evt++;
      fillTrack(*it,ntuple);
      /*
      //add nhits
      int ngood_before_x1(-2), ngood_before_u(-2),ngood_before_v(-2),ngood_before_x2(-2);
      int ntot_before_x1(-2), ntot_before_u(-2),ntot_before_v(-2),ntot_before_x2(-2);
      int ngood_after_x1(-2), ngood_after_u(-2),ngood_after_v(-2),ngood_after_x2(-2);
      int ntot_after_x1(-2), ntot_after_u(-2),ntot_after_v(-2),ntot_after_x2(-2);
      //end info added by AD
      
      ngood_before_x1=(*it)->info(1810,-1);
      ngood_before_u=(*it)->info(1811,-1);
      ngood_before_v=(*it)->info(1812,-1);
      ngood_before_x2=(*it)->info(1813,-1);
      
      ntot_before_x1=(*it)->info(1820,-1);
      ntot_before_u=(*it)->info(1821,-1);
      ntot_before_v=(*it)->info(1822,-1);
      ntot_before_x2=(*it)->info(1823,-1);
      
      ngood_after_x1=(*it)->info(1830,-1);
      ngood_after_u=(*it)->info(1831,-1);
      ngood_after_v=(*it)->info(1832,-1);
      ngood_after_x2=(*it)->info(1833,-1);
      
      ntot_after_x1=(*it)->info(1840,-1);
      ntot_after_u=(*it)->info(1841,-1);
      ntot_after_v=(*it)->info(1842,-1);
      ntot_after_x2=(*it)->info(1843,-1);
      debug()<<"added info val in prdownstreamchecker :"<<endmsg;
      debug()<<"ngood_before_x1 = "<<ngood_before_x1<<endmsg;
   debug()<<"ngood_before_u = "<<ngood_before_u<<endmsg;
   debug()<<"ngood_before_v = "<<ngood_before_v<<endmsg;
   debug()<<"ngood_before_x2 = "<<ngood_before_x2<<endmsg;

   debug()<<"ntot_before_x1 = "<<ntot_before_x1<<endmsg;
   debug()<<"ntot_before_u = "<<ntot_before_u<<endmsg;
   debug()<<"ntot_before_v = "<<ntot_before_v<<endmsg;
   debug()<<"ntot_before_x2 = "<<ntot_before_x2<<endmsg;

   debug()<<"ngood_after_x1 = "<<ngood_after_x1<<endmsg;
   debug()<<"ngood_after_u = "<<ngood_after_u<<endmsg;
   debug()<<"ngood_after_v = "<<ngood_after_v<<endmsg;
   debug()<<"ngood_after_x2 = "<<ngood_after_x2<<endmsg;

   debug()<<"ntot_after_x1 = "<<ntot_after_x1<<endmsg;
   debug()<<"ntot_after_u = "<<ntot_after_u<<endmsg;
   debug()<<"ntot_after_v = "<<ntot_after_v<<endmsg;
   debug()<<"ntot_after_x2 = "<<ntot_after_x2<<endmsg;

   

      ntuple->column("beforeFAR_good_x1",ngood_before_x1);
      ntuple->column("beforeFAR_good_u",ngood_before_u);
      ntuple->column("beforeFAR_good_v",ngood_before_v);
      ntuple->column("beforeFAR_good_x2",ngood_before_x2);

      ntuple->column("beforeFAR_total_x1",ntot_before_x1);
      ntuple->column("beforeFAR_total_u",ntot_before_u);
      ntuple->column("beforeFAR_total_v",ntot_before_v);
      ntuple->column("beforeFAR_total_x2",ntot_before_x2);

      ntuple->column("afterFAR_good_x1",ngood_after_x1);
      ntuple->column("afterFAR_good_u",ngood_after_u);
      ntuple->column("afterFAR_good_v",ngood_after_v);
      ntuple->column("afterFAR_good_x2",ngood_after_x2);

      ntuple->column("afterFAR_total_x1",ntot_after_x1);
      ntuple->column("afterFAR_total_u",ntot_after_u);
      ntuple->column("afterFAR_total_v",ntot_after_v);
      ntuple->column("afterFAR_total_x2",ntot_after_x2);
    */
      //added AD 1-29-14
      ntuple->column("n_DSCandidates_per_seed",(*it)->info(2000,-1));
      ntuple->column("presel_count",(*it)->info(2001,-1));
      ntuple->column("candid_nhits",(*it)->info(2002,-1));
      ntuple->column("n_Replaced_3hits",(*it)->info(2003,-1));
      ntuple->column("n_Replaced_4hits",(*it)->info(2004,-1));
      ntuple->column("n_Replaced_5ormore_hits",(*it)->info(2005,-1));
      ntuple->column("n_mc_clones",(*it)->info(2006,-1));

      ntuple->column("mean_ghost_chi2",(*it)->info(2007,-1.) );
      ntuple->column("best_ghost_chi2",(*it)->info(2008,-1.) );
      ntuple->column("n_hits_best_ghost", (*it)->info(2009,-1) );
      ntuple->column("best_mc_chi2", (*it)->info(2010,-1.) );
      
      ntuple->column("presel_forward_match_x", (*it)->info(2011,-1.) );
      ntuple->column("presel_forward_match_uv", (*it)->info(2012,-1.) );
      ntuple->column("candid_forward_match_x", (*it)->info(2013,-1.) );
      ntuple->column("candid_forward_match_uv", (*it)->info(2014,-1.) );

      //m_directTable_seed 
      AsctTool::DirectType::Range range = m_directTable_updatedseed->relations(*it);
      bool ghost = range.empty();
      
      const LHCb::MCParticle* particle(0);
      if (!ghost) {
	particle = range.begin()->to();
	if( particle && particle->particleID().threeCharge()==0 ) {
	  particle = 0;
	}
      }
      /*
      MCTrackInfo trackInfo( evtSvc(), msgSvc() );
        
      // information for the tuple
      //      bool reconstructible_asLong(false);
      bool B_child(false);
      bool D_child(false);
      bool isLong(false);
      bool isDown(false);
      bool isInVelo(false);
      bool over5(false);
      int parentABSID = -999;//add to MC partcilefill. AD 1-12-14
      //int nHits = -999;
      int has_all_info(-1);
      
      if(particle){
        //need to make sure we have full track info.
        has_all_info = 
        trackInfo.fullInfo( particle );  */
        /*reconstructible_asLong = 
          m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,particle);*/
      /*        isLong = 
          trackInfo.hasVeloAndT( particle );
        isDown = 
          trackInfo.hasT( particle ) &&  trackInfo.hasTT( particle );
        over5 = 
          5000. < fabs( (particle)->p() );
        isInVelo = 
          trackInfo.hasVelo( particle );
        if ( 0 != (particle)->originVertex() ) {
          const LHCb::MCParticle* mother =  (particle)->originVertex()->mother();
          if ( 0 != mother ) {
            if ( 0 != mother->originVertex() ) {
              //double rOrigin = mother->originVertex()->position().rho();
              parentABSID = abs( mother->particleID().pid() );
            }
          }
        }
        
        
        if(bAncestor(particle)){
          B_child = true;
        }
        if(dAncestor(particle)){
          D_child = true;
        }
        
      }
      
      /// good candidate
      
      
      /// fill the tuple with MC information
      ntuple->column("parentABSID", parentABSID);
      ntuple->column("fullInfo", has_all_info);
      //      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      ntuple->column("isLong", isLong);
      ntuple->column("isDown", isDown);
      ntuple->column("over5", over5);
      ntuple->column("isInVelo", isInVelo);
      ntuple->column("B_child", B_child);
      ntuple->column("D_child", D_child);
      */
      fillMCParticle(particle,ntuple);
      /*    
      bool reconstructible_asLong(false);
      reconstructible_asLong = 
	m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,particle);
      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      */ ///this is now in the fill MC particle. AD 1-15-14
      bool reconstructible_asDown(false);
      reconstructible_asDown =
        m_selector->isReconstructibleAs(IMCReconstructible::ChargedDownstream,particle);
      ntuple->column("reconstructible_asDown",reconstructible_asDown);//AD 1-10-14
      
      LHCb::Track* seedTr = new LHCb::Track;
      LHCb::Track* seedUTTr = *it;
      
      SmartRefVector<LHCb::Track>& ancestor = seedUTTr->ancestors();
      for( SmartRefVector<LHCb::Track>::iterator trIt = ancestor.begin();
	   ancestor.end() != trIt; trIt++) {
	seedTr = *trIt;
      }
      
      AsctTool::DirectType::Range range_seed = m_directTable_updatedseed->relations(seedTr);
      bool ghost_seed = range_seed.empty();
      ntuple->column("ghost_seed", ghost_seed);

      //add reconstructible seed
      AsctTool::InverseType::Range range_inv = m_inverseTable_updatedseed->relations(particle);
      bool reconstructed_seed = !(range_inv.empty());
      ntuple->column("reconstructed_seed",reconstructed_seed);

      //get all the info we added about the residuals
      //ntuple->column("n_x1_mc_hits", (*it)->info(2100,-999.));
      int nx1_mc = (*it)->info(21000,-999);
      int nu_mc = (*it)->info(22000,-999);
      int nv_mc = (*it)->info(23000,-999);
      int nx2_mc = (*it)->info(24000,-999);

      int nx1_final_candid = (*it)->info(31000,-999);
      int nu_final_candid = (*it)->info(32000,-999);
      int nv_final_candid = (*it)->info(33000,-999);
      int nx2_final_candid = (*it)->info(34000,-999);
      
      //initialize all the std::vectors we need.
      //mc hits
      hacked_state_x.clear();hacked_state_y.clear();
      hacked_state_tx.clear();hacked_state_ty.clear();
      
      hacked_state_zpUT_x.clear();hacked_state_zpUT_y.clear();
      hacked_state_zpUT_tx.clear();hacked_state_zpUT_ty.clear();

      hacked_EndVelo_state_x.clear();hacked_EndVelo_state_y.clear();
      hacked_EndVelo_state_tx.clear();hacked_EndVelo_state_ty.clear();

      hacked_midTT_state_x.clear();hacked_midTT_state_y.clear();
      hacked_midTT_state_tx.clear();hacked_midTT_state_ty.clear();

      hacked_begT_state_x.clear();hacked_begT_state_y.clear();
      hacked_begT_state_tx.clear();hacked_begT_state_ty.clear();

      
      x1_linehit_mc_x.clear();x1_linehit_mc_y.clear();
      presel_track_at_x1_mc_x.clear();presel_track_at_x1_mc_y.clear();
      final_track_at_x1_mc_x.clear(); final_track_at_x1_mc_y.clear();
      x1_truth_mc_x.clear(); x1_truth_mc_y.clear();x1_truth_mc_z.clear();
      
      u_linehit_mc_x.clear();u_linehit_mc_y.clear();
      presel_track_at_u_mc_x.clear();presel_track_at_u_mc_y.clear();
      final_track_at_u_mc_x.clear();final_track_at_u_mc_y.clear();
      u_truth_mc_x.clear();u_truth_mc_y.clear();u_truth_mc_z.clear();
      
      v_linehit_mc_x.clear();v_linehit_mc_y.clear();
      presel_track_at_v_mc_x.clear();presel_track_at_v_mc_y.clear();
      final_track_at_v_mc_x.clear();final_track_at_v_mc_y.clear();
      v_truth_mc_x.clear();v_truth_mc_y.clear();v_truth_mc_z.clear();
      
      x2_linehit_mc_x.clear();x2_linehit_mc_y.clear();
      presel_track_at_x2_mc_x.clear();presel_track_at_x2_mc_y.clear();
      final_track_at_x2_mc_x.clear();final_track_at_x2_mc_y.clear();
      x2_truth_mc_x.clear();x2_truth_mc_y.clear();x2_truth_mc_z.clear();
      //final track candid hits
      x1_linehit_final_candid_x.clear();x1_linehit_final_candid_y.clear();
      presel_track_at_x1_final_candid_x.clear();presel_track_at_x1_final_candid_y.clear();
      final_track_at_x1_final_candid_x.clear();final_track_at_x1_final_candid_y.clear();
      x1_truth_final_candid_x.clear();x1_truth_final_candid_y.clear();
      
      u_linehit_final_candid_x.clear();u_linehit_final_candid_y.clear();
      presel_track_at_u_final_candid_x.clear();presel_track_at_u_final_candid_y.clear();
      final_track_at_u_final_candid_x.clear();final_track_at_u_final_candid_y.clear();
      u_truth_final_candid_x.clear();u_truth_final_candid_y.clear();
      
      v_linehit_final_candid_x.clear();v_linehit_final_candid_y.clear();
      presel_track_at_v_final_candid_x.clear(); presel_track_at_v_final_candid_y.clear();
      final_track_at_v_final_candid_x.clear();final_track_at_v_final_candid_y.clear();
      v_truth_final_candid_x.clear();v_truth_final_candid_y.clear();
      
      x2_linehit_final_candid_x.clear();x2_linehit_final_candid_y.clear();
      presel_track_at_x2_final_candid_x.clear();presel_track_at_x2_final_candid_y.clear();
      final_track_at_x2_final_candid_x.clear();final_track_at_x2_final_candid_y.clear();
      x2_truth_final_candid_x.clear();x2_truth_final_candid_y.clear();
      new_chi2s.clear();
      new_ndfs.clear();
      new_chi2_is_mc_cand.clear();
      //new_chi2_is_matched_2_store.clear();
      
      
      /*
      std::vector<double>x1_linehit_mc_x,x1_linehit_mc_y;
      std::vector<double>presel_track_at_x1_mc_x, presel_track_at_x1_mc_y;
      std::vector<double>final_track_at_x1_mc_x, final_track_at_x1_mc_y;
      std::vector<double>x1_truth_mc_x, x1_truth_mc_y;
      
      std::vector<double>u_linehit_mc_x,u_linehit_mc_y;
      std::vector<double>presel_track_at_u_mc_x, presel_track_at_u_mc_y;
      std::vector<double>final_track_at_u_mc_x, final_track_at_u_mc_y;
      std::vector<double>u_truth_mc_x, u_truth_mc_y;
      
      std::vector<double>v_linehit_mc_x,v_linehit_mc_y;
      std::vector<double>presel_track_at_v_mc_x, presel_track_at_v_mc_y;
      std::vector<double>final_track_at_v_mc_x, final_track_at_v_mc_y;
      std::vector<double>v_truth_mc_x, v_truth_mc_y;
      
      std::vector<double>x2_linehit_mc_x,x2_linehit_mc_y;
      std::vector<double>presel_track_at_x2_mc_x, presel_track_at_x2_mc_y;
      std::vector<double>final_track_at_x2_mc_x, final_track_at_x2_mc_y;
      std::vector<double>x2_truth_mc_x, x2_truth_mc_y;
      //final track candid hits
      std::vector<double>x1_linehit_final_candid_x,x1_linehit_final_candid_y;
      std::vector<double>presel_track_at_x1_final_candid_x, presel_track_at_x1_final_candid_y;
      std::vector<double>final_track_at_x1_final_candid_x, final_track_at_x1_final_candid_y;
      std::vector<double>x1_truth_final_candid_x, x1_truth_final_candid_y;
      
      std::vector<double>u_linehit_final_candid_x,u_linehit_final_candid_y;
      std::vector<double>presel_track_at_u_final_candid_x, presel_track_at_u_final_candid_y;
      std::vector<double>final_track_at_u_final_candid_x, final_track_at_u_final_candid_y;
      std::vector<double>u_truth_final_candid_x, u_truth_final_candid_y;
      
      std::vector<double>v_linehit_final_candid_x,v_linehit_final_candid_y;
      std::vector<double>presel_track_at_v_final_candid_x, presel_track_at_v_final_candid_y;
      std::vector<double>final_track_at_v_final_candid_x, final_track_at_v_final_candid_y;
      std::vector<double>v_truth_final_candid_x, v_truth_final_candid_y;
      
      std::vector<double>x2_linehit_final_candid_x,x2_linehit_final_candid_y;
      std::vector<double>presel_track_at_x2_final_candid_x, presel_track_at_x2_final_candid_y;
      std::vector<double>final_track_at_x2_final_candid_x, final_track_at_x2_final_candid_y;
      std::vector<double>x2_truth_final_candid_x, x2_truth_final_candid_y;
*/      

        
      //get all the extra info about the track and the hits from the extra info
      //AD 5-10-14
      int nhit_info = 21000;
      for(int i=0; i<nx1_mc;++i){
        nhit_info++;
        x1_linehit_mc_x.push_back((*it)->info(nhit_info,-999));
        x1_linehit_mc_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_x1_mc_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_x1_mc_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_x1_mc_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_x1_mc_y.push_back((*it)->info(nhit_info+100,-999));
        x1_truth_mc_x.push_back((*it)->info(nhit_info+120, -999));
        x1_truth_mc_y.push_back((*it)->info(nhit_info+140, -999));
        x1_truth_mc_z.push_back((*it)->info(nhit_info+160, -999));
        //hacked state
        hacked_state_x.push_back((*it)->info(nhit_info+20000, -999));
        hacked_state_y.push_back((*it)->info(nhit_info+20020, -999));
        hacked_state_tx.push_back((*it)->info(nhit_info+20040, -999));
        hacked_state_ty.push_back((*it)->info(nhit_info+20060, -999));

        hacked_state_zpUT_x.push_back((*it)->info(nhit_info+20400, -999));
        hacked_state_zpUT_y.push_back((*it)->info(nhit_info+20420, -999));
        hacked_state_zpUT_tx.push_back((*it)->info(nhit_info+20440, -999));
        hacked_state_zpUT_ty.push_back((*it)->info(nhit_info+20460, -999));
        
        hacked_EndVelo_state_x.push_back((*it)->info(nhit_info+20100, -999));
        hacked_EndVelo_state_y.push_back((*it)->info(nhit_info+20120, -999));
        hacked_EndVelo_state_tx.push_back((*it)->info(nhit_info+20140, -999));
        hacked_EndVelo_state_ty.push_back((*it)->info(nhit_info+20160, -999));

        hacked_midTT_state_x.push_back((*it)->info(nhit_info+20200, -999));
        hacked_midTT_state_y.push_back((*it)->info(nhit_info+20220, -999));
        hacked_midTT_state_tx.push_back((*it)->info(nhit_info+20240, -999));
        hacked_midTT_state_ty.push_back((*it)->info(nhit_info+20260, -999));

        hacked_begT_state_x.push_back((*it)->info(nhit_info+20300, -999));
        hacked_begT_state_y.push_back((*it)->info(nhit_info+20320, -999));
        hacked_begT_state_tx.push_back((*it)->info(nhit_info+20340, -999));
        hacked_begT_state_ty.push_back((*it)->info(nhit_info+20360, -999));


        
        //        hacked_state_qop.push_back((*it)->info(nhit_info+20080, -999));
        /*info()<<"Got hacked state params"<<endmsg
              <<"x   = "<<(*it)->info(nhit_info+20000, -999)<<endmsg
              <<"y   = "<<(*it)->info(nhit_info+20020, -999)<<endmsg
              <<"tx  = "<<(*it)->info(nhit_info+20040, -999)<<endmsg
              <<"ty  = "<<(*it)->info(nhit_info+20060, -999)<<endmsg
              <<"qop = "<<(*it)->info(nhit_info+20080, -999)<<endmsg;
        */
      }
      //u
      nhit_info=22000;
      for(int i=0; i<nu_mc;++i){
        nhit_info++;
        u_linehit_mc_x.push_back((*it)->info(nhit_info,-999));
        u_linehit_mc_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_u_mc_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_u_mc_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_u_mc_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_u_mc_y.push_back((*it)->info(nhit_info+100,-999));
        u_truth_mc_x.push_back((*it)->info(nhit_info+120, -999));
        u_truth_mc_y.push_back((*it)->info(nhit_info+140, -999));
        u_truth_mc_z.push_back((*it)->info(nhit_info+160, -999));
      }
      //v
      nhit_info=23000;
      for(int i=0; i<nv_mc;++i){
        nhit_info++;
        v_linehit_mc_x.push_back((*it)->info(nhit_info,-999));
        v_linehit_mc_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_v_mc_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_v_mc_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_v_mc_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_v_mc_y.push_back((*it)->info(nhit_info+100,-999));
        v_truth_mc_x.push_back((*it)->info(nhit_info+120, -999));
        v_truth_mc_y.push_back((*it)->info(nhit_info+140, -999));
        v_truth_mc_z.push_back((*it)->info(nhit_info+160, -999));
      }
      
      //x2
      nhit_info=24000;
      for(int i=0; i<nx2_mc;++i){
        nhit_info++;
        x2_linehit_mc_x.push_back((*it)->info(nhit_info,-999));
        x2_linehit_mc_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_x2_mc_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_x2_mc_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_x2_mc_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_x2_mc_y.push_back((*it)->info(nhit_info+100,-999));
        x2_truth_mc_x.push_back((*it)->info(nhit_info+120, -999));
        x2_truth_mc_y.push_back((*it)->info(nhit_info+140, -999));
        x2_truth_mc_z.push_back((*it)->info(nhit_info+160, -999));
      }

      //final track
nhit_info = 31000;
      for(int i=0; i<nx1_final_candid;++i){
        nhit_info++;
        x1_linehit_final_candid_x.push_back((*it)->info(nhit_info,-999));
        x1_linehit_final_candid_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_x1_final_candid_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_x1_final_candid_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_x1_final_candid_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_x1_final_candid_y.push_back((*it)->info(nhit_info+100,-999));
        x1_truth_final_candid_x.push_back((*it)->info(nhit_info+120, -999));
        x1_truth_final_candid_y.push_back((*it)->info(nhit_info+140, -999));
      }
      //u
      nhit_info=32000;
      for(int i=0; i<nu_final_candid;++i){
        nhit_info++;
        u_linehit_final_candid_x.push_back((*it)->info(nhit_info,-999));
        u_linehit_final_candid_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_u_final_candid_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_u_final_candid_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_u_final_candid_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_u_final_candid_y.push_back((*it)->info(nhit_info+100,-999));
        u_truth_final_candid_x.push_back((*it)->info(nhit_info+120, -999));
        u_truth_final_candid_y.push_back((*it)->info(nhit_info+140, -999));
      }
      //v
      nhit_info=33000;
      for(int i=0; i<nv_final_candid;++i){
        nhit_info++;
        v_linehit_final_candid_x.push_back((*it)->info(nhit_info,-999));
        v_linehit_final_candid_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_v_final_candid_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_v_final_candid_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_v_final_candid_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_v_final_candid_y.push_back((*it)->info(nhit_info+100,-999));
        v_truth_final_candid_x.push_back((*it)->info(nhit_info+120, -999));
        v_truth_final_candid_y.push_back((*it)->info(nhit_info+140, -999));
      }
      
      //x2
      nhit_info=34000;
      for(int i=0; i<nx2_final_candid;++i){
        nhit_info++;
        x2_linehit_final_candid_x.push_back((*it)->info(nhit_info,-999));
        x2_linehit_final_candid_y.push_back((*it)->info(nhit_info+20,-999));
        presel_track_at_x2_final_candid_x.push_back((*it)->info(nhit_info+40,-999));
        presel_track_at_x2_final_candid_y.push_back((*it)->info(nhit_info+60,-999));
        final_track_at_x2_final_candid_x.push_back((*it)->info(nhit_info+80,-999));
        final_track_at_x2_final_candid_y.push_back((*it)->info(nhit_info+100,-999));
        x2_truth_final_candid_x.push_back((*it)->info(nhit_info+120, -999));
        x2_truth_final_candid_y.push_back((*it)->info(nhit_info+140, -999));
      }



      //tuple x1 MC hits
      ntuple->farray("x1_linehit_mc_x",x1_linehit_mc_x.begin(),x1_linehit_mc_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("x1_linehit_mc_y",x1_linehit_mc_y.begin(),x1_linehit_mc_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("presel_track_at_x1_mc_x",presel_track_at_x1_mc_x.begin(),presel_track_at_x1_mc_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("presel_track_at_x1_mc_y",presel_track_at_x1_mc_y.begin(),presel_track_at_x1_mc_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("final_track_at_x1_mc_x",final_track_at_x1_mc_x.begin(),final_track_at_x1_mc_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("final_track_at_x1_mc_y",final_track_at_x1_mc_y.begin(),final_track_at_x1_mc_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("x1_truth_mc_x",x1_truth_mc_x.begin(),x1_truth_mc_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("x1_truth_mc_y", x1_truth_mc_y.begin(), x1_truth_mc_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("x1_truth_mc_z", x1_truth_mc_z.begin(), x1_truth_mc_z.end(),"n_x1_mc_hits",20);
      //u
      ntuple->farray("u_linehit_mc_x",u_linehit_mc_x.begin(),u_linehit_mc_x.end(),"n_u_mc_hits",20);
      ntuple->farray("u_linehit_mc_y",u_linehit_mc_y.begin(),u_linehit_mc_y.end(),"n_u_mc_hits",20);
      ntuple->farray("presel_track_at_u_mc_x",presel_track_at_u_mc_x.begin(),presel_track_at_u_mc_x.end(),"n_u_mc_hits",20);
      ntuple->farray("presel_track_at_u_mc_y",presel_track_at_u_mc_y.begin(),presel_track_at_u_mc_y.end(),"n_u_mc_hits",20);
      ntuple->farray("final_track_at_u_mc_x",final_track_at_u_mc_x.begin(),final_track_at_u_mc_x.end(),"n_u_mc_hits",20);
      ntuple->farray("final_track_at_u_mc_y",final_track_at_u_mc_y.begin(),final_track_at_u_mc_y.end(),"n_u_mc_hits",20);
      ntuple->farray("u_truth_mc_x",u_truth_mc_x.begin(),u_truth_mc_x.end(),"n_u_mc_hits",20);
      ntuple->farray("u_truth_mc_y", u_truth_mc_y.begin(), u_truth_mc_y.end(),"n_u_mc_hits",20);
      ntuple->farray("u_truth_mc_z", u_truth_mc_z.begin(), u_truth_mc_z.end(),"n_u_mc_hits",20);
      
      //v
      ntuple->farray("v_linehit_mc_x",v_linehit_mc_x.begin(),v_linehit_mc_x.end(),"n_v_mc_hits",20);
      ntuple->farray("v_linehit_mc_y",v_linehit_mc_y.begin(),v_linehit_mc_y.end(),"n_v_mc_hits",20);
      ntuple->farray("presel_track_at_v_mc_x",presel_track_at_v_mc_x.begin(),presel_track_at_v_mc_x.end(),"n_v_mc_hits",20);
      ntuple->farray("presel_track_at_v_mc_y",presel_track_at_v_mc_y.begin(),presel_track_at_v_mc_y.end(),"n_v_mc_hits",20);
      ntuple->farray("final_track_at_v_mc_x",final_track_at_v_mc_x.begin(),final_track_at_v_mc_x.end(),"n_v_mc_hits",20);
      ntuple->farray("final_track_at_v_mc_y",final_track_at_v_mc_y.begin(),final_track_at_v_mc_y.end(),"n_v_mc_hits",20);
      ntuple->farray("v_truth_mc_x",v_truth_mc_x.begin(),v_truth_mc_x.end(),"n_v_mc_hits",20);
      ntuple->farray("v_truth_mc_y", v_truth_mc_y.begin(), v_truth_mc_y.end(),"n_v_mc_hits",20);
      ntuple->farray("v_truth_mc_z", v_truth_mc_z.begin(), v_truth_mc_z.end(),"n_v_mc_hits",20);
      //x2
      
            ntuple->farray("x2_linehit_mc_x",x2_linehit_mc_x.begin(),x2_linehit_mc_x.end(),"n_x2_mc_hits",20);
      ntuple->farray("x2_linehit_mc_y",x2_linehit_mc_y.begin(),x2_linehit_mc_y.end(),"n_x2_mc_hits",20);
      ntuple->farray("presel_track_at_x2_mc_x",presel_track_at_x2_mc_x.begin(),presel_track_at_x2_mc_x.end(),"n_x2_mc_hits",20);
      ntuple->farray("presel_track_at_x2_mc_y",presel_track_at_x2_mc_y.begin(),presel_track_at_x2_mc_y.end(),"n_x2_mc_hits",20);
      ntuple->farray("final_track_at_x2_mc_x",final_track_at_x2_mc_x.begin(),final_track_at_x2_mc_x.end(),"n_x2_mc_hits",20);
      ntuple->farray("final_track_at_x2_mc_y",final_track_at_x2_mc_y.begin(),final_track_at_x2_mc_y.end(),"n_x2_mc_hits",20);
      ntuple->farray("x2_truth_mc_x",x2_truth_mc_x.begin(),x2_truth_mc_x.end(),"n_x2_mc_hits",20);
      ntuple->farray("x2_truth_mc_y", x2_truth_mc_y.begin(), x2_truth_mc_y.end(),"n_x2_mc_hits",20);
      ntuple->farray("x2_truth_mc_z", x2_truth_mc_z.begin(), x2_truth_mc_z.end(),"n_x2_mc_hits",20);

      //tuple final hits
      ntuple->farray("x1_linehit_final_candid_x",x1_linehit_final_candid_x.begin(),x1_linehit_final_candid_x.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("x1_linehit_final_candid_y",x1_linehit_final_candid_y.begin(),x1_linehit_final_candid_y.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("presel_track_at_x1_final_candid_x",presel_track_at_x1_final_candid_x.begin(),presel_track_at_x1_final_candid_x.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("presel_track_at_x1_final_candid_y",presel_track_at_x1_final_candid_y.begin(),presel_track_at_x1_final_candid_y.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("final_track_at_x1_final_candid_x",final_track_at_x1_final_candid_x.begin(),final_track_at_x1_final_candid_x.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("final_track_at_x1_final_candid_y",final_track_at_x1_final_candid_y.begin(),final_track_at_x1_final_candid_y.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("x1_truth_final_candid_x",x1_truth_final_candid_x.begin(),x1_truth_final_candid_x.end(),"n_x1_final_candid_hits",20);
      ntuple->farray("x1_truth_final_candid_y", x1_truth_final_candid_y.begin(), x1_truth_final_candid_y.end(),"n_x1_final_candid_hits",20);
      //u
      ntuple->farray("u_linehit_final_candid_x",u_linehit_final_candid_x.begin(),u_linehit_final_candid_x.end(),"n_u_final_candid_hits",20);
      ntuple->farray("u_linehit_final_candid_y",u_linehit_final_candid_y.begin(),u_linehit_final_candid_y.end(),"n_u_final_candid_hits",20);
      ntuple->farray("presel_track_at_u_final_candid_x",presel_track_at_u_final_candid_x.begin(),presel_track_at_u_final_candid_x.end(),"n_u_final_candid_hits",20);
      ntuple->farray("presel_track_at_u_final_candid_y",presel_track_at_u_final_candid_y.begin(),presel_track_at_u_final_candid_y.end(),"n_u_final_candid_hits",20);
      ntuple->farray("final_track_at_u_final_candid_x",final_track_at_u_final_candid_x.begin(),final_track_at_u_final_candid_x.end(),"n_u_final_candid_hits",20);
      ntuple->farray("final_track_at_u_final_candid_y",final_track_at_u_final_candid_y.begin(),final_track_at_u_final_candid_y.end(),"n_u_final_candid_hits",20);
      ntuple->farray("u_truth_final_candid_x",u_truth_final_candid_x.begin(),u_truth_final_candid_x.end(),"n_u_final_candid_hits",20);
      ntuple->farray("u_truth_final_candid_y", u_truth_final_candid_y.begin(), u_truth_final_candid_y.end(),"n_u_final_candid_hits",20);
      
      //v
      ntuple->farray("v_linehit_final_candid_x",v_linehit_final_candid_x.begin(),v_linehit_final_candid_x.end(),"n_v_final_candid_hits",20);
      ntuple->farray("v_linehit_final_candid_y",v_linehit_final_candid_y.begin(),v_linehit_final_candid_y.end(),"n_v_final_candid_hits",20);
      ntuple->farray("presel_track_at_v_final_candid_x",presel_track_at_v_final_candid_x.begin(),presel_track_at_v_final_candid_x.end(),"n_v_final_candid_hits",20);
      ntuple->farray("presel_track_at_v_final_candid_y",presel_track_at_v_final_candid_y.begin(),presel_track_at_v_final_candid_y.end(),"n_v_final_candid_hits",20);
      ntuple->farray("final_track_at_v_final_candid_x",final_track_at_v_final_candid_x.begin(),final_track_at_v_final_candid_x.end(),"n_v_final_candid_hits",20);
      ntuple->farray("final_track_at_v_final_candid_y",final_track_at_v_final_candid_y.begin(),final_track_at_v_final_candid_y.end(),"n_v_final_candid_hits",20);
      ntuple->farray("v_truth_final_candid_x",v_truth_final_candid_x.begin(),v_truth_final_candid_x.end(),"n_v_final_candid_hits",20);
      ntuple->farray("v_truth_final_candid_y", v_truth_final_candid_y.begin(), v_truth_final_candid_y.end(),"n_v_final_candid_hits",20);
      //x2
      
            ntuple->farray("x2_linehit_final_candid_x",x2_linehit_final_candid_x.begin(),x2_linehit_final_candid_x.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("x2_linehit_final_candid_y",x2_linehit_final_candid_y.begin(),x2_linehit_final_candid_y.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("presel_track_at_x2_final_candid_x",presel_track_at_x2_final_candid_x.begin(),presel_track_at_x2_final_candid_x.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("presel_track_at_x2_final_candid_y",presel_track_at_x2_final_candid_y.begin(),presel_track_at_x2_final_candid_y.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("final_track_at_x2_final_candid_x",final_track_at_x2_final_candid_x.begin(),final_track_at_x2_final_candid_x.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("final_track_at_x2_final_candid_y",final_track_at_x2_final_candid_y.begin(),final_track_at_x2_final_candid_y.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("x2_truth_final_candid_x",x2_truth_final_candid_x.begin(),x2_truth_final_candid_x.end(),"n_x2_final_candid_hits",20);
      ntuple->farray("x2_truth_final_candid_y", x2_truth_final_candid_y.begin(), x2_truth_final_candid_y.end(),"n_x2_final_candid_hits",20);


      //hacked state
      ntuple->farray("hacked_state_x",hacked_state_x.begin(), hacked_state_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_y",hacked_state_y.begin(), hacked_state_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_tx",hacked_state_tx.begin(), hacked_state_tx.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_ty",hacked_state_ty.begin(), hacked_state_ty.end(),"n_x1_mc_hits",20);

      ntuple->farray("hacked_state_zpUT_x",hacked_state_zpUT_x.begin(), hacked_state_zpUT_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_zpUT_y",hacked_state_zpUT_y.begin(), hacked_state_zpUT_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_zpUT_tx",hacked_state_zpUT_tx.begin(), hacked_state_zpUT_tx.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_state_zpUT_ty",hacked_state_zpUT_ty.begin(), hacked_state_zpUT_ty.end(),"n_x1_mc_hits",20);

      ntuple->farray("hacked_EndVelo_state_x",hacked_EndVelo_state_x.begin(), hacked_EndVelo_state_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_EndVelo_state_y",hacked_EndVelo_state_y.begin(), hacked_EndVelo_state_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_EndVelo_state_tx",hacked_EndVelo_state_tx.begin(), hacked_EndVelo_state_tx.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_EndVelo_state_ty",hacked_EndVelo_state_ty.begin(), hacked_EndVelo_state_ty.end(),"n_x1_mc_hits",20);

      ntuple->farray("hacked_midTT_state_x",hacked_midTT_state_x.begin(), hacked_midTT_state_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_midTT_state_y",hacked_midTT_state_y.begin(), hacked_midTT_state_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_midTT_state_tx",hacked_midTT_state_tx.begin(), hacked_midTT_state_tx.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_midTT_state_ty",hacked_midTT_state_ty.begin(), hacked_midTT_state_ty.end(),"n_x1_mc_hits",20);

      ntuple->farray("hacked_begT_state_x",hacked_begT_state_x.begin(), hacked_begT_state_x.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_begT_state_y",hacked_begT_state_y.begin(), hacked_begT_state_y.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_begT_state_tx",hacked_begT_state_tx.begin(), hacked_begT_state_tx.end(),"n_x1_mc_hits",20);
      ntuple->farray("hacked_begT_state_ty",hacked_begT_state_ty.begin(), hacked_begT_state_ty.end(),"n_x1_mc_hits",20);
      //      ntuple->farray("hacked_state_qop",hacked_state_qop.begin(), hacked_state_qop.end(),"n_x1_mc_hits",20);

      //reference state at 2320 mm
      ntuple->column("ref_state_x",(*it)->info(35000,-999.));
      ntuple->column("ref_state_y",(*it)->info(35001,-999.));
      ntuple->column("ref_state_tx",(*it)->info(35002,-999.));
      ntuple->column("ref_state_ty",(*it)->info(35003,-999.));
      ntuple->column("ref_state_qop",(*it)->info(35004,-999.));
      //new chi2 stuff
      debug()<<"starting chi2 check."<<endmsg;
      
      int size_of_track_candidates = (*it)->info(36000,-999.);
      if(size_of_track_candidates >100)
      {
        size_of_track_candidates=100;}
      
      debug()<<"size_of_track_candidates = "<<size_of_track_candidates<<endmsg;
      ntuple->column("n_track_candidates_considered",size_of_track_candidates);
      
      for(int i=1; i<size_of_track_candidates;++i){
        debug()<<"i = "<<i<<endmsg;
        debug()<<"info 36000 + i = "<<(*it)->info(36000+i,-9999.)<<endmsg;
        if(i>100)continue;
        new_chi2s.push_back((*it)->info(36000+i,-9999.));
        new_ndfs.push_back((*it)->info(36500+i,-9999.));
        new_chi2_is_mc_cand.push_back((*it)->info(37000+i,-1));
        //new_chi2_is_matched_2_store.push_back((*it)->info(38000+i,-1));
        
        debug()<<"new_chi2s[i] = "<<new_chi2s[i]<<endmsg;
      }
      debug()<<"makking array"<<endmsg;
      ntuple->farray("new_chi2",new_chi2s.begin(), new_chi2s.end(),"size_of_track_candidates",100);
      ntuple->farray("new_ndf",new_ndfs.begin(), new_ndfs.end(),"size_of_track_candidates",100);
      ntuple->farray("new_chi2_is_good_mc_cand",new_chi2_is_mc_cand.begin(),new_chi2_is_mc_cand.end(),
                     "size_of_track_candidates",100);
      //ntuple->farray("new_chi2_is_final_track",new_chi2_is_matched_2_store.begin(), new_chi2_is_matched_2_store.end(),
      //             "size_of_track_candidates",100);
      debug()<<"array done"<<endmsg;
      ntuple->column("saved_track_new_chi2", (*it)->info(38000,-1.) );
      ntuple->column("saved_track_new_ndf", (*it)->info(39000,-1.) );
      ntuple->column("ft_begin_true_xmid",(*it)->info(99000,-99999.));
      ntuple->column("ft_begin_true_ymid",(*it)->info(99001,-99999.));
      ntuple->column("ft_begin_true_zmid",(*it)->info(99002,-99999.));
      ntuple->column("ft_begin_true_tx",(*it)->info(99003,-99999.));
      ntuple->column("ft_begin_true_ty",(*it)->info(99004,-99999.));
      ntuple->column("ft_begin_true_qop",(*it)->info(99005,-99999.));
      ntuple->column("ft_begin_true_entryx",(*it)->info(99006,-99999.));
      ntuple->column("ft_begin_true_entryy",(*it)->info(99007,-99999.));
      ntuple->column("ft_begin_true_entryz",(*it)->info(99008,-99999.));

      //hack states.
      /*
      double true_endvelo_x(-9999.), true_endvelo_y(-9999.),true_endvelo_tx(-9999.), true_endvelo_ty(-9999.);
      double true_midUT_x(-9999.), true_midUT_y(-9999.),true_midUT_tx(-9999.), true_midUT_ty(-9999.);
      double true_ftstation_x(-9999.), true_ftstation_y(-9999.),true_ftstation_tx(-9999.), true_ftstation_ty(-9999.);
      

      //get states from PRDownstream as before
      //TODO.
      std::cout<<"true_endvelo_x = "<<true_endvelo_x<<std::endl;
      std::cout<<"true_endvelo_y = "<<true_endvelo_y<<std::endl;
      std::cout<<"true_endvelo_tx = "<<true_endvelo_tx<<std::endl;
      std::cout<<"true_endvelo_ty = "<<true_endvelo_ty<<std::endl;
      
      
      std::cout<<"true_midUT_x = "<<true_midUT_x<<std::endl;
      std::cout<<"true_midUT_y = "<<true_midUT_y<<std::endl;
      std::cout<<"true_midUT_tx = "<<true_midUT_tx<<std::endl;
      std::cout<<"true_midUT_ty = "<<true_midUT_ty<<std::endl;
      
      std::cout<<"true_ftstation_x = "<<true_ftstation_x<<std::endl;
      std::cout<<"true_ftstation_y = "<<true_ftstation_y<<std::endl;
      std::cout<<"true_ftstation_tx = "<<true_ftstation_tx<<std::endl;
      std::cout<<"true_ftstation_ty = "<<true_ftstation_ty<<std::endl;
      
      
      ntuple->column("true_endvelo_x",true_endvelo_x);
      ntuple->column("true_endvelo_y",true_endvelo_y);
      ntuple->column("true_endvelo_tx",true_endvelo_tx);
      ntuple->column("true_endvelo_ty",true_endvelo_ty);

      ntuple->column("true_midUT_x",true_midUT_x);
      ntuple->column("true_midUT_y",true_midUT_y);
      ntuple->column("true_midUT_tx",true_midUT_tx);
      ntuple->column("true_midUT_ty",true_midUT_ty);


      ntuple->column("true_ftstation_x",true_ftstation_x);
      ntuple->column("true_ftstation_y",true_ftstation_y);
      ntuple->column("true_ftstation_tx",true_ftstation_tx);
      ntuple->column("true_ftstation_ty",true_ftstation_ty);

      */
      StatusCode sc = ntuple->write();
      if (sc.isFailure()) return Error("Failed to fill the ntuple", /*StatusCode::FAILURE*/StatusCode::SUCCESS);
      
    }//end of loop over tracks
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Fill ntuple with MC particle's information
//=============================================================================
void PrDownstreamChecker::fillMCParticle(const LHCb::MCParticle* p, Tuples::Tuple ntuple)
{
  double true_p(-999.99), true_px(-999.99), true_py(-999.99), true_pz(-999.99),true_pt(-999.99);
  double true_theta(0),true_eta(0),true_phi(0);
  double true_x(-999.99), true_y(-999.99), true_z(-999.99);
  double true_tx(-999.99), true_ty(-999.99);
  bool ghost(true);
  int pdgcode(0), nUT_hits(0), nUT_layers(0);
  double true_qop(0.00);

  /*
  double true_endvelo_x(-9999.), true_endvelo_y(-9999.),true_endvelo_tx(-9999.), true_endvelo_ty(-9999.);
  double true_midUT_x(-9999.), true_midUT_y(-9999.),true_midUT_tx(-9999.), true_midUT_ty(-9999.);
  double true_ftstation_x(-9999.), true_ftstation_y(-9999.),true_ftstation_tx(-9999.), true_ftstation_ty(-9999.);
*/
  //int parentABSID(-999);
      
  // information for the tuple. do it in fillmcparticle from 1-15-14
  bool reconstructible_asLong(false);
  bool reconstructible_asDown(false);
  bool B_child(false);
  bool D_child(false);
  bool isLong(false);
  bool isDown(false);
  bool isInVelo(false);
  bool over5(false);
  int parentABSID = -999;//add to MC partcilefill. AD 1-12-14
  int has_all_info(-1);
  
  MCTrackInfo trackInfos( evtSvc(), msgSvc() );
   
  //int nHits = -999;
  //need to make sure we have full track info.
  
  // mc level quantities
  if (p) {
    ghost = false;
    true_p = p->p();
    true_px = p->momentum().x();
    true_py = p->momentum().y();
    true_pz = p->momentum().z();
    true_pt = p->pt();
    true_theta = p->momentum().theta();
    true_eta = p->momentum().Eta();
    true_phi = p->momentum().Phi();
    true_x = p->originVertex()->position4vector().x();
    true_y = p->originVertex()->position4vector().y();
    true_z = p->originVertex()->position4vector().z();
    true_tx = p->momentum().x()/p->momentum().z();
    true_ty = p->momentum().y()/p->momentum().z();
    pdgcode = p->particleID().pid();
    nUT_hits = getUT_hits(p);
    nUT_layers = getUT_layers(p);
    true_qop = (p->particleID().threeCharge()/3.0)/(p->p());    
    /*
    //fill some additional info on the state
    LHCb::State *state_end_velo,*state_midUT, *state_FT;

    m_stateCreator->createState(p,StateParameters::ZEndVelo,*state_end_velo);
    m_stateCreator->createState(p,StateParameters::ZMidTT,*state_midUT);
    m_stateCreator->createState(p,StateParameters::ZBegT,*state_FT);
    //here. Move this to another function.
    state_end_velo = IdealStateHacker(p,StateParameters::ZEndVelo);
    state_midUT=IdealStateHacker(p,StateParameters::ZMidTT);
    state_FT = IdealStateHacker(p,StateParameters::ZBegT);
    
    if(state_end_velo==NULL){
      debug()<<"no end velo ideal state created"<<endmsg;
    }
    if(state_midUT==NULL){
      debug()<<"no mid UT ideal state created"<<endmsg;
    }
    if(state_FT==NULL){
      debug()<<"no FT ideal state created"<<endmsg;
    }
    
    
    true_endvelo_x =state_end_velo->x();
    true_endvelo_y = state_end_velo->y();
    true_endvelo_tx = state_end_velo->tx();
    true_endvelo_ty = state_end_velo->ty();

    true_midUT_x =state_midUT->x();
    true_midUT_y = state_midUT->y();
    true_midUT_tx = state_midUT->tx();
    true_midUT_ty = state_midUT->ty();

    true_ftstation_x =state_FT->x();
    true_ftstation_y = state_FT->y();
    true_ftstation_tx = state_FT->tx();
    true_ftstation_ty = state_FT->ty();
    
    std::cout<<"true_endvelo_x = "<<true_endvelo_x<<std::endl;
    std::cout<<"true_endvelo_y = "<<true_endvelo_y<<std::endl;
    std::cout<<"true_endvelo_tx = "<<true_endvelo_tx<<std::endl;
    std::cout<<"true_endvelo_ty = "<<true_endvelo_ty<<std::endl;


    std::cout<<"true_midUT_x = "<<true_midUT_x<<std::endl;
    std::cout<<"true_midUT_y = "<<true_midUT_y<<std::endl;
    std::cout<<"true_midUT_tx = "<<true_midUT_tx<<std::endl;
    std::cout<<"true_midUT_ty = "<<true_midUT_ty<<std::endl;

    std::cout<<"true_ftstation_x = "<<true_ftstation_x<<std::endl;
    std::cout<<"true_ftstation_y = "<<true_ftstation_y<<std::endl;
    std::cout<<"true_ftstation_tx = "<<true_ftstation_tx<<std::endl;
    std::cout<<"true_ftstation_ty = "<<true_ftstation_ty<<std::endl;

*/
    
      //need to make sure we have full track info.
      has_all_info = 
        trackInfos.fullInfo( p );

      reconstructible_asLong = 
        m_selector->isReconstructibleAs(IMCReconstructible::ChargedLong,p);

      reconstructible_asDown =
        m_selector->isReconstructibleAs(IMCReconstructible::ChargedDownstream,p);
      
      isLong = 
        trackInfos.hasVeloAndT( p );
      
      isDown = 
        trackInfos.hasT( p ) &&  trackInfos.hasTT( p );
      
      over5 = 
        5000. < fabs( (p)->p() );
      
      isInVelo = 
        trackInfos.hasVelo( p );

      if ( 0 != (p)->originVertex() ) {
        const LHCb::MCParticle* mother =  (p)->originVertex()->mother();
        if ( 0 != mother ) {
          if ( 0 != mother->originVertex() ) {
            //double rOrigin = mother->originVertex()->position().rho();
            parentABSID = abs( mother->particleID().pid() );
          }
        }
      }

      
      if(bAncestor(p)){
        B_child = true;
      }
      if(dAncestor(p)){
        D_child = true;
      }

  }
  
  
  ntuple->column("evt", evt);
  ntuple->column("n_in_evt", n_in_evt);
  ntuple->column("ghost", ghost);
  ntuple->column("pdgcode", pdgcode);
  ntuple->column("true_p", true_p);
  ntuple->column("true_px", true_px);
  ntuple->column("true_py", true_py);
  ntuple->column("true_pz", true_pz);
  ntuple->column("true_pt", true_pt);
  ntuple->column("true_theta", true_theta);
  ntuple->column("true_eta", true_eta);
  ntuple->column("true_phi", true_phi);
  ntuple->column("true_x", true_x);
  ntuple->column("true_y", true_y);
  ntuple->column("true_z", true_z);
  ntuple->column("true_tx", true_tx);
  ntuple->column("true_ty", true_ty);
  ntuple->column("true_nUT_hits", nUT_hits);
  ntuple->column("true_nUT_layers", nUT_layers);
  ntuple->column("true_qop", true_qop);

      ntuple->column("parentABSID", parentABSID);
      ntuple->column("fullInfo", has_all_info);
      ntuple->column("reconstructible_asLong", reconstructible_asLong);
      ntuple->column("isLong", isLong);
      ntuple->column("isDown", isDown);
      ntuple->column("over5", over5);
      ntuple->column("isInVelo", isInVelo);
      ntuple->column("B_child", B_child);
      ntuple->column("D_child", D_child);
      /*
      ntuple->column("true_endvelo_x",true_endvelo_x);
      ntuple->column("true_endvelo_y",true_endvelo_y);
      ntuple->column("true_endvelo_tx",true_endvelo_tx);
      ntuple->column("true_endvelo_ty",true_endvelo_ty);


      ntuple->column("true_midUT_x",true_midUT_x);
      ntuple->column("true_midUT_y",true_midUT_y);
      ntuple->column("true_midUT_tx",true_midUT_tx);
      ntuple->column("true_midUT_ty",true_midUT_ty);


      ntuple->column("true_ftstation_x",true_ftstation_x);
      ntuple->column("true_ftstation_y",true_ftstation_y);
      ntuple->column("true_ftstation_tx",true_ftstation_tx);
      ntuple->column("true_ftstation_ty",true_ftstation_ty);
      */
  
}

//=============================================================================
// Get MC TT hist for particle p
//=============================================================================
int PrDownstreamChecker::getUT_layers(const LHCb::MCParticle* p) const
{
  LinkedFrom<LHCb::STCluster,LHCb::MCParticle> ttLink(evtSvc(),msgSvc(), LHCb::STClusterLocation::UTClusters);
  if (ttLink.notFound()) {
    return 0;
  }
  int nhits = 0;
  std::vector<int> fired;
  const LHCb::STCluster* TTCluster = ttLink.first(p);
  for ( ; 0 != TTCluster; TTCluster = ttLink.next()) {
    if ( !TTCluster->isUT() ) continue;
    nhits++;
    int layer = TTCluster->layer() + 10*TTCluster->station();
    if ( fired.empty() || std::find(fired.begin(), fired.end(), layer) == fired.end() )
      fired.push_back(layer);
  }
  int nlayers = fired.size();
  
  return nlayers;
}

//=============================================================================
// Get MC TT hist for particle p
//=============================================================================
int PrDownstreamChecker::getUT_hits(const LHCb::MCParticle* p) const
{
  LinkedFrom<LHCb::STCluster,LHCb::MCParticle> ttLink(evtSvc(),msgSvc(), LHCb::STClusterLocation::UTClusters);
  if (ttLink.notFound()){
    return 0;
  }
  
  int nhits = 0;
  std::vector<int> fired;
  const LHCb::STCluster* TTCluster = ttLink.first(p);
  for ( ; 0 != TTCluster; TTCluster = ttLink.next()) {
    if ( !TTCluster->isUT() ) continue;
    nhits++;
  }
  
  /*unsigned int station() const;
  /// shortcut for layer
  unsigned int layer() const;
  /// short cut for detRegion
  unsigned int detRegion() const;
  /// short cut for sector
  unsigned int sector() const;
  /// short cut for strip
  unsigned int strip() const;
  */
  return nhits;
}

//=============================================================================
// Look if particle mcPart comes from a b-hadron decay
//=============================================================================
bool PrDownstreamChecker::bAncestor(const LHCb::MCParticle* mcPart) const
{
   bool fromB = false;
   const LHCb::MCParticle* mother = mcPart->mother();
   while ( mother !=0 && fromB == false) {
      fromB = mother->particleID().hasBottom();
      mother = mother->mother();
   }
   return fromB;
}

//=============================================================================
// Fill ntuple with track's information
//=============================================================================
void PrDownstreamChecker::fillTrack(const LHCb::Track *track, Tuples::Tuple ntuple)
{
  // reconstructed quantities
  int type(-1), ndf(-1), q(0);
  double p(-999.99), px(-999.99), py(-999.99), pz(-999.99), pt(-999.99);
  double x(-999.99), y(-999.99), z(-999.99);
  double x_endVelo(-999.99), y_endVelo(-999.99), z_endVelo(-999.99);
  double x_atUT(-999.99), y_atUT(-999.99), z_atUT(-999.99);
  double tx_endVelo(-999.99), ty_endVelo(-999.99), tx_atUT(-999.99), ty_atUT(-999.99);
  double x_atFT(-999.99), y_atFT(-999.99), z_atFT(-999.99);
  double tx_atFT(-999.99), ty_atFT(-999.99);
  
  double chi2(-999.99),chi2PerDoF(-999.99);
  double qop(0.00);
  double eta(0.00),phi(0.00);
  double ds_track_p(0.);
  //double p_corrected;
  int n_UThits(-1);
  int n_FThits(-1);
  int nMeasurementsRemoved(-1);
  int nLHCbIDs(-1);
  //fill info from PrDownstream MC Matching. AD 1-7-14
  int getPreSelection(-1),findMatchingHits(-1),fitXProjection(-1),addUVHits(-1);
  int fitAndRemove(-1),acceptCandidate(-1),beforeStore(-1);
  //HERE. ADD STATES

  double x_zpUT(-9999.0),y_zpUT(-9999.0),tx_zpUT(-9999.0),ty_zpUT(-9999.0);
  double x_zpreUT(-9999.0),y_zpreUT(-9999.0),tx_zpreUT(-9999.0),ty_zpreUT(-9999.0);
  
  if (track) {
    type = track->type();
    px = track->momentum().x();
    py = track->momentum().y();
    pz = track->momentum().z();
    p = track->p();
    pt = track->pt();
    eta = track->pseudoRapidity();
    phi = track->phi();
    q = track->charge();
    chi2 = track->chi2();
    chi2PerDoF = track->chi2PerDoF();
    ndf = track->nDoF();
    x = track->firstState().x();
    y = track->firstState().y();
    z = track->firstState().z();

    //    if(type == 5){
      getPreSelection = track->info(1801,-1);
      findMatchingHits = track->info(1802,-1);
      fitXProjection = track->info(1803,-1);
      addUVHits = track->info(1804,-1);
      fitAndRemove = track->info(1805,-1);
      acceptCandidate = track->info(1806,-1);
      beforeStore = track->info(1807,-1);
      ds_track_p = track->info(1808,-1);
      
      //projected state at zpUT
      x_zpUT = track->info(1809,-1);
      y_zpUT = track->info(1810,-1);
      tx_zpUT = track->info(1811,-1);
      ty_zpUT = track->info(1812,-1);
      //pre ut
      x_zpreUT = track->info(1813,-1);
      y_zpreUT = track->info(1814,-1);
      tx_zpreUT = track->info(1815,-1);
      ty_zpreUT = track->info(1816,-1);
      
      //    }else{
      //      debug()<<"not a downsteam track"<<endmsg;
      //    }
      
    debug() << "chi2 = " << track->chi2() << endmsg;
    
    const LHCb::State *state = track->stateAt(LHCb::State::EndVelo);
    if (state) {
      x_endVelo = state->x();
      y_endVelo = state->y();
      z_endVelo = state->z();
      tx_endVelo = state->tx();
      ty_endVelo = state->ty();
      //p_corrected = p;
      //if(pt<500.){
      //p_corrected = 500./sqrt(tx_endVelo*tx_endVelo + ty_endVelo*ty_endVelo);
      //}
      //if(p_corrected<3000.) p_corrected = 3000.;
    }else{
      debug () << "no velo state" << endmsg;
    }
    
    /*const LHCb::State * state_zmid = &(track->closestState(StateParameters::ZMidTT));
    if(state_zmid){
      debug() << "x_atUT = " << state_zmid->x() << endmsg;
    }else{
      debug() << "no zmid state" << endmsg;
      }*/
      
    //const LHCb::State *state_tt = track->stateAt(LHCb::State::AtTT);
    const LHCb::State * state_tt = &(track->closestState(StateParameters::ZMidTT));
    
    if (state_tt) {
      debug () << "good TT state" << endmsg;
      x_atUT = state_tt->x();
      debug () << "x_atUT = " << x_atUT << endmsg;
      y_atUT = state_tt->y();
      z_atUT = state_tt->z();
      tx_atUT = state_tt->tx();
      ty_atUT = state_tt->ty();
      qop = state_tt->qOverP();
    }else{
      debug() << "no TT state" << endmsg;
    }

    //T state
    const LHCb::State * state_ft = &(track->closestState(StateParameters::ZBegT));
    
    if (state_ft) {
      debug () << "good T state" << endmsg;
      x_atFT = state_ft->x();
      debug () << "x_atFT = " << x_atFT << endmsg;
      y_atFT = state_ft->y();
      z_atFT = state_ft->z();
      tx_atFT = state_ft->tx();
      ty_atFT = state_ft->ty();

    }else{
      debug() << "no TT state" << endmsg;
    }
    
    // Loop over the LHCbIDs of the Track
    n_UThits = 0;
    n_FThits = 0;
    for( std::vector<LHCb::LHCbID>::const_iterator iId = track->lhcbIDs().begin();
	 track->lhcbIDs().end() != iId; ++iId ){
      //if( (*iId).isVelo() || (*iId).isVP() || (*iId).isVL() ) {
      if ( (*iId).isTT() || (*iId).isUT() ) {
	n_UThits++;
      } else if ( (*iId).isIT() || (*iId).isOT() || (*iId).isFT() ) {
	n_FThits++;
      }
    }
    
    nMeasurementsRemoved = track->nMeasurementsRemoved();
    nLHCbIDs = track->nLHCbIDs();
    //std::vector<LHCb::LHCbID> IDs = track->LHCbIDContainer();
    //debug() << "IDs.size() = " << IDs.size() << endmsg;
    /*LinkedFrom<LHCb::STCluster,LHCb::Track> UTlink(evtSvc(),msgSvc(), LHCb::STClusterLocation::UTClusters);
    if (!UTlink.notFound()){
      n_UThits = 0;
      std::vector<int> fired;
      const LHCb::STCluster* UTCluster = UTlink.first(track);
      for ( ; 0 != UTCluster; UTCluster = UTlink.next()) {
	if ( !UTCluster->isUT() ) continue;
	debug() << "layer = " << UTCluster->layer() << endmsg;
	debug() << "station = " << UTCluster->station() << endmsg;
	n_UThits++;
      }
      }*/
    /*LinkedFrom<LHCb::FTCluster,LHCb::Track> FTlink(evtSvc(),msgSvc(), LHCb::FTClusterLocation::Default);
    if (!FTlink.notFound()){
      int nhits = 0;
      std::vector<int> fired;
      const LHCb::FTCluster* FTCluster = FTlink.first(track);
      for ( ; 0 != FTCluster; FTCluster = FTlink.next()) {
	if ( !FTCluster->isFT() ) continue;
	n_FThits++;
      }
      }*/
  }

  ntuple->column("evt",evt);
  ntuple->column("n_in_evt",n_in_evt);
  ntuple->column("nMeasurementsRemoved", nMeasurementsRemoved);
  ntuple->column("nLHCbIDs", nLHCbIDs);
  ntuple->column("n_UThits", n_UThits);
  ntuple->column("n_FThits", n_FThits);
  ntuple->column("type", type);
  ntuple->column("p", p);
  //ntuple->column("p_corrected", p_corrected);
  ntuple->column("px", px);
  ntuple->column("py", py);
  ntuple->column("pz", pz);
  ntuple->column("pt", pt);
  ntuple->column("eta",eta);
  ntuple->column("phi",phi);
  ntuple->column("x", x);
  ntuple->column("y", y);
  ntuple->column("z", z);
  ntuple->column("q", q);
  ntuple->column("qop", qop);
  ntuple->column("x_endVelo", x_endVelo);
  ntuple->column("y_endVelo", y_endVelo);
  ntuple->column("z_endVelo", z_endVelo);
  ntuple->column("tx_endVelo", tx_endVelo);
  ntuple->column("ty_endVelo", ty_endVelo);
  ntuple->column("x_atUT", x_atUT);
  ntuple->column("y_atUT", y_atUT);
  ntuple->column("z_atUT", z_atUT);
  ntuple->column("tx_atUT", tx_atUT);
  ntuple->column("ty_atUT", ty_atUT);

  ntuple->column("x_atFT", x_atFT);
  ntuple->column("y_atFT", y_atFT);
  ntuple->column("z_atFT", z_atFT);
  ntuple->column("tx_atFT", tx_atFT);
  ntuple->column("ty_atFT", ty_atFT);

  
  ntuple->column("chi2", chi2);
  ntuple->column("chi2PerDoF", chi2PerDoF);
  ntuple->column("ndf", ndf);
  ntuple->column("nDownstreamTracks",nDownstreamTracks);
  ntuple->column("nSeedTracks",nSeedTracks);
    
      ntuple->column("getPreSelection",getPreSelection);
      ntuple->column("findMatchingHits",findMatchingHits);
      ntuple->column("fitXProjection",fitXProjection);
      ntuple->column("addUVHits",addUVHits);
      ntuple->column("fitAndRemove",fitAndRemove);
      ntuple->column("acceptCandidate",acceptCandidate);
      ntuple->column("beforeStore",beforeStore);
      ntuple->column("ds_track_p",ds_track_p);
      ntuple->column("x_zpUT_propagated",x_zpUT);
      ntuple->column("y_zpUT_propagated",y_zpUT);
      ntuple->column("tx_zpUT_propagated",tx_zpUT);
      ntuple->column("ty_zpUT_propagated",ty_zpUT);

      ntuple->column("x_zpreUT_propagated",x_zpreUT);
      ntuple->column("y_zpreUT_propagated",y_zpreUT);
      ntuple->column("tx_zpreUT_propagated",tx_zpreUT);
      ntuple->column("ty_zpreUT_propagated",ty_zpreUT);
      
      //HERE add columns for good and total before and after fitandremove
      /*
      int ngood_before_x1(-2), ngood_before_u(-2),ngood_before_v(-2),ngood_before_x2(-2);
      int ntot_before_x1(-2), ntot_before_u(-2),ntot_before_v(-2),ntot_before_x2(-2);
      int ngood_after_x1(-2), ngood_after_u(-2),ngood_after_v(-2),ngood_after_x2(-2);
      int ntot_after_x1(-2), ntot_after_u(-2),ntot_after_v(-2),ntot_after_x2(-2);
      */
      //add the nhits stuff
      
}


#endif //PrDownstreamChecker_CPP

