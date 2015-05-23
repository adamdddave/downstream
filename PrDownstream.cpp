// Include files

#include <algorithm>
#include <cmath>
#include <assert.h>
// from boost
#include <boost/assign/list_of.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>

// from ROOT
#include <Math/CholeskyDecomp.h>
using ROOT::Math::CholeskyDecomp;

// from Gaudi
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/SystemOfUnits.h"
#include "GaudiKernel/Point3DTypes.h"

#include "STDet/DeSTDetector.h"
#include "STDet/DeUTSector.h"

#include "Event/Track.h"

#include "PrKernel/IPrDebugUTTool.h"
#include "TfKernel/RecoFuncs.h"
#include "TfKernel/IUTHitCreator.h"
#include "Event/StateParameters.h"

// local
#include "PrDownTrack.h"
#include "PrDownTrack2.h"
#include "PrDownstream.h"

//-----------------------------------------------------------------------------
// Implementation file for class : PrDownstream
//
// 2005-11-18 : Olivier Callot
// 2009-10-07 : Sascha Stahl
// 2015-05-18 : Adam Davis
//-----------------------------------------------------------------------------

DECLARE_ALGORITHM_FACTORY( PrDownstream )

/*
bool sortByChi2Tr( PrDownTrack lhs, PrDownTrack rhs) {return lhs.chisq()< rhs.chisq();}
bool sortByHitsTr(PrDownTrack lhs, PrDownTrack rhs) {return lhs.hits().size()> rhs.hits().size();}*/
bool sortByChi2Tr( PrDownTrack2 lhs, PrDownTrack2 rhs) {return lhs.chisq()< rhs.chisq();}
bool sortByHitsTr(PrDownTrack2 lhs, PrDownTrack2 rhs) {return lhs.hits().size()> rhs.hits().size();}
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrDownstream::PrDownstream( const std::string& name,
                      ISvcLocator* pSvcLocator)
  : GaudiAlgorithm ( name , pSvcLocator ),
    m_downTime( 0 ),
    m_printing( false ),
    m_magnetOff( false ),
    m_magFieldSvc(NULL),
    m_utHitManager(NULL),
    m_debugTool(NULL),
    m_timerTool(NULL)
{
  declareProperty( "InputLocation" , m_inputLocation  = ""    );
  declareProperty( "OutputLocation", m_outputLocation = LHCb::TrackLocation::Downstream );
  declareProperty( "ForwardLocation", m_forwardLocation = LHCb::TrackLocation::Forward );
  declareProperty( "MatchLocation", m_matchLocation = LHCb::TrackLocation::Match );
  declareProperty( "PrintTracks"   , m_printTracks    = false );
  declareProperty( "TimingMeasurement", m_doTiming = false);
  declareProperty( "deltaP"        , m_deltaP        = 2.0  );
  declareProperty( "xPredTol"      , m_xPredTol      = 150. * Gaudi::Units::mm * Gaudi::Units::GeV  );
  declareProperty( "xPredTol2"      , m_xPredTol2      = 20. * Gaudi::Units::mm );
  declareProperty( "TolMatch"      , m_tolX          = 1.5 *  Gaudi::Units::mm );
  declareProperty( "TolUV"         , m_tolUV         = 2.0 *  Gaudi::Units::mm );
  declareProperty( "TolMomentum"   , m_tolMomentum   = 20000.    );
  declareProperty( "maxWindowSize" , m_maxWindow     = 10.0 * Gaudi::Units::mm );
  declareProperty( "MaxDistance"   , m_maxDistance   = 0.30 * Gaudi::Units::mm );
  declareProperty( "MaxChisq"      , m_maxChisq      = 20.        );
  // Delta Chisq wrt the best track if multiple tracks are considered
  declareProperty( "DeltaChisq"      , m_deltaChisq      = 3.5        );
  declareProperty( "errorZMagnet"  , m_errZMag       = 30. *  Gaudi::Units::mm );
  declareProperty( "minUTx"        , m_minUTx        = 35. *  Gaudi::Units::mm );
  declareProperty( "minUTy"        , m_minUTy        = 35. *  Gaudi::Units::mm );
  // Define parameters for MC09 field, zState = 9410
  {
    std::vector<double> tmp = boost::assign::list_of(5376.8)(-3895.12)(309.877)(85527.9);
    declareProperty( "zMagnetParams" , m_zMagnetParams  = tmp);
  }
  {
    std::vector<double> tmp = boost::assign::list_of(1148.65)(961.786)(5326.81);
    declareProperty( "momentumParams", m_momentumParams = tmp);
  }
  {
    std::vector<double> tmp =  boost::assign::list_of(5.)(2000.);
    declareProperty( "yParams"       , m_yParams        = tmp);
  }

  declareProperty( "zUT"           , m_zUT           = 2485.* Gaudi::Units::mm );
  declareProperty( "zUTa"          , m_zUTa          = 2350.* Gaudi::Units::mm );

  declareProperty( "StateErrorX2"  , m_stateErrorX2  =   4.0);
  declareProperty( "StateErrorY2"  , m_stateErrorY2  = 400.);
  declareProperty( "StateErrorTX2" , m_stateErrorTX2 =  6.e-5);
  declareProperty( "StateErrorTY2" , m_stateErrorTY2 =  1.e-4);
  declareProperty( "StateErrorP"   , m_stateErrorP   =  0.15);
  declareProperty( "MinPt"                 , m_minPt                 =  0. * Gaudi::Units::MeV );
  declareProperty( "MinMomentum"           , m_minMomentum           =  0. * Gaudi::Units::GeV );

  // Change this in order to remove hits and T-tracks used for longtracks.
  // RemoveAll configures that everything is removed.
  // If false only hits and T-tracks from good longtracks are removed.
  // The criterion for this is the Chi2 of the longtracks from the fit.
  declareProperty( "RemoveUsed"    , m_removeUsed    = false     );
  declareProperty( "RemoveAll"    ,  m_removeAll = false      );
  declareProperty( "LongChi2"      , m_longChi2     = 1.5       );

  // Modes the algorithm can run in
  // true == original Downstream tracking 
  // false == New approach to downstream tracking
  declareProperty( "Mode"           , m_legacy     = true       );
  // HitRemoval in fitAndRemove
  declareProperty( "RemoveHits"    ,  m_removeHits = true      );
  // Number of tracks per seed to keep if mode == true
  declareProperty( "TracksToKeep"   , m_tracksToKeep     = 3       );

  //== debugging options
  declareProperty( "SeedKey"       , m_seedKey       = -1    );
  declareProperty( "WithDebugTool" , m_withDebugTool = false );
  declareProperty( "DebugTool"     , m_debugToolName = "PrDebugUTTruthTool" );
  declareProperty("NumberIterationsNewChi2",m_number_Iterations_new_chi2=0);
  declareProperty("PrintVerbose", m_printVerbose = false);
  
  //declareProperty( "PrintExtraInformation",m_print_extra_info = false);
  
}
//=============================================================================
// Destructor
//=============================================================================
PrDownstream::~PrDownstream() {}

//=============================================================================
// Initialization
//=============================================================================
StatusCode PrDownstream::initialize() {
  StatusCode sc = GaudiAlgorithm::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Initialize" << endmsg;

  m_utHitManager = tool<Tf::UTStationHitManager <PrUTHit> >("PrUTStationHitManager");

  m_debugTool = 0;
  if ( m_withDebugTool ) {
    m_debugTool = tool<IPrDebugUTTool>( m_debugToolName );
  }
  

  info() << "========================================"  << endmsg
         << " deltaP             = " << m_deltaP        << endmsg
         << " xPredTol           = " << m_xPredTol      << endmsg
         << " xPredTol2           = " << m_xPredTol2      << endmsg
         << " TolMomentum        = " << m_tolMomentum   << endmsg
         << " TolMatch           = " << m_tolX          << endmsg
         << " TolUV              = " << m_tolUV         << endmsg
         << " MaxDistance        = " << m_maxDistance   << endmsg
         << " MaxChisq           = " << m_maxChisq      << endmsg
         << " MaxWindowSize      = " << m_maxWindow     << endmsg
         << " RemoveUsed         = " << m_removeUsed    << endmsg
         << " RemoveAll          = " << m_removeAll     << endmsg
         << " LongChisq          = " << m_longChi2      << endmsg
         << "========================================"  << endmsg;

  info() << "zMagnetParams ";
  for ( unsigned int kk = 0; m_zMagnetParams.size() > kk ; kk++) {
    info() << m_zMagnetParams[kk] << " ";
  }
  info() << endmsg << "momentumParams ";
  for ( unsigned int kk = 0; m_momentumParams.size() > kk ; kk++) {
    info() << m_momentumParams[kk] << " ";
  }
  info() << endmsg ;
  if ( 3 > m_zMagnetParams.size() ) {
    return Warning( "Not enough zMagnetParams" );
  }
  if ( 3 > m_momentumParams.size() ) {
    return Warning( "Not enough momentumParams" );
  }

  m_magFieldSvc = svc<ILHCbMagnetSvc>( "MagneticFieldSvc", true );

  if ( m_doTiming) {
    m_timerTool = tool<ISequencerTimerTool>( "SequencerTimerTool" );
    m_timerTool->increaseIndent();
    m_downTime = m_timerTool->addTimer( "Internal PrDownstream" );
    m_timerTool->decreaseIndent();
  }

  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode PrDownstream::execute() {

  if ( m_doTiming ) m_timerTool->start( m_downTime );

  m_printing = msgLevel( MSG::DEBUG );
  if (m_printing) debug() << "==> Execute" << endmsg;

  //== if no specified input location, get seed if it exists, else Tsa.
  if ( m_inputLocation.empty() ) {
    if ( exist<LHCb::Tracks>( LHCb::TrackLocation::Seed ) ) {
      m_inputLocation = LHCb::TrackLocation::Seed;
    } else {
      m_inputLocation = LHCb::TrackLocation::Tsa;
    }
  }
  //==========================================================================
  // Prepare hits in UT, optional: remove used in PrForward.
  //==========================================================================

  utCoordCleanup();

  //==========================================================================
  // Prepare the tracks, optional: removing the seends already used in Match.
  //==========================================================================
  LHCb::Tracks* inTracks   = get<LHCb::Tracks>( m_inputLocation );

  //== Local container of good track.
  std::vector<LHCb::Track*> myInTracks;
  myInTracks.reserve(inTracks->size());

  std::vector< PrDownTrack2 > finalUVtracks;
  finalUVtracks.reserve( inTracks->size() );

  //=== Optional: Remove the seeds already used in Match, check by ancestor of Match tracks
  //=== Prepare T-Seeds
  if(m_legacy) prepareSeeds( inTracks, myInTracks );
  else prepareSeedsNew( inTracks, myInTracks );


  if ( m_printing ) debug() << "-- Started from " << inTracks->size() << " Seed, left with "
      << myInTracks.size() << " tracks after Matched removal." << endmsg;

  //==========================================================================
  // Get the output container ( created by PrDataStore )
  //==========================================================================
  LHCb::Tracks* prefinalTracks = new LHCb::Tracks();
  prefinalTracks->reserve(100);
  LHCb::Tracks* finalTracks = new LHCb::Tracks();
  finalTracks->reserve(100);
  put( finalTracks, m_outputLocation);
  //add for updated TSeeds for debugging
  LHCb::Tracks* updated_TSeeds = new LHCb::Tracks();
  if(m_withDebugTool){
    updated_TSeeds->reserve(myInTracks.size());
    put( updated_TSeeds, m_outputLocation+"UpdatedSeed" );
  }
  
  //==========================================================================
  // Main loop on tracks
  //==========================================================================
  Tf::UTStationHitManager<PrUTHit>::HitRange utCoords = m_utHitManager->hits();

  for ( LHCb::Tracks::const_iterator itT = myInTracks.begin();
          myInTracks.end() != itT; ++itT ) {
      LHCb::Track* tr = (*itT);
      if(m_withDebugTool){//AD add debug information for plots, etc.
        if(m_printing){ 
          if(m_printVerbose)m_debugTool->printTSeedInfo(tr);
          else{m_debugTool->printTSeedInfoShort(tr);}
        }
        
        m_debugTool->recordStepInProcess("reconstructible_asDown",
                                         m_debugTool->isTrackReconstructible(*tr));  
      }
      
      if ( 0 <= m_seedKey && m_seedKey == tr->key() ) m_printing = true;
      
      const double magScaleFactor = m_magFieldSvc->signedRelativeCurrent() ;

      if( std::abs(magScaleFactor) > 1e-6 ){
          m_magnetOff = false;
      } else m_magnetOff = true;
      //PrDownTrack track( tr, m_zUT, m_zMagnetParams, m_momentumParams, m_yParams, m_errZMag, magScaleFactor*(-1) );
      PrDownTrack2 track(tr,m_zMagnetParams,m_momentumParams,magScaleFactor*(-1));
      

      if(m_legacy) {
        //Y. Xie: get rid of particles from beampipe
        /*
        const double xAtUT = track.xAtZ( m_zUTa );
        const double yAtUT = track.yAtZ( m_zUTa );
        if( m_minUTx > fabs(xAtUT) &&
        m_minUTy > fabs(yAtUT)       ) continue;*/
      }

      const double deltaP = track.moment() * track.state()->qOverP() - 1.;
      
      if ( m_printing ) {
          for ( auto itH : utCoords ){
              PrUTHit* hit = itH;
              if (hit->hit()->ignore()) continue;
              const double yTrack = track.yAtZ( 0. );
              const double tyTr   = track.slopeY();
              updateUTHitForTrack( hit, yTrack, tyTr );
              //if(m_withDebugTool)m_debugTool->PrintHitTable(hit);//this is the wrong place to do this. AD 5-22-15
          }

          info() << "Track " << tr->key()
            /*
                 << format( " [%7.2f %7.2f GeV] x,y(UTa)%7.1f%7.1f dp/p%7.3f errXMag%6.1f YUT%6.1f",
                            .001/track.state()->qOverP(), .001*track.moment(),
                            track.xAtZ( m_zUTa ), track.yAtZ( m_zUTa ), deltaP,
                            track.errXMag(), track.errYMag() )*/
                 << endmsg;
            info() << format( " Y slope %8.6f computed %8.6f", track.state()->ty(), track.slopeY() )
              << endmsg;
          info()<<"calling debugutclusterontrack"<<endmsg;
          if ( m_debugTool ) m_debugTool->debugUTClusterOnTrack( tr, utCoords.begin(), utCoords.end() );
          info()<<"done calling debugutclusterontrack"<<endmsg;
      }
      

      // check for compatible momentum
      if ( m_deltaP < fabs(deltaP) ) {
          if ( m_printing ) info() << "   --- deltaP " << deltaP << " -- rejected" << endmsg;
          if ( !m_magnetOff ) continue;
          continue;
      }

      // Get hits in UT around a first track estimate
      getPreSelection( track );
      
      
      if(m_withDebugTool){
        std::vector<LHCb::LHCbID> x1lhcbids,ulhcbids,vlhcbids,x2lhcbids;
        for(auto xhit:m_xHits){
          if(xhit->planeCode()==0)x1lhcbids.push_back(xhit->hit()->lhcbID());
          if(xhit->planeCode()==3)x2lhcbids.push_back(xhit->hit()->lhcbID());
        }
        for(auto uvhit:m_uvHits){
          if(uvhit->planeCode()==1)ulhcbids.push_back(uvhit->hit()->lhcbID());
          if(uvhit->planeCode()==2)vlhcbids.push_back(uvhit->hit()->lhcbID());
        }
        m_debugTool->recordStepInProcess("getPreselection",
                                         m_debugTool->isPreselGood(x1lhcbids,ulhcbids,vlhcbids,x2lhcbids,tr));
        x1lhcbids.clear();ulhcbids.clear();vlhcbids.clear();x2lhcbids.clear();
      }
      if(m_withDebugTool){
        info()<<"calling verbose printing"<<endmsg;
        
        for(auto hit: m_xHits){
          if(m_printVerbose) m_debugTool->PrintHitTable(hit);
          else{m_debugTool->PrintHitTableShort(hit);}
        }//this is the right place to do this. AD 5-22-15
        for(auto hit: m_uvHits){
          if(m_printVerbose) m_debugTool->PrintHitTable(hit);
          else{m_debugTool->PrintHitTableShort(hit);}
        }
        
      }
      
      //now the preselection is there, we can print everything
      
      int nbOK = 0;
      //PrDownTrack bestTrack( track );
      PrDownTrack2 bestTrack( track );

      //std::vector<PrDownTrack> trackCands;
      std::vector<PrDownTrack2> trackCands;
      trackCands.reserve(300);

      PrUTHits MatchingXHits;
      MatchingXHits.reserve(16);
      int maxPoints = 0;
      double minChisq = m_maxChisq + 10.;

      if(m_legacy) {
          PrUTHits::const_iterator itH;
          if ( 2 < ( m_xHits.size() +  m_uvHits.size() ) ) {
              //==============================================================
              // Try to find a candidate: X first, then UV.
              //==============================================================

              for ( itH = m_xHits.begin(); m_xHits.end() != itH; itH++ ) {
                  PrUTHit* myHit = *itH;
                  if ( myHit->hit()->testStatus( Tf::HitBase::UsedByPatDownstream ) ) continue;
                  double meanZ = myHit->z();
                  double posX  = myHit->x( );
                  int myPlane  = myHit->planeCode();

                  track.startNewCandidate();

                  // Create track estimate with one x hit
                  //double slopeX = (track.xMagnet() - posX + track.sagitta( meanZ)) / (track.zMagnet() - meanZ);
                  double slopeX = (track.xMagnet()- posX)/(track.zMagnet() - meanZ);
                  track.setSlopeX( slopeX );
                  
                  if ( m_printing ) {
                      const double tolMatch = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolX)) ?
                          m_maxWindow : (m_tolX + m_tolMomentum / track.state()->p());
                      info() << endmsg 
                          << format( "... start plane %1d x%8.2f z%8.1f slope%8.2f tolMatch%7.3f",
                                  myPlane, posX, meanZ, 1000. * slopeX, tolMatch )
                          << endmsg;
                  }        

                  // Fit x projection
                  findMatchingHits( MatchingXHits, track, myPlane );
                  if(m_withDebugTool){
                    std::vector<int> xids,uvids;
                    for(auto hit: MatchingXHits){xids.push_back(hit->hit()->lhcbID().lhcbID());}
                    for(auto hit: m_uvHits){uvids.push_back(hit->hit()->lhcbID().lhcbID());}
                    m_debugTool->recordStepInProcess("findMatchingHits",
                                                     m_debugTool->CheckMCinIntermediateHits(xids,uvids));
                  }
                  
                  //fitXProjection( MatchingXHits, track, myHit );
                  fitXProjection2( MatchingXHits, track, myHit );
                  if(m_withDebugTool){
                    std::vector<int> xids,uvids;
                    for(auto hit: track.hits()){xids.push_back(hit->hit()->lhcbID().lhcbID());}
                    for(auto hit: m_uvHits){uvids.push_back(hit->hit()->lhcbID().lhcbID());}
                    m_debugTool->recordStepInProcess("fitXProjection",
                                                     m_debugTool->CheckMCinIntermediateHits(xids,uvids)); 
                  }
                  // Fit 3D track
                  //addUVHits( track );
                  addUVHits2( track );
                  if(m_withDebugTool){
                    std::vector<int> ids;
                    for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                    m_debugTool->recordStepInProcess("addUVHits",
                                                     m_debugTool->CheckMCinTrackHits(ids)); 
                  }
                  
                  if ( 3 > track.hits().size() ) {
                      if ( m_printing ) info() << " === not enough hits" << endmsg;
                      continue;
                  }
                  fitAndRemove2( track );
                  if(m_withDebugTool){
                    std::vector<int> ids;
                    for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                    m_debugTool->recordStepInProcess("fitAndRemove",
                                                     m_debugTool->CheckMCinTrackHits(ids)); 
                  }
                  //fitAndRemove( track );
                  //here now for update

                  // Check if candidate is better than the old one
                  if ( !acceptCandidate( track, bestTrack, maxPoints, minChisq ) ) continue;
                  if(m_withDebugTool){
                    std::vector<int> ids;
                    for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                    m_debugTool->recordStepInProcess("acceptCandidate",
                                                     m_debugTool->CheckMCinTrackHits(ids)); 
                  }
                  // New best track
                  bestTrack = track;
                  nbOK++;
              }
          }

          //== Debug the track.
          if ( 0 == nbOK ) {
              if ( m_printing ) {
                  info() << format( "Track %3d P=%7.2f GeV --- discarded, not enough planes",
                          tr->key(),  .001*track.moment() ) << endmsg;
              }
          } else {

              // Best one  
              track = bestTrack;

              // add hits in x layers in overlap regions of UT
              addOverlapRegions( track );     
              if(m_withDebugTool){
                std::vector<int> ids;
                for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                m_debugTool->recordStepInProcess("beforeStore",
                                                 m_debugTool->CheckMCinTrackHits(ids)); 
              }
              //=== Store the track
              storeTrack( track, finalTracks, tr,updated_TSeeds );

          }

      } else {
        if ( 2 < ( m_xHits.size() +  m_uvHits.size() ) ) {
          //==============================================================
          // Try to find a candidate: X first, then UV.
          //==============================================================

          for ( auto itH : m_xHits ) {
            if ( m_doTiming ) m_timerTool->start( m_buildTime );
            PrUTHit* myHit = itH;
            if ( myHit->hit()->testStatus( Tf::HitBase::UsedByPatDownstream ) ) continue;
            double meanZ = myHit->z();
            double posX  = myHit->x( );
            int myPlane  = myHit->planeCode();

            if ( 1!=myHit->mask() ) continue;

            track.startNewCandidate();

            //Beampipe cut
            //TODO: Hardcoded for now - we may want to have this as a tuning parameter
            if(std::sqrt(track.xMagnet()*track.xMagnet()+track.yMagnet()*track.yMagnet()) < 80) continue;
            // Create track estimate with one x hit
            //double slopeX = (track.xMagnet() - posX + track.sagitta( meanZ)) / (track.zMagnet() - meanZ);
            double slopeX = (track.xMagnet() - posX) / (track.zMagnet() - meanZ);
            track.setSlopeX( slopeX );
            
            if ( m_printing ) {
              const double tolMatch = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolX)) ?
                m_maxWindow : (m_tolX + m_tolMomentum / track.state()->p());
              info() << endmsg
                     << format( "... start plane %1d x%8.2f z%8.1f slope%8.2f tolMatch%7.3f",
                                myPlane, posX, meanZ, 1000. * slopeX, tolMatch )
                     << endmsg;
            }

            // Fit x projection
            findMatchingHits( MatchingXHits, track, myPlane );
            if(m_withDebugTool){
              std::vector<int> xids,uvids;
              for(auto hit: MatchingXHits){xids.push_back(hit->hit()->lhcbID().lhcbID());}
              for(auto hit: m_uvHits){uvids.push_back(hit->hit()->lhcbID().lhcbID());}
              m_debugTool->recordStepInProcess("findMatchingHits",
                                               m_debugTool->CheckMCinIntermediateHits(xids,uvids));
            }
                  
            //fitXProjection( MatchingXHits, track, myHit );
            fitXProjection2( MatchingXHits, track, myHit );
            if(m_withDebugTool){
              std::vector<int> xids,uvids;
              for(auto hit: track.hits()){xids.push_back(hit->hit()->lhcbID().lhcbID());}
              for(auto hit: m_uvHits){uvids.push_back(hit->hit()->lhcbID().lhcbID());}
              m_debugTool->recordStepInProcess("fitXProjection",
                                               m_debugTool->CheckMCinIntermediateHits(xids,uvids)); 
            }

            std::vector< PrDownTrack2 > UVtracks;

            //addUVHits( track, &UVtracks );
            addUVHits2( track, &UVtracks );
            if(m_withDebugTool){
              std::vector<int> ids;
              for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
              m_debugTool->recordStepInProcess("addUVHits",
                                               m_debugTool->CheckMCinTrackHits(ids)); 
            }
                  
            if(0==UVtracks.size()) continue;

            for(auto trc : UVtracks) {
                    
              fitAndRemove2( trc );
              //fitAndRemove( trc );
              if(m_withDebugTool){
                std::vector<int> ids;
                for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                m_debugTool->recordStepInProcess("fitAndRemove",
                                                 m_debugTool->CheckMCinTrackHits(ids)); 
              }
              //Beampipe cut - once again after the fit
              //TODO: Hardcoded for now - we may want to have this as a tuning parameter
              //TODO: We might want to reconsider if this is necessary twice
              if(std::sqrt(track.xMagnet()*track.xMagnet()+track.yMagnet()*track.yMagnet()) < 80) continue;

              if(legitCandidate( trc ) ){
                trackCands.push_back( trc );
                if(m_withDebugTool){
                  std::vector<int> ids;
                  for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
                  m_debugTool->recordStepInProcess("acceptCandidate",
                                                   m_debugTool->CheckMCinTrackHits(ids)); 
                }//debug tool
              }//legit candidate
            }//uvtracks
          }//xhits
        }
        
        //== Debug the track.
        if ( 0==trackCands.size() ) {
          if ( m_printing ) {
            info() << format( "Track %3d P=%7.2f GeV --- discarded, not enough planes",
                              tr->key(),  .001*track.moment() ) << endmsg;
          }
        } else {

          std::stable_sort(trackCands.begin(), trackCands.end(), sortByChi2Tr);
          std::stable_sort(trackCands.begin(), trackCands.end(), sortByHitsTr);
          std::vector<PrDownTrack2> trackCandsClean;
          trackCandsClean.reserve(trackCands.size());


          for(unsigned count = 0; count < trackCands.size(); count++) {

            //TODO: Hardcoded for now - we may want to have this as a tuning parameter
            //TODO: We might want to reconsider the meaning of this cut
            if(std::sqrt(track.slopeX()*track.slopeX()+track.slopeY()*track.slopeY()) < 0.010) continue;
            track = trackCands[count];

            bool clone = false;
            for ( std::vector<PrDownTrack2>::const_iterator itT1 = trackCandsClean.begin(); trackCandsClean.end() != itT1; ++itT1 ) {
              PrDownTrack2 testTrack = *itT1;

              PrUTHits::const_iterator itHa;
              PrUTHits::const_iterator itHb;

              int hitsmin = std::min(testTrack.hits().size(), track.hits().size());
              int hitscommon = 0;
              for ( itHa = track.hits().begin(); track.hits().end() != itHa; itHa++ ) {
                for ( itHb = testTrack.hits().begin(); testTrack.hits().end() != itHb; itHb++ ) {
                  if((*itHa)->hit()->lhcbID()==(*itHb)->hit()->lhcbID())
                    hitscommon++;
                }
              }
                      
              double common = (double)hitscommon/(double)hitsmin;
              
              //TODO: Hardcoded for now - we may want to have this as a tuning parameter
              if(common>0.66) {
                clone = true;
                break;
              }
            }
            if(clone) break;

            trackCandsClean.push_back(track);
          }
          
          std::stable_sort(trackCandsClean.begin(), trackCandsClean.end(), sortByChi2Tr);
          unsigned maxcount = 0;
          if(trackCandsClean.size()> m_tracksToKeep)
            maxcount = m_tracksToKeep;
          else
            maxcount = trackCandsClean.size();

          if(m_withDebugTool){
            std::vector<int> ids;
            for(auto hit:track.hits()){ids.push_back(hit->hit()->lhcbID().lhcbID());}
            m_debugTool->recordStepInProcess("beforeStore",
                                             m_debugTool->CheckMCinTrackHits(ids)); 
          }
          for(unsigned count = 0; count < maxcount; count++) {
            track = trackCandsClean[count];

            if(count == 0) {
              track.sortFinalHits();
              
              //=== Store the track
              storeTrack( track, finalTracks, tr,updated_TSeeds );
            }

            if(count > 0) {
              //std::cout << trackCandsClean[count].chisq()-trackCandsClean[0].chisq() << std::endl;
              if(trackCandsClean[count].chisq()-trackCandsClean[0].chisq()<m_deltaChisq) {
                // add hits in x layers in overlap regions of UT
                //////addOverlapRegions( track );
                //////if ( m_maxChisq < track.chisq() ) continue;
                //////if ( m_maxDisplY < std::abs(track.displY()) ) continue;



                track.sortFinalHits();

                //=== Store the track
                storeTrack( track, finalTracks, tr, updated_TSeeds );
              }
            }
          }
        }
      }
  }
  


  ///might have to put that back in when we not to decide to veto long tracks
  if( UNLIKELY( msgLevel(MSG::DEBUG) ) )
    debug() << "Found " << finalTracks->size() << " tracks." << endmsg;

  if ( m_doTiming ) {
      m_timerTool->stop( m_downTime );
  }


  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode PrDownstream::finalize() {

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Finalize" << endmsg;

  return GaudiAlgorithm::finalize();  // must be called after all other actions
}

//=========================================================================
//  Cleanup already used T-Seeds - New implementation
//=========================================================================
void PrDownstream::prepareSeedsNew(LHCb::Tracks* inTracks, std::vector<LHCb::Track*>& myInTracks){
    LHCb::Tracks* match = getIfExists<LHCb::Tracks>( m_matchLocation );
    LHCb::Tracks* forward = getIfExists<LHCb::Tracks>( m_forwardLocation );
    if ( NULL != match ) {
        if (!m_removeUsed) {
            myInTracks.insert( myInTracks.end(), inTracks->begin(), inTracks->end() );
        } else {
            if (m_printing) debug()<<"Remove seeds and ut hits from Match tracks"<<endmsg;
            for ( LHCb::Tracks::const_iterator itT = inTracks->begin();
                    inTracks->end() != itT; itT++ ) {
                LHCb::Track* tr = (*itT);
                bool forwd = true;
                
                for( auto forwardTr : *forward ) {
                    int ft1 = 0;
                    int ft2 = 0;
                    int ftc = 0;
                    bool hasUT = false;

                    for (std::vector< LHCb::LHCbID>::const_iterator id1 = forwardTr->lhcbIDs().begin();
                            forwardTr->lhcbIDs().end() != id1; ++id1 ) {
                        if(!(*id1).isFT()) continue;
                        for (std::vector< LHCb::LHCbID>::const_iterator id2 = tr->lhcbIDs().begin();
                                tr->lhcbIDs().end() != id2; ++id2 ) {
                            if(!(*id2).isFT()) continue;
                            if((*id1)==(*id2)) ftc++;
                        }
                    }
                    for (std::vector< LHCb::LHCbID>::const_iterator id1 = forwardTr->lhcbIDs().begin();
                            forwardTr->lhcbIDs().end() != id1; ++id1 ) {
                        if((*id1).isFT()) ft1++;
                        if((*id1).isUT()) hasUT=true;
                    }
                    for (std::vector< LHCb::LHCbID>::const_iterator id2 = tr->lhcbIDs().begin();
                            tr->lhcbIDs().end() != id2; ++id2 ) {
                        if((*id2).isFT()) ft2++;
                    }
                    
                    if(ft1!=0 && ft2!=0 && (double)ftc/(double)std::min(ft1, ft2)>=1 && forwardTr->chi2PerDoF() < m_longChi2 && hasUT) {
                        forwd = false;
                        tagUsedUT( forwardTr );
                    }
                    if ( !forwd ) break;
                }

                //TODO: Keep this in case we want to use both information from the Match and the Forward in the offline tracking
                /*
                bool mat =true;
                if ( m_printing ) debug() << "Seed " << tr->key();
                for( auto matchTr : *match ) {
                    const SmartRefVector<LHCb::Track>& ancestors = matchTr->ancestors();
                    for ( SmartRefVector<LHCb::Track>::const_iterator itA = ancestors.begin();
                            ancestors.end() != itA; ++itA ) {
                        const LHCb::Track* pt = (*itA);
                        if ( tr == pt ) {
                            if ( m_printing ) debug() << " is used in match " << matchTr->key();
                            if ( m_removeAll || (matchTr->chi2PerDoF() < m_longChi2 && !forwd) ) {
                                if ( m_printing ) debug() << " good longtrack " << matchTr->key()<<endmsg;
                                mat = false;
                                tagUsedUT( matchTr );
                                break;
                            }
                            break;
                        }
                    }
                    if ( !mat ) break;
                }

                if (mat) {
                    myInTracks.push_back( tr );
                    if ( m_printing ) debug() << " will be processed";
                }
                if ( m_printing ) debug() << endmsg;
                */
                if (forwd) {
                    myInTracks.push_back( tr );
                    if ( m_printing ) debug() << " will be processed";
                }
                if ( m_printing ) debug() << endmsg;

            }
            }
            } else {  //== Copy tracks without ancestor...
                for(auto tr : *inTracks) {
                    //== Ignore tracks with ancestor = forward...
                    if ( m_removeUsed && 0 < tr->ancestors().size() ) continue;
                    myInTracks.push_back( tr );
                }
            }
            }

//=========================================================================
//  Cleanup already used T-Seeds
//=========================================================================
void PrDownstream::prepareSeeds(LHCb::Tracks* inTracks, std::vector<LHCb::Track*>& myInTracks){
  LHCb::Tracks* match = getIfExists<LHCb::Tracks>( m_matchLocation );
  if ( NULL != match ) {
    if (!m_removeUsed) {
      myInTracks.insert( myInTracks.end(), inTracks->begin(), inTracks->end() );
    } else {
      if (m_printing) debug()<<"Remove seeds and ut hits from Match tracks"<<endmsg;
      for ( LHCb::Tracks::const_iterator itT = inTracks->begin();
            inTracks->end() != itT; itT++ ) {
        LHCb::Track* tr = (*itT);
        bool store =true;
        if ( m_printing ) debug() << "Seed " << tr->key();
        for( auto matchTr : *match ) {
          const SmartRefVector<LHCb::Track>& ancestors = matchTr->ancestors();
          for ( SmartRefVector<LHCb::Track>::const_iterator itA = ancestors.begin();
                ancestors.end() != itA; ++itA ) {
            const LHCb::Track* pt = (*itA);
            if ( tr == pt ) {
              if ( m_printing ) debug() << " is used in match " << matchTr->key();
              if ( m_removeAll || (matchTr->chi2PerDoF() < m_longChi2 ) ) {
                if ( m_printing ) debug() << " good longtrack " << matchTr->key()<<endmsg;
                store = false;
                tagUsedUT( matchTr );
                break;
              }
              break;
            }
          }
          if ( !store ) break;
        }
        if (store) {
          myInTracks.push_back( tr );
          if ( m_printing ) debug() << " will be processed";
        }
        if ( m_printing ) debug() << endmsg;
      }
    }
  } else {  //== Copy tracks without ancestor...
    for(auto tr : *inTracks) {
      //== Ignore tracks with ancestor = forward...
      if ( m_removeUsed && 0 < tr->ancestors().size() ) continue;
      myInTracks.push_back( tr );
    }
  }
}



//=========================================================================
//  Cleanup already used UT coord
//=========================================================================
void PrDownstream::utCoordCleanup ( ) {
  Tf::UTStationHitManager<PrUTHit>::HitRange utCoords = m_utHitManager->hits();
  for(auto hit : utCoords)
    hit->hit()->setStatus( Tf::HitBase::UsedByPatDownstream, false );

  //== Tag hit used in forward
  if ( m_removeUsed ) {
    if ( exist<LHCb::Tracks>( m_forwardLocation ) ) {
      if (m_printing) debug()<<"Remove UT hits from Forward tracks from location "
                             <<m_forwardLocation <<endmsg;
      LHCb::Tracks* tracks = get<LHCb::Tracks>( m_forwardLocation );
      for(auto tr : *tracks) {
        if (m_removeAll || (tr->chi2PerDoF()<m_longChi2)) tagUsedUT( tr );
      }
    }
  }
}

//=========================================================================
//  Tag used UT hits on this track.
//=========================================================================
void PrDownstream::tagUsedUT( const LHCb::Track* tr ) {
  Tf::UTStationHitManager<PrUTHit>::HitRange utCoords = m_utHitManager->hits();
  for(auto id : tr->lhcbIDs()) {
    if ( !id.isUT() ) continue;
    for(auto hit : utCoords) {
      if ( hit->hit()->lhcbID() == id ) {
        if (m_printing) debug()<<"tag hit as used "<<hit->hit()->lhcbID()<<endmsg;
        hit->hit()->setStatus( Tf::HitBase::UsedByPatMatch, true );
        break;
      }
    }
  }
}

//=========================================================================
//  Get the PreSelection
//=========================================================================
void PrDownstream::getPreSelection( PrDownTrack2& track ) {
  // Max Pt around 100 MeV for strange particle decay -> maximum displacement is in 1/p.
  double xPredTol = m_xPredTol2;
  if (std::abs(track.moment()) >  1e-6)
      xPredTol = m_xPredTol / fabs( track.moment() ) + m_xPredTol2;  // P dependance + overal tol.
  double yTol = xPredTol;

  m_xHits.clear();
  m_uvHits.clear();

  //== Collect all hits compatible with the extrapolation, region by region.
  if ( m_printing ) info() << "-- collect hits with tolerance " << xPredTol << endmsg;

  for ( int kSta = 0; m_utHitManager->maxStations() > kSta; ++kSta ) {
    for ( int kLay = 0; m_utHitManager->maxLayers() > kLay; ++kLay ) {
      for ( int kReg = 0; m_utHitManager->maxRegions() > kReg; ++kReg ) {
        const Tf::IUTHitCreator::STRegion* reg = m_utHitManager->region( kSta, kLay, kReg );
        const double yReg = track.yAtZ( reg->z() );
        if ( !reg->isYCompatible( yReg, yTol ) ) continue;
        const double xReg = track.xAtZ( reg->z() );
        if ( !reg->isXCompatible( xReg, xPredTol ) ) continue;
        Tf::UTStationHitManager<PrUTHit>::HitRange coords = m_utHitManager->hits( kSta, kLay, kReg );
        for(auto hit : coords) {
          if (hit->hit()->ignore()) continue;
          if ( hit->hit()->testStatus( Tf::HitBase::UsedByPatMatch ) ) {
            if (m_printing) debug()<<"Skip hit "<<hit->hit()->lhcbID()<<endmsg;
            continue;
          }

          const double yPos   = track.yAtZ( hit->z() );
          if ( !hit->hit()->isYCompatible( yPos, yTol ) ) continue;

          const double yTrack = track.yAtZ( 0. );
          const double tyTr   = track.slopeY();
          updateUTHitForTrack( hit, yTrack, tyTr );
          const double pos    = track.xAtZ( hit->z() );
          if ( xPredTol < fabs( pos - hit->x() ) ) continue;
          hit->hit()->setStatus( Tf::HitBase::UsedByPatDownstream, false );
          hit->setProjection( fabs( hit->x()-pos ) );
          // check if x or stereo hit
          if ( hit->hit()->lhcbID().stID().layer() ==
               hit->hit()->lhcbID().stID().station() ) {
            m_xHits.push_back( hit );
          } else {
            m_uvHits.push_back( hit );
          }
          if ( m_printing ) {
            info() << format( "  plane%2d z %8.2f x %8.2f pos %8.2f High%2d dist %8.2f",
                              hit->planeCode(), hit->z(), hit->x(), pos,
                              hit->hit()->sthit()->cluster().highThreshold(), hit->x() - pos);
            info()<<"calling debug tool now"<<endmsg;
            if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
            info()<<"done calling debug tool now"<<endmsg;
            info() << endmsg;
          }
        }
      }
    }
  }
  std::sort( m_xHits.begin(),  m_xHits.end(),  Tf::increasingByProjection<PrUTHit>() );
  std::sort( m_uvHits.begin(), m_uvHits.end(), Tf::increasingByProjection<PrUTHit>() );
}


//=========================================================================
//  Fit and remove the worst hit, as long as over tolerance
//=========================================================================
void PrDownstream::fitAndRemove ( PrDownTrack& track ) {

  if ( 2 > track.hits().size() ) return;  // no fit if single point !

  bool again;
  do {
    again = false;
    //== Fit, using the magnet point as constraint.
    double mat[6], rhs[3];
    mat[0] = 1./( track.errXMag() * track.errXMag() );
    mat[1] = 0.;
    mat[2] = 0.;
    mat[3] = 0.;
    mat[4] = 0.;
    mat[5] = 0.;
    rhs[0] = mat[0] * track.dxMagnet();//( m_magnetSave.x() - m_magnet.x() );
    rhs[1] = 0.;
    rhs[2] = 0.;
    int nbUV = 0;

    for(auto hit : track.hits()) {
      double yTrack = track.yAtZ( 0. );
      double tyTr   = track.slopeY();
      updateUTHitForTrack( hit, yTrack, tyTr );
      double dz   = hit->z() - track.zMagnet();
      double dist = track.distance( hit );
      double w    = hit->hit()->weight();
      double t    = hit->hit()->sinT();
      mat[0] += w;
      mat[1] += w * dz;
      mat[2] += w * dz * dz;
      mat[3] += w * t;
      mat[4] += w * dz * t ;
      mat[5] += w * t  * t ;
      rhs[0] += w * dist;
      rhs[1] += w * dist * dz;
      rhs[2] += w * dist * t ;
      if ( hit->hit()->lhcbID().stID().layer() !=
           hit->hit()->lhcbID().stID().station() ) nbUV++;
      if ( m_printing ) {
        info() << format( "   Plane %2d x %7.2f dist %6.3f ",
                          hit->planeCode(), hit->x(), dist );
        
        if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
        info() << endmsg;
      }
    }

    if (0 != nbUV) {
      CholeskyDecomp<double, 3> decomp(mat);
      if (!decomp) {
	track.setChisq(1e42);
	return;
      } else {
	decomp.Solve(rhs);
      }
    } else {
      CholeskyDecomp<double, 2> decomp(mat);
      if (!decomp) {
	track.setChisq(1e42);
	return;
      } else {
	decomp.Solve(rhs);
      }
      rhs[2] = 0.;
    }

    const double dx = rhs[0], dsl = rhs[1], dy = rhs[2];

    if ( m_printing ) {
      info() << format( "  dx %7.3f dsl %7.6f dy %7.3f, displY %7.2f",
                        dx, dsl, dy, track.displY() ) << endmsg;
    }

    if ( 4 > nbUV ) track.updateX( dx, dsl );
    track.setDisplY( track.displY() + dy );

    //== Remove worst hit and retry, if too far.

    double chisq = track.initialChisq();

    double maxDist = -1.;
    PrUTHits::iterator worst;
    for (PrUTHits::iterator itH = track.hits().begin();
         itH != track.hits().end(); itH++){
      PrUTHit* hit = *itH;
      double yTrack = track.yAtZ( 0. );
      double tyTr   = track.slopeY();
      updateUTHitForTrack( hit, yTrack, tyTr );
      double dist = track.distance( hit );
      if ( m_printing ) {
        info() << format( "   Plane %2d x %7.2f dist %6.3f ",
                          hit->planeCode(), hit->x(), dist );
        if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
        info() << endmsg;
      }
      chisq += dist * dist * hit->hit()->weight();
      if ( maxDist < fabs(dist) ) {
        maxDist = fabs( dist );
        worst = itH;
      }
    }
    if ( 2 < track.hits().size() ) chisq /= (track.hits().size() - 2 );
    track.setChisq( chisq );
 if(m_removeHits) {
    if ( 5. * m_maxDistance > maxDist && 0 < nbUV ) {
      // remove if Y is incompatible
      for (PrUTHits::iterator itH = track.hits().begin();
           itH != track.hits().end(); itH++){
        PrUTHit* hit = *itH;
        double yTrack = track.yAtZ( hit->z() );
        if ( yTrack > hit->hit()->yMin() && yTrack < hit->hit()->yMax() ) continue;
        if ( m_printing ) {
          info() << "   remove Y incompatible hit measure = " << hit->x()
              << " : y " << yTrack << " min " << hit->hit()->yMin()
              << " max " << hit->hit()->yMax() << endmsg;
        }
        track.hits().erase( itH );
        maxDist = 0.;  // avoid deleting worst, may be already gone...
        if ( 2 < track.hits().size() ) again = true;
        break;
      }
      if ( again ) continue;
    }

    if ( m_maxDistance < maxDist ) {
      track.hits().erase( worst );
      if ( 2 < track.hits().size() ) {
        again = true;
        if ( m_printing ) info() << "   remove worst and retry" << endmsg;
        continue;
      }
    }

    if ( m_printing ) {
      info() << format( "  ---> chi2 %7.2f maxDist %7.3f tol %7.3f", chisq, maxDist, m_maxDistance ) << endmsg;
    }
  }
  } while (again);
}


//=========================================================================
//  Collect the hits in the other x layer
//=========================================================================
void PrDownstream::findMatchingHits( PrUTHits& MatchingXHits, PrDownTrack2& track, int plane ) {
  MatchingXHits.clear();
  // reconsider this; HERE
  //search window = const1/momentum + const2
  const double tol = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolX)) ?
    m_maxWindow : (m_tolX + m_tolMomentum / track.state()->p());
  for(auto hit : m_xHits) {
      if ( !m_legacy && 8!=hit->mask() ) continue ;
    if ( plane != hit->planeCode() ) {
      const double adist = std::abs(track.distance( hit ));
      if ( adist < tol ) MatchingXHits.push_back( hit );
    }
  }
  std::sort( MatchingXHits.begin(), MatchingXHits.end(),
	  Tf::increasingByZ<PrUTHit>() );
}

//=========================================================================
//  Add the UV hits.
//=========================================================================
void PrDownstream::addUVHits ( PrDownTrack& track, std::vector< PrDownTrack >* UVtracks ) {
  //search window = const1/momentum + const2

  //first, clone the old track
  const double tol = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolUV)) ?
    m_maxWindow : (m_tolUV + m_tolMomentum / track.state()->p());
  PrUTHits uhits;
  PrUTHits vhits;

  //TODO: We might reconsider if we can vetoe some of the obviously wrong
  //combinations right here to save time spend in combinatorics later
  for(auto hit : m_uvHits) {
    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();
    updateUTHitForTrack( hit, yTrack, tyTr );
    const double adist = std::abs( track.distance( hit ) );


    if ( adist < tol && 2==hit->mask() ) uhits.push_back( hit );
    if ( adist < tol && 4==hit->mask() ) vhits.push_back( hit );
  }

  for(auto uhit : uhits) {
      for(auto vhit : vhits) {
          //first, clone the old track
          PrDownTrack newTrack( track );
          newTrack.hits().push_back( uhit );
          newTrack.hits().push_back( vhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

  if(uhits.size()==0) {
      for(auto vhit : vhits) {
          //first, clone the old track
          PrDownTrack newTrack( track );
          newTrack.hits().push_back( vhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

  if(vhits.size()==0) {
      for(auto uhit : uhits) {
          //first, clone the old track
          PrDownTrack newTrack( track );
          newTrack.hits().push_back( uhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

}

void PrDownstream::addUVHits2 ( PrDownTrack2& track, std::vector< PrDownTrack2 >* UVtracks ) {
  //search window = const1/momentum + const2

  //first, clone the old track
  const double tol = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolUV)) ?
    m_maxWindow : (m_tolUV + m_tolMomentum / track.state()->p());
  PrUTHits uhits;
  PrUTHits vhits;

  //TODO: We might reconsider if we can vetoe some of the obviously wrong
  //combinations right here to save time spend in combinatorics later
  for(auto hit : m_uvHits) {
    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();
    updateUTHitForTrack( hit, yTrack, tyTr );
    const double adist = std::abs( track.distance( hit ) );


    if ( adist < tol && 2==hit->mask() ) uhits.push_back( hit );
    if ( adist < tol && 4==hit->mask() ) vhits.push_back( hit );
  }

  for(auto uhit : uhits) {
      for(auto vhit : vhits) {
        //first, clone the old track
        PrDownTrack2 newTrack( track );
          newTrack.hits().push_back( uhit );
          newTrack.hits().push_back( vhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

  if(uhits.size()==0) {
      for(auto vhit : vhits) {
          //first, clone the old track
          PrDownTrack2 newTrack( track );
          newTrack.hits().push_back( vhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

  if(vhits.size()==0) {
      for(auto uhit : uhits) {
          //first, clone the old track
          PrDownTrack2 newTrack( track );
          newTrack.hits().push_back( uhit );
          newTrack.sortFinalHits();
          UVtracks->push_back( newTrack );
      }
  }

}

//=========================================================================
//  Add the UV hits.
//=========================================================================
void PrDownstream::addUVHits ( PrDownTrack& track ) {
  //search window = const1/momentum + const2
  const double tol = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolUV)) ?
    m_maxWindow : (m_tolUV + m_tolMomentum / track.state()->p());
  for(auto hit : m_uvHits) {
    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();
    updateUTHitForTrack( hit, yTrack, tyTr );
    const double adist = std::abs( track.distance( hit ) );
    if ( adist < tol ) track.hits().push_back( hit );
  }
  track.sortFinalHits();
}

void PrDownstream::addUVHits2 ( PrDownTrack2& track ) {
  //search window = const1/momentum + const2
  const double tol = (std::abs(track.state()->p() / m_tolMomentum) < 1. / (m_maxWindow - m_tolUV)) ?
    m_maxWindow : (m_tolUV + m_tolMomentum / track.state()->p());
  for(auto hit : m_uvHits) {
    const double yTrack = track.yAtZ( 0. );
    const double tyTr   = track.slopeY();
    updateUTHitForTrack( hit, yTrack, tyTr );
    const double adist = std::abs( track.distance( hit ) );
    if ( adist < tol ) track.hits().push_back( hit );
  }
  track.sortFinalHits();
  
  
}


//=========================================================================
//  Check if the new candidate is better than the old one
//=========================================================================
bool PrDownstream::legitCandidate( PrDownTrack2& track ){

  const int nbMeasureOK = track.hits().size();

  //== Enough mesures to have Chi2/ndof.
  if ( 3 > nbMeasureOK ) {
    if ( m_printing ) info() << " === not enough points" << endmsg;
    return false;
  }

  //== Good enough Chi2/ndof
  if ( m_maxChisq < track.chisq() ) {
    if ( m_printing ) info() << " === Chisq too big " << track.chisq() << endmsg;
    return false;
  }
  /*
  const double deltaP = track.moment() * track.state()->qOverP() - 1.;
  //== Compatible momentum
  if ( m_deltaP < fabs(deltaP) ) {
    if ( m_printing ) info() << " === Deltap too big " << deltaP << endmsg;
    if ( !m_magnetOff )return false;
    return false;
  }
  */
  //== Count if enough with high threshold
  int nbHigh = 0;
  int nbUsed = 0;

  for(auto hit : track.hits()) {
    if ( hit->hit()->sthit()->cluster().highThreshold() ) ++nbHigh;
    if ( hit->hit()->testStatus( Tf::HitBase::UsedByPatDownstream )  ) ++nbUsed;
  }
  if ( 2 > nbHigh ) {
    if ( m_printing ) info() << " === not enough high threshold points" << endmsg;
    return false;
  }
  if ( nbMeasureOK == nbUsed ) {
    if ( m_printing ) info() << " === is a clone" << endmsg;
    //return false;
  }

  const double momentum = std::abs(track.moment());
  const double tx2 = track.slopeX()*track.slopeX();
  const double ty2 = track.slopeY()*track.slopeY();
  const double sinTrack = sqrt( 1. - 1./(1.+tx2 + ty2) );
  const double pt = sinTrack*momentum;

  if (momentum<m_minMomentum) return false;
  if (pt<m_minPt) return false;

  return true;
}

//=========================================================================
//  Check if the new candidate is better than the old one
//=========================================================================
bool PrDownstream::acceptCandidate( PrDownTrack2& track, PrDownTrack2&  bestTrack,
                                     int& maxPoints, double& minChisq ){

  const int nbMeasureOK = track.hits().size();

  //== Enough mesures to have Chi2/ndof.
  if ( 3 > nbMeasureOK ) {
    if ( m_printing ) info() << " === not enough points" << endmsg;
    return false;
  }

  //== Good enough Chi2/ndof
  if ( m_maxChisq < track.chisq() ) {
    if ( m_printing ) info() << " === Chisq too big " << track.chisq() << endmsg;
    return false;
  }
  /*
  const double deltaP = track.moment() * track.state()->qOverP() - 1.;
  //== Compatible momentum
  if ( m_deltaP < fabs(deltaP) ) {
    if ( m_printing ) info() << " === Deltap too big " << deltaP << endmsg;
    if ( !m_magnetOff )return false;
    return false;
  }
  */
  //== Longest -> Keeep it
  if ( maxPoints <= nbMeasureOK ) {
    //== Same : Keep if better Chi2
    if ( minChisq < track.chisq() ) {
      if ( m_printing ) info() << " === Chisq bigger than previous " << track.chisq() << endmsg;
      return false;
    }
  } else {
    if ( m_printing ) info() << " === less points than previous" << endmsg;
    return false;
  }

  //== Count if enough with high threshold
  int nbHigh = 0;
  int nbUsed = 0;

  for(auto hit : track.hits()) {
    if ( hit->hit()->sthit()->cluster().highThreshold() ) ++nbHigh;
    if ( hit->hit()->testStatus( Tf::HitBase::UsedByPatDownstream )  ) ++nbUsed;
  }
  if ( 2 > nbHigh ) {
    if ( m_printing ) info() << " === not enough high threshold points" << endmsg;
    return false;
  }
  if ( nbMeasureOK == nbUsed ) {
    if ( m_printing ) info() << " === is a clone" << endmsg;
    return false;
  }

  if ( m_printing ) {
    info() << format( "  *** Good candidate ***  slope%8.2f displX%8.2f Y%8.2f Chi2%8.2f",
                      //1000.*track.slopeX(), track.displX(), track.displY(), track.chisq() );
                      1000.*track.slopeX(), track.xAtZ(0), track.yAtZ(0), track.chisq() );
    if ( 0 == bestTrack.hits().size() ) {
      info() << " OK"  << endmsg;
    } else {
      info() << " Better" << endmsg;
    }
  }

  //== Better candidate.
  minChisq = track.chisq();
  maxPoints = nbMeasureOK;
  if ( maxPoints > 4 ) maxPoints = 4;

  //== calculate pt and p

  const double momentum = std::abs(track.moment());
  const double tx2 = track.slopeX()*track.slopeX();
  const double ty2 = track.slopeY()*track.slopeY();
  const double sinTrack = sqrt( 1. - 1./(1.+tx2 + ty2) );
  const double pt = sinTrack*momentum;

  if (momentum<m_minMomentum) return false;
  if (pt<m_minPt) return false;



  for( auto hit : bestTrack.hits() ) {
    hit->hit()->setStatus( Tf::HitBase::UsedByPatDownstream, false );
  }

  track.sortFinalHits();

  for( auto hit : track.hits() ) {
    hit->hit()->setStatus( Tf::HitBase::UsedByPatDownstream, true );
  }


  if ( m_printing ) {
    for(auto hit : track.hits()) {
      LHCb::STChannelID icID = hit->hit()->lhcbID().stID();
      double xCoord = hit->x() ;
      info() << "      UT Clus "
             << format( "(S%1d,L%1d,R%2d,S%2d,S%3d) x%7.1f  dist%7.3f High %1d",
                        icID.station(), icID.layer(), icID.detRegion(),
                        icID.sector(), icID.strip(), xCoord,
                        track.distance( hit ), hit->hit()->sthit()->cluster().highThreshold() ) ;
      if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
      info() << endmsg;
    }
  }
  return true;
}

//=========================================================================
//  Store Track
//=========================================================================
void PrDownstream::storeTrack( PrDownTrack2& track, LHCb::Tracks* finalTracks, LHCb::Track* tr, LHCb::Tracks* updatedSeeds){
  //=== Store the tracks
  //== new LHCb::Track
  info()<<"starting store track"<<endmsg;
  
  LHCb::Track* dTr = new LHCb::Track();
  //== add ancestor
  dTr->addToAncestors( tr );
  //== Adjust flags
  dTr->setType(         LHCb::Track::Downstream  );
  dTr->setHistory(      LHCb::Track::PrDownstream   );
  dTr->setPatRecStatus( LHCb::Track::PatRecIDs   );
  //== add Seed LHCbIDs
  dTr->addSortedToLhcbIDs(  tr->lhcbIDs()      );
  //== add new LHCbIDs
  for ( PrUTHits::iterator itH = track.hits().begin();
        track.hits().end() != itH; ++itH ){
    dTr->addToLhcbIDs( (*itH)->hit()->lhcbID() );
  }

  dTr->setChi2PerDoF(track.chisq());

  //== add states
  // S. Stahl added 3 T-States

  // check for FPE and magnet off
  double QOverP = 1e-5;
  if (std::abs(track.moment()) >  1e-6){
      QOverP = 1.0 / track.moment();
  }


  //== create a state at zUTa
  LHCb::State utState;
  utState.setState( track.xAtZ( m_zUTa ),
                    track.yAtZ( m_zUTa ),
                    m_zUTa,
                    track.slopeX(),
                    track.slopeY(),
                    QOverP );
  Gaudi::TrackSymMatrix cov;
  cov(0,0) = m_stateErrorX2;
  cov(1,1) = m_stateErrorY2;
  cov(2,2) = m_stateErrorTX2;
  cov(3,3) = m_stateErrorTY2;
  double errQOverP = m_stateErrorP * QOverP;
  cov(4,4) = errQOverP * errQOverP;

  utState.setCovariance( cov );
  dTr->addToStates( utState );

  //== add seed states
  std::vector<LHCb::State*> newstates;
  newstates.reserve(3);
  newstates.push_back(tr->closestState(StateParameters::ZBegT).clone());
  newstates.push_back(tr->closestState(StateParameters::ZMidT).clone());
  // make sure we don't include same state twice
  if (std::abs(newstates[newstates.size() - 2]->z() -
               newstates.back()->z()) < 300.) {
    delete newstates.back();
    newstates.pop_back();
  }
  newstates.push_back(tr->closestState(StateParameters::ZEndT).clone());
  // make sure we don't include same state twice
  if (std::abs(newstates[newstates.size() - 2]->z() -
               newstates.back()->z()) < 300.) {
    delete newstates.back();
    newstates.pop_back();
  }

  // adjust q/p and its uncertainty
  //BOOST_FOREACH(LHCb::State* st, newstates) {
  for(auto st : newstates) {
    st->covariance()(4,4) = errQOverP * errQOverP;
    st->setQOverP(QOverP);
  }
  dTr->addToStates(newstates);
  //write all the extra info to the track before storing
  if(m_withDebugTool){
    std::cout<<"storing track information"<<std::endl;
    m_debugTool->writeExtraInfoToDownstreamTrack(*dTr);
    m_debugTool->writeExtraInfoToSeed(*dTr,*tr);
  }
  
  
  //== Save track
  finalTracks->insert( dTr );
  if(m_withDebugTool){
    updatedSeeds->insert(tr);
  }
  
}



//=============================================================================
// Fit the projection in the zx plane, one hit in each x layer
//=============================================================================
void PrDownstream::fitXProjection( PrUTHits& xHits, PrDownTrack& track, PrUTHit* firstHit ){
  // Catch if there is no second hit in other station
  track.hits().push_back( firstHit );
  PrDownTrack bestTrack( track );
  PrDownTrack newTrack( track );
  double minChisq = 1000000;
  for( auto hit : xHits) {
    newTrack.startNewXCandidate( firstHit );
    newTrack.hits().push_back( hit );
    fitAndRemove( newTrack );
    if ( newTrack.chisq() < minChisq ) {
      minChisq = newTrack.chisq();
      bestTrack = newTrack;
    }
  }
  track = bestTrack;
}
//added this, update the header
void PrDownstream::fitXProjection2( PrUTHits& xHits, PrDownTrack2& track, PrUTHit* firstHit ){
  //change this to instead fit a line to two x points, then compare to the xmagnet position.
  // I use a simple chi2 for comparison.
  // should save time in matrix inversion, etc.
  //catch if there is no second hit, as before
  track.hits().push_back(firstHit);
  PrDownTrack2 bestTrack(track);
  PrDownTrack2 newTrack(track);
  double distSave(10000.);
  double alpha_x2 = -0.0062195;
  double beta_x2 = 1843.;
  double q = (track.moment()>0)?1:-1;
  
  double tolerance = 30;//mm, just a guess for now.
  for( auto hit: xHits ){
    //adjust x2 hit to be on the new line.
    double beginX = firstHit->x();
    double endX = hit->x();
    if(hit->planeCode()==3){
      endX-= alpha_x2 *q + beta_x2/newTrack.moment();
    }
    else{
      beginX-= alpha_x2 *q + beta_x2/newTrack.moment();
    }
    
    //form x candidate
    
    double slopex = (beginX - endX)/(firstHit->z() - hit->z());
    newTrack.setSlopeX(slopex);
    double x0 = (beginX + endX - slopex*(firstHit->z() + hit->z()))/2.;
    newTrack.setXat0(x0);
    double xAtZmagTest = x0+ slopex*track.zMagnet();
    double dist = xAtZmagTest - track.xMagnet();
    double int_chi2 = (dist*dist)/(tolerance*tolerance);
    
    if(int_chi2<10 && fabs(dist)<fabs(distSave)){
      distSave = dist;
      newTrack.hits().push_back(hit);
      bestTrack = newTrack;
    }
  }
  track = bestTrack;
}

//=============================================================================
// This is needed for tracks which have more than one x hit in one layer
//=============================================================================

void PrDownstream::addOverlapRegions( PrDownTrack2& track ){
  bool hitAdded = false;
  for( auto hit : m_xHits ) {
    if ( hit->hit()->testStatus( Tf::HitBase::UsedByPatDownstream )) continue;
    if ( m_maxDistance > std::abs( track.distance( hit ) ) ) {
      double yTrack = track.yAtZ( hit->z() );
      if ( yTrack < hit->hit()->yMin() || yTrack > hit->hit()->yMax() ) continue;
      track.hits().push_back( hit );
      hitAdded = true;
    }
  }
  if ( hitAdded ) {
      track.sortFinalHits();
      //fitAndRemove( track );
      fitAndRemove2( track );
  }
}



//=========================================================================
//  Fit and remove the worst hit, as long as over tolerance
//=========================================================================
void PrDownstream::fitAndRemove2 ( PrDownTrack2& track ) {

  if ( 2 > track.hits().size() ) return;  // no fit if single point !
  double ty_beta = 1.869;
  double y_beta = -3678;
  double zput = 2320.;
  double tyT = track.state()->ty();
  tyT-=ty_beta/track.moment();
  double seed_yat0 = track.state()->y() - tyT*track.state()->z();
  double yT = seed_yat0 + tyT*zput;
  yT-=y_beta/track.moment();//signed momentum.
  
  double alphaY =11.87;
  double betaY = -6.4035e4;
  double gammaY = 4.8315e8;
  
  double alphaTY = 0.0043405;
  double betaTY = -28.455;
  double gammaTY = 1.7005e5;  
  
  double error_yy = alphaY + betaY/fabs(track.moment()) + gammaY/(track.moment()*track.moment());
  error_yy*=error_yy;
  double error_tyty = alphaTY + betaTY/fabs(track.moment()) + gammaTY/(track.moment()*track.moment());
  error_tyty*=error_tyty;
  double rho = -0.766;
  double error_tyy = rho*sqrt(error_tyty*error_yy);
  
  double guessY = yT - tyT*zput +0.001;//add an offset so they're not exactly the same
  double guessSY = tyT + 0.0001;//same here.
  track.setSlopeY(guessSY);
  track.setYat0(guessY);
  
  //upstream or downstream of center of UT
  double x_upstream = 0;
  double z_upstream = 2700.;
  
  double x_downstream = 0;
  double z_downstream = 2300.;
  double internal_chi2 = -1.;
  double q = (track.moment()>0)? 1: -1;
  //for adjustment of positions to fall on line.
  double alpha_v = -0.00418;
  double beta_v = 1309.;
  double alpha_x2 = -0.0062195;
  double beta_x2 = 1843.;

  double alpha_sigmav = 0.0070945;
  double beta_sigmav = 213.7;
  double gamma_sigmav = 1.0178e6;
  
  double alpha_sigmax2 = 0.009321;
  double beta_sigmax2 = 249.95;
  double gamma_sigmax2 = 1.34e6;
  //use the longest lever arm we can for the initial x guess.
  for(auto hit: track.hits()){
    if(hit->planeCode()<2){
      x_upstream = hit->z()< z_upstream ? hit->x() : x_upstream;
      z_upstream = hit->z()< z_upstream ? hit->z() : z_upstream;
    }
    else{
      //x_downstream = hit->z()> z_downstream ? hit->x() : x_downstream;
      if(hit->planeCode()==2){
        x_downstream = hit->z()>z_downstream? 
          (hit->hit()->xMin()*hit->hit()->cosT() + hit->hit()->yMax()*hit->hit()->sinT() 
           -( alpha_v*q+ beta_v/track.moment()))*hit->hit()->cosT():
          x_downstream;
      }
      else{
        x_downstream = hit->z()>z_downstream?
          (hit->hit()->x() - (alpha_x2 * q + beta_x2/track.moment())):
          x_downstream;
      }
      z_downstream = hit->z()> z_downstream ? hit->z() : z_downstream;
    }
    
  }
  //std::cout<<"z_upstream = "<<z_upstream<<", x_upstream = "<<x_upstream<<std::endl;
  //  std::cout<<"z_downstream = "<<z_downstream<<",x_downstream = "<<x_downstream<<std::endl;
  //redo guesses, especially if no x2
  double guessSX= (x_downstream - x_upstream)/(z_downstream - z_upstream);
  double guessX =0.5 * ( (x_downstream + x_upstream) - guessSX * (z_upstream+z_downstream));
  track.setSlopeX(guessSX);
  track.setXat0(guessX);  
  bool again;
  //std::cout<<"will loop over "<< track.hits().size()<<" hits"<<std::endl;
  int nit = 0;
  do {
    again = false;
    //== Fit, using the t seed as a measurement.
    double mat[10], rhs[4];
    internal_chi2 = 0;
      
    /*
      The matrix which is being filled looks like
      --           --
      |  0          |
      |  1  2       |
      |  3  4  5    |
      |  6  7  8  9 |
      --           --
    */
    
    mat[0] = 0.;
    mat[1] = 0.;
    mat[2] = 0.;
    mat[3] = 0.;
    mat[4] = 0.;
    mat[5] = (-1)*error_tyty/(error_tyy*error_tyy - error_tyty*error_yy);
    mat[6] = 0.;
    mat[7] = 0.;
    mat[8] = (error_tyy - zput*error_tyty)/(error_tyy*error_tyy - error_tyty*error_yy);
    mat[9] = ((-2)*error_tyy*zput + zput*zput*error_tyty + error_yy)/((-1)*error_tyy*error_tyy + error_tyty*error_yy);
    
    rhs[0] = 0.;
    rhs[1] = 0.;
    //implemented as in python checks.
    rhs[2]=(yT - guessY - guessSY*zput)*error_tyty/(error_tyty*error_yy - error_tyy*error_tyy) + 
      (tyT - guessSY)*error_tyy/(error_tyty*error_yy - error_tyy*error_tyy);
    rhs[3]=(yT - guessY - guessSY*zput)*(zput*error_tyty - error_tyy)/(error_tyty*error_yy - error_tyy*error_tyy)+
      (guessSY- tyT)*(zput*error_tyy - error_yy)/(error_tyty*error_yy - error_tyy*error_tyy);
    
    //std::cout<<"yT = "<<yT<<", tyT = "<<tyT<<std::endl;
    
    int nbUV = 0;
    double dum = 0.;
    for(auto hit : track.hits()) {
      double yTrack = seed_yat0;
      double tyTr   = tyT;//track.slopeY();
      updateUTHitForTrack( hit, yTrack, tyTr );
      //double dz   = hit->z() - track.zMagnet();
      double z = hit->z();
      double st = hit->hit()->sinT();
      double ct = hit->hit()->cosT();
      double meas_pos = 0;
      double w    = hit->hit()->weight();

      switch (hit->planeCode()){
      case 0:
        meas_pos = hit->hit()->x();
        break;
      case 1:
        meas_pos = hit->hit()->xMin()*ct + hit->hit()->yMin()*st;
        break;
      case 2:
        meas_pos = hit->hit()->xMin()*ct + hit->hit()->yMax()*st;//theta is different for u and v layers.
        //std::cout<<"modifying mean position "<<meas_pos<<" by "<<alpha_v * q + beta_v/track.moment()<<"to be"<<std::endl;
        meas_pos-= alpha_v*q+ beta_v/track.moment();
        //std::cout<<meas_pos<<std::endl;
        dum = alpha_sigmav + beta_sigmav/track.moment() + gamma_sigmav/track.moment()/track.moment();
        w = 1./(dum * dum);
        break;
      case 3:
        meas_pos = hit->hit()->x();
        //std::cout<<"modifying mean position "<<meas_pos<<" by "<<alpha_x2 * q + beta_x2/track.moment()<<"to be"<<std::endl;
        meas_pos -= alpha_x2 * q + beta_x2/track.moment();
        //std::cout<<meas_pos<<std::endl;
        
        dum = alpha_sigmax2 + beta_sigmax2/track.moment() + gamma_sigmax2/track.moment()/track.moment();;
        w = 1./(dum * dum);
        break;
      default:
        info()<<"something wrong with setting the hit specific info"<<endmsg;
      }
               
      double guess_pos = (guessX + guessSX*z)*ct + (guessY + guessSY*z)*st;
      double dist = meas_pos - guess_pos;//track.distance( hit );
      
      mat[0] += ct*ct*w;
      mat[1] += z*ct*ct*w;
      mat[2] += ct*ct*z*z*w;
      mat[3] += ct*st*w;
      mat[4] += ct*st*z*w;
      mat[5] += st*st*w;
      mat[6] += z*ct*st*w;
      mat[7] += ct*st*z*z*w;
      mat[8] += st*st*z*w;
      mat[9] += st*st*z*z*w;

      rhs[0] += dist*w;
      rhs[1] += z*dist*w;
      rhs[2] += st*dist*w;
      rhs[3] += z*st*dist*w;
      if ( hit->hit()->lhcbID().stID().layer() !=
           hit->hit()->lhcbID().stID().station() ) nbUV++;
      if ( m_printing ) {
        info() << format( "   Plane %2d x %7.2f dist %6.3f ",
                          hit->planeCode(), hit->x(), dist );
        if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
        info() << endmsg;
      }
    }//loop on hits
    
    if (0 != nbUV) {
      CholeskyDecomp<double, 4> decomp(mat);
      //std::cout<<"!decomp = "<<!decomp<<std::endl;
      if (!decomp) {
        track.setChisq(1e42);
        //info()<<"! decomp in 4x4 matrix."<<endmsg;
        return;
      } else {
        decomp.Solve(rhs);
      }
    } 
    else {
      /* CholeskyDecomp<double, 2> decomp(mat);
      if (!decomp) {
        track.setChisq(1e42);
        info()<<"! decomp in 2x2 matrix."<<endmsg;
        return;
      } else {
        decomp.Solve(rhs);
      }
      rhs[2] = 0.;
      rhs[3] = 0.;
      */
      return;
    }
    
    const double dx = rhs[0], dsl = rhs[1], dy = rhs[2], dsly = rhs[3];
    /*info() << format( "  dx %7.3f dsl %7.6f dy %7.3f, dSY %7.6f",
      dx, dsl, dy,dsly ) << endmsg;*/
    
    double newx = guessX + dx;
    double newsx = guessSX + dsl;
    double newy = guessY + dy;
    double newsy = guessSY + dsly;
    //std::cout<<"guessX = "<<guessX<<", guessSX = "<<guessSX<<", guessY = "<<guessY<<"guessSY = "<<guessSY<<std::endl;
    //std::cout<<"newx = "<<newx<<", newsx = "<<newsx<<", newy = "<<newy<<", newsy = "<<newsy<<std::endl;
    /*
    track.setSlopeX(guessSX);
    track.setSlopeY(guessSY);
    track.setXat0(guessX);
    track.setYat0(guessY);*/
    //track.setDisplY(track.displY()+dy);
    //track.setDisplX(track.displX()+dx);
    
    track.updateTrack(dx,dsl,dy,dsly);
    
    
    if ( m_printing ) {
      info() << format( "  dx %7.3f dsl %7.6f dy %7.3f, displY %7.2f, displ SY %7.2f",
                        dx, dsl, dy, track.yAtZ(0),dsly ) << endmsg;
    }
    
    //if ( 4 > nbUV ) track.updateX( dx, dsl );
    //track.setDisplY( track.displY() + dy );

    //== Remove worst hit and retry, if too far.
    
    //double chisq = track.initialChisq();
    internal_chi2 = (tyT - newsy)*((yT - newy - newsy*zput)*error_tyty/(error_tyy*error_tyy - error_tyty*error_yy) +
                                   (tyT - newsy)*error_yy/(error_tyty*error_yy - error_tyy*error_tyy))+
      (yT - newy - newsy*zput)*((yT - newy - newsy*zput)*error_tyty/(error_tyty*error_yy - error_tyy*error_tyy) +
                                (tyT - newsy)*error_tyy/(error_tyy*error_tyy - error_tyty*error_yy));
    debug()<<"chi2 from tseed = "<<internal_chi2<<endmsg;
    
    
    //all the info from the tseed.
    double maxDist = -1.;
    PrUTHits::iterator worst;
    
    for (PrUTHits::iterator itH = track.hits().begin();
         itH != track.hits().end(); itH++){
      PrUTHit* hit = *itH;
      
      // double yTrack = track.yAtZ( 0. );
      // double tyTr   = track.slopeY();
      // updateUTHitForTrack( hit, yTrack, tyTr );
      // double dist = track.distance( hit );
      updateUTHitForTrack(hit,newy,newsy);
      double meas_pos=0.;
      double w = hit->hit()->weight();
      double z = hit->z();
      double st = hit->hit()->sinT();
      double ct = hit->hit()->cosT();
      switch (hit->planeCode()){
      case 0:
        meas_pos = hit->hit()->x();
        break;
      case 1:
        meas_pos = hit->hit()->xMin()*ct + hit->hit()->yMin()*st;
        break;
      case 2:
        meas_pos = hit->hit()->xMin()*ct + hit->hit()->yMax()*st;
        meas_pos-= alpha_v*q+ beta_v/track.moment();
        dum = alpha_sigmav + beta_sigmav/track.moment() + gamma_sigmav/track.moment()/track.moment();
        w = 1./(dum * dum);
        break;
      case 3:
        meas_pos = hit->hit()->x();
        meas_pos -= alpha_x2 * q + beta_x2/track.moment();
        dum = alpha_sigmax2 + beta_sigmax2/track.moment() + gamma_sigmax2/track.moment()/track.moment();;
        w = 1./(dum * dum);
        break;
      default:
        info()<<"something wrong with setting the hit specific info"<<endmsg;
      }
      
      double guess_pos = (newx + newsx*z)*ct + (newy + newsy*z)*st;
      double dist = meas_pos - guess_pos;//track.distance( hit );
      //info()<<"dist for hit "<<hit->planeCode()<<"="<<dist<<endmsg;
      //double z_mag = 5393. + (-1.363e5)/fabs(track.moment()) + (1.102e9)/track.moment()/track.moment();
      //double tseed_x0=(track.state.x() - track.state.z()*track.state.tx()) ;
      //double diff_at_xmag_new = (track.);

      if ( m_printing ) {
        info() << format( "   Plane %2d x %7.2f dist %6.3f ",
                          hit->planeCode(), hit->x(), dist );
        if ( m_debugTool ) m_debugTool->debugUTCluster( info(), hit );
        info() << endmsg;
      }
      debug()<<"for hit in plane "<<hit->planeCode()<<", adding chi2 "<< dist*dist*w<<endmsg;
      
      internal_chi2 += dist * dist * w;
      
      if ( maxDist < fabs(dist) ) {
        maxDist = fabs( dist );
        worst = itH;
      }
    }
    debug()<<"Total chi2 = "<<internal_chi2<<endmsg;
    nit++;
    guessX = newx;
    guessSX = newsx;
    guessY = newy;
    guessSY = newsy;
    again = (nit >m_number_Iterations_new_chi2 ) ? false:true;
    
  
    
    //    if ( 2 < track.hits().size() ) internal_chi2 /= (track.hits().size() - 2 );
    //use chi2 not chi2/ndf
    track.setChisq( internal_chi2 );
    
    if(m_removeHits) {
      if ( 5. * m_maxDistance > maxDist && 0 < nbUV ) {
        // remove if Y is incompatible
        for (PrUTHits::iterator itH = track.hits().begin();
             itH != track.hits().end(); itH++){
          PrUTHit* hit = *itH;
          //double yTrack = track.yAtZ( hit->z() );
          double yTrack = newy + newsy*hit->z();
          if ( yTrack > hit->hit()->yMin() && yTrack < hit->hit()->yMax() ) continue;
          if ( m_printing ) {
            info() << "   remove Y incompatible hit measure = " << hit->x()
                   << " : y " << yTrack << " min " << hit->hit()->yMin()
                   << " max " << hit->hit()->yMax() << endmsg;
          }
          track.hits().erase( itH );
          maxDist = 0.;  // avoid deleting worst, may be already gone...
          if ( 2 < track.hits().size() ) again = true;
          break;
        }
        if ( again ) continue;
      }

      if ( m_maxDistance < maxDist ) {
        track.hits().erase( worst );
        if ( 2 < track.hits().size() ) {
          again = true;
          if ( m_printing ) info() << "   remove worst and retry" << endmsg;
          continue;
        }
      }

      if ( m_printing ) {
        info() << format( "  ---> chi2 %7.2f maxDist %7.3f tol %7.3f", /*chisq*/internal_chi2, maxDist, m_maxDistance ) << endmsg;
      }
    }
    
  } while (again);
  
  return ;
  
}

