// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 
#include "GaudiKernel/IRegistry.h"

#include "Linker/LinkedTo.h"

#include "Event/MCParticle.h"
#include "Event/STCluster.h"
#include "Event/Track.h"

#include "STDet/DeSTDetector.h"

// local
#include "PrDebugUTTruthTool.h"
#include <map>
#include <iostream>

//AD
#include "Linker/AllLinks.h"
#include "Event/MCTrackInfo.h"
#include "MCInterfaces/ILHCbIDsToMCParticles.h"
#include "MCInterfaces/ILHCbIDsToMCHits.h"
#include "PrKernel/PrUTHit.h"
#include "Event/MCHit.h"
#include "Event/StateParameters.h"
//-----------------------------------------------------------------------------
// Implementation file for class : PrDebugUTTruthTool
//
// 2007-10-22 : Olivier Callot
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( PrDebugUTTruthTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrDebugUTTruthTool::PrDebugUTTruthTool( const std::string& type,
                                          const std::string& name,
                                          const IInterface* parent )
  : GaudiTool ( type, name , parent )
{
  declareInterface<IPrDebugUTTool>(this);

}


StatusCode PrDebugUTTruthTool::initialize()
{
  StatusCode sc = GaudiTool::initialize();
  if (sc.isFailure()) return Error("Failed to initialize", sc);

  m_tracker = getDet<DeSTDetector>(DeSTDetLocation::UT);
  m_flags["reconstructible_asLong"] = false;
  m_flags["reconstructible_asDown"] = false;
  m_flags["getPreSelection"] = false;
  m_flags["findMatchingHits"] = false;
  m_flags["fitXProjection"] = false;
  m_flags["addUVHits"] = false;
  m_flags["fitAndRemove"] = false;
  m_flags["acceptCandidate"] = false;
  m_flags["beforeStore"] = false;
  
  m_lhcbid2mcparticles = tool<ILHCbIDsToMCParticles>("LHCbIDsToMCParticles", this);//AD
  m_lhcbid2mchits = tool<ILHCbIDsToMCHits>("LHCbIDsToMCHits", this);//AD  

  return sc;
}

//=============================================================================
// Destructor
//=============================================================================
PrDebugUTTruthTool::~PrDebugUTTruthTool() {} 
//=========================================================================
//  Print the true UT clusters associated to the specified track
//=========================================================================
void PrDebugUTTruthTool::debugUTClusterOnTrack (  const LHCb::Track* track, 
                                                  const PrUTHits::const_iterator beginCoord,
                                                  const PrUTHits::const_iterator endCoord   ) {
  //== Get the truth on this track
  std::cout<<"starting debug UT cluster on track"<<std::endl;
  std::string containerName = track->parent()->registry()->identifier();
  std::cout<<"containerName = "<<containerName<<std::endl;
  std::string linkerName = "Link/"+containerName;
  std::cout<<"linkerName = "<<linkerName<<std::endl;
  
  if ( "/Event/" == containerName.substr(0,7) ) {
    linkerName = "Link/" + containerName.substr(7);
  }
  std::cout<<"linker name after check for event is "<<linkerName<<std::endl;
  std::cout<<"exist<LHCb::LinksByKey>( linkerName ) = "<<exist<LHCb::LinksByKey>( linkerName )<<std::endl;
  
  if ( exist<LHCb::LinksByKey>( linkerName ) ) {
    LinkedTo<LHCb::MCParticle, LHCb::Track> trLink( evtSvc(), msgSvc(), containerName );
    
    LinkedTo<LHCb::MCParticle> itLink( evtSvc(), msgSvc(), LHCb::STClusterLocation::UTClusters );
    if(trLink.notFound() ||  itLink.notFound()){
      debug()<<"Cannot load linkers"<<endmsg;
      return;
    }
    
    LHCb::MCParticle* part = trLink.first( track->key() );
    
    while ( 0 != part ) {
      double momentum = part->momentum().R();
      info() << format( " MC Key %4d PID %6d  P %8.3f GeV tx %9.5f ty %9.5f",
                        part->key(), part->particleID().pid(), momentum / Gaudi::Units::GeV,
                        part->momentum().x()/ part->momentum().z(),
                        part->momentum().y()/ part->momentum().z() ) << endreq;
      for ( PrUTHits::const_iterator itH = beginCoord; endCoord != itH; ++itH ){
        const PrUTHit* hit = *itH;
        LHCb::STChannelID id = hit->hit()->lhcbID().stID();
        bool found = false;
        for( unsigned int kk = 0; hit->hit()->sthit()->cluster().pseudoSize() > kk; ++kk ) {
	  LHCb::MCParticle* clusPart = 0;
          if ( id != LHCb::STChannelID(0)) clusPart = itLink.first( id );
          while ( 0 != clusPart ) {
            if ( clusPart->key() == part->key() ) found = true;
            clusPart = itLink.next();
          }
          id = m_tracker->nextRight(id);
        }
        if ( found ) {
          double xCoord = hit->hit()->x() ;
          info() << "      UT Clus " 
                 << format( "(S%1d,L%1d,R%2d,S%2d,S%3d) x%7.1f High %1d", 
                            id.station(), id.layer(), id.detRegion(), 
                            id.sector(), id.strip(), xCoord, 
                            (*itH)->hit()->sthit()->cluster().highThreshold() ) << endreq;
        }
      }
      part = trLink.next();
    }
    
    }

}

//=========================================================================
//  Print the MC keys associated to this cluster
//=========================================================================
void PrDebugUTTruthTool::debugUTCluster( MsgStream& msg, const PrUTHit* hit ) {
  
  LinkedTo<LHCb::MCParticle> itLink( evtSvc(), msgSvc(), LHCb::STClusterLocation::UTClusters );
  std::string value = "";
  
  LHCb::STChannelID id = hit->hit()->lhcbID().stID();
  int lastKey = -1;
  for( unsigned int kk = 0; hit->hit()->sthit()->cluster().pseudoSize() > kk; ++kk ) {
    LHCb::MCParticle* part = 0;
    if (id != LHCb::STChannelID(0)) part = itLink.first( id );
    while ( 0 != part ) {
      if ( lastKey != part->key() ) {
        lastKey = part->key();
        msg << " " << lastKey;
      }
      part = itLink.next();
    }
    id = m_tracker->nextRight(id);
  }
}

//=============================================================================
// Added information by Adam to try to clean up the ntuple coding
//=============================================================================

void PrDebugUTTruthTool::recordStepInProcess(std::string step,bool result){
  m_flags[step] |= result;//don't change things if we already have the right answer.
  info()<<"Recorded for step "<<step<<" result"<<result<<endmsg;
}
//=============================================================================
// Extra method from PrAddUTHits to check the hits form the UT if it's on the track
// Added AD 12-11-13
//=============================================================================

bool PrDebugUTTruthTool::isIDOnMCParticle(LHCb::LHCbID id, LHCb::Track& track)
{
  // Method here: Get the MC Particle associated with the ID and the track.
  // If we have more than one particle associated to a track or hit, take the one
  // with the most associations. Do this by flipping the map. There are easier 
  // ways, this is just what I wrote first.
  bool printing = false;
  bool isTrue = false;
  std::map<LHCb::MCParticle*, unsigned int> idLinkMap;
  std::map<unsigned int, LHCb::MCParticle*> idLinkMapFlip;
  std::map<LHCb::MCParticle*, unsigned int> trackLinkMap;

  // -- Link the id to the MCParticle and the track to the MCParticle
  m_lhcbid2mcparticles->link(id, idLinkMap);
  if(printing){ 
    info()<<"idLinkMap = "<<endmsg;
    for(auto map_ent : idLinkMap){
      info()<<"first->p() = "<<map_ent.first->p() <<" second -->"<<map_ent.second<<endmsg;
    }
  }
  m_lhcbid2mcparticles->link(track, trackLinkMap);
  if(printing){ 
    info()<<"trackLinkMap = "<<endmsg;
    for(auto map_ent : trackLinkMap){
      info()<<"first->p() = "<<map_ent.first->p() <<" second -->"<<map_ent.second<<endmsg;
    }
}

  // -- Get the MCParticle the track corresponds to (with most of the hits)
  std::map<LHCb::MCParticle*, unsigned int>::iterator it = trackLinkMap.begin();
  
  LHCb::MCParticle* trackParticle = (*it).first;
  if(!trackParticle){ debug()<<"the track particle doesn't exits!!"<<endmsg;
    return false;
  }  
  if(printing){ 
    info()<<"original trackParticle->p() = "<<trackParticle->p()<<endmsg;
  }
  
  //take the track particle with the most associations if there are more than 2
  
  if( trackLinkMap.size()>1 ){
    std::map<unsigned int, LHCb::MCParticle*> trackLinkMapFlip;
    //flip the link map. we have more than one particle. automatically take the one with the most associations
    std::map<LHCb::MCParticle*,unsigned int>::iterator itr;    
    for(itr= trackLinkMap.begin(); itr!=trackLinkMap.end();++itr){
      trackLinkMapFlip.insert(std::pair<unsigned int, LHCb::MCParticle*>((*itr).second,(*itr).first));
    }
    std::map< unsigned int, LHCb::MCParticle*>::reverse_iterator itr2 = trackLinkMapFlip.rbegin();//sorted in increasing order
    
    trackParticle = (*itr2).second;
    if(printing){
      
      info()<<"after flip, trackParticle->p() = "<<trackParticle->p();
    }
    
    debug()<<"printing the trackParticle reversed map"<<endmsg;
    for(std::map< unsigned int, LHCb::MCParticle*>::reverse_iterator itr3 = trackLinkMapFlip.rbegin(); 
        itr3!=trackLinkMapFlip.rend();++itr3){
      debug()<<"trackLinkMapFlip.first = "<<(*itr3).first<<", trackLinkMapFlip.second = "<<(*itr3).second<<endmsg;  
    }
  }
  

  // -- Get the MCParticle the id corresponds to
  std::map<LHCb::MCParticle*, unsigned int>::iterator it2 = idLinkMap.begin();
  LHCb::MCParticle* idParticle = (*it2).first;
  if(printing){
    info()<<"original idParticle->p() = "<<idParticle->p()<<endmsg;
  }
  if(idLinkMap.size() > 1 ){
    //flip the link map. we have more than one particle. automatically take the one with the most associations
    std::map<LHCb::MCParticle*,unsigned int>::iterator itr;    
    for(itr= idLinkMap.begin(); itr!=idLinkMap.end();++itr){
      idLinkMapFlip.insert(std::pair<unsigned int, LHCb::MCParticle*>((*itr).second,(*itr).first));
    }
    std::map< unsigned int, LHCb::MCParticle*>::reverse_iterator it3 = idLinkMapFlip.rbegin();//sorted in increasing order
    idParticle = (*it3).second;//get the particle with the most MC associations.
    if(printing){
      info()<<"after flip, idParticle->p() = "<<idParticle->p();
    }
    
    debug()<< "idlinkmap has size > 1, taking particle with most associations"<<endmsg;//ad
  }

  if(!idParticle){ debug()<<"the id particle doesn't exits!!"<<endmsg;
    return false;
  }
  
  
  if(idParticle == NULL && trackParticle == NULL) return false;// AD. Bug fix for if both particles don't exist
  //AD 12-13-13
  
  const unsigned int numberAssocIDs = (*it).second;
  // -- if MCParticle corresponding to id and to track with hightest weight agree, then call it associated
  if( idParticle->key() == trackParticle->key()) isTrue = true;
 
  return isTrue;
}


bool PrDebugUTTruthTool::isTrackReconstructible(LHCb::Track& track)///there's a problem here. Come fix this first.
{
  std::map<LHCb::MCParticle*, unsigned int> trackLinkMap;
  std::map<unsigned int, LHCb::MCParticle*> trackLinkMapFlip;
  m_lhcbid2mcparticles->link(track, trackLinkMap);
  
  std::map<LHCb::MCParticle*, unsigned int>::iterator it = trackLinkMap.begin();
  if((*it).first == NULL) return false;//no particle in the map
  
  LHCb::MCParticle* part =  (*it).first;
  
  if( trackLinkMap.size()>1 ){
    //flip the link map. we have more than one particle. automatically take the one with the most associations
    std::map<LHCb::MCParticle*,unsigned int>::iterator itr;    
    for(itr= trackLinkMap.begin(); itr!=trackLinkMap.end();++itr){
      trackLinkMapFlip.insert(std::pair<unsigned int, LHCb::MCParticle*>((*itr).second,(*itr).first));
    }
    std::map< unsigned int, LHCb::MCParticle*>::reverse_iterator it2 = trackLinkMapFlip.rbegin();//sorted in increasing order
    part = (*it2).second;
  }
  //now let's double check that the MC particle is a downstream track.
  MCTrackInfo trackInfo( evtSvc(), msgSvc() );
  if ( 0 == trackInfo.fullInfo( part ) ) {return false;}
  if(!part) return false;
  bool isdown = trackInfo.hasT( part ) &&  trackInfo.hasTT( part ) ;
  bool notvelo =(!trackInfo.hasVelo(part));
  
  isdown = isdown &&(abs(part->particleID().pid())!=11);//kill the electrons
  notvelo = notvelo&&(abs(part->particleID().pid())!=11);//kill the electrons
  //all the rest of this is for checks which can be made additionally.
  bool fromKs = false;
  //stolen from PrChecker
  if ( 0 != part->originVertex() ) {
    const LHCb::MCParticle* mother =  part->originVertex()->mother();
    if ( 0 != mother ) {
      if ( 0 != mother->originVertex() ) {
        double rOrigin = mother->originVertex()->position().rho();
        if ( fabs( rOrigin ) < 5. ) {//min sep is hard coded. this is hackish
          int pid = abs( mother->particleID().pid() );
          if(310 == pid){
            fromKs = true;
          }
        }
      }
    }
  }
  bool over5 = 5000. < fabs( part->p() );
  bool over300pt = 300. < fabs( part->pt() );
  bool eta25 = (part->pseudoRapidity() > 2 && part->pseudoRapidity() < 5);
    
  //if (isdown/* && over5&&over300pt && eta25 && fromKs*/){
  //  if(setMCParticle){m_mc_part = part;}
  //}
  return (isdown && notvelo /*&& over5 && over300pt && eta25 && fromKs*/);

}

std::vector<double> PrDebugUTTruthTool::TrueHitXYZ(LHCb::LHCbID id,bool extended){
  std::vector<double> ret;
  //get out the info
  // *xmid
  // *ymid
  // *zmid
  //and with extended
  // *tx
  // *ty
  // *qop
  // *entryx
  // *entryy
  // *entryz

  ILHCbIDsToMCHits::LinkMap testMap;
  m_lhcbid2mchits->link(id,testMap);
  if(testMap.size()<1){
    ret.push_back(-999999.);//badx
    ret.push_back(-999999.);//bady
    ret.push_back(-999999.);//badz
    if(extended){
      ret.push_back(-999999.);//badtx
      ret.push_back(-999999.);//badty
      ret.push_back(-999999.);//badqop
      ret.push_back(-999999.);//badentryx
      ret.push_back(-999999.);//badentryy
      ret.push_back(-999999.);//badentryz
      ret.push_back(-999999.);//badp
    }
    
    return ret;
  }

  ILHCbIDsToMCHits::LinkMap::iterator it = testMap.begin();
  LHCb::MCHit* hit = it->first;
  if(!hit)
  {
    ret.push_back(-999999.);//badx
    ret.push_back(-999999.);//bady
    ret.push_back(-999999.);//badz
    if(extended){
      ret.push_back(-999999.);//badtx
      ret.push_back(-999999.);//badty
      ret.push_back(-999999.);//badqop
      ret.push_back(-999999.);//badentryx
      ret.push_back(-999999.);//badentryy
      ret.push_back(-999999.);//badentryz
      ret.push_back(-999999.);//badp
    }

    return ret;
  }
  
  ret.push_back(hit->midPoint().x());
  ret.push_back(hit->midPoint().y());
  ret.push_back(hit->midPoint().z());
  if(extended){
    ret.push_back(hit->dxdz());//tx
    ret.push_back(hit->dydz());//ty
    //get q over p
    const double charge = (hit->mcParticle()->particleID().threeCharge() )/3.;
    const double phit = hit->p();
    const double ppart = hit->mcParticle()->p();
    double qop = (phit < TrackParameters::lowTolerance ? 
                  (ppart < TrackParameters::lowTolerance ? 0.0 : charge/ppart) : charge/phit);
    ret.push_back(qop);//qop
    ret.push_back(hit->entry().x());//x entry position for state
    ret.push_back(hit->entry().y());//y entry position for state
    ret.push_back(hit->entry().z());//z entry position for state
    ret.push_back(phit);
    
  }
  return ret;
  
}

bool PrDebugUTTruthTool::CheckMCinIntermediateHits(std::vector<int> containerx,
                                                   std::vector<int> containeruv)
{
  int ngood_x1(0),ngood_x2(0),ngood_u(0),ngood_v(0);
  
  for(auto lhcbid : containerx){
    if(m_x1_matching_map[lhcbid] ==1) ngood_x1++;
    if(m_x2_matching_map[lhcbid] ==1) ngood_x2++;
  }
  for(auto lhcbid: containeruv){
    if(m_u_matching_map[lhcbid] ==1) ngood_u++;
    if(m_v_matching_map[lhcbid] ==1) ngood_v++;
  }
  //require 3 planes
  bool the_ans = (ngood_x1+ngood_x2+ngood_u >2 //x1,u,x2 fired
                  || ngood_x1+ngood_x2+ngood_v>2//x1,v, x2 fired
                  || ngood_x2+ngood_u+ngood_v >2//u,v,x2 fired
                  || ngood_x1+ngood_u+ngood_v >2)?true:false;//x1,u,v fired
  info()<<"In check intermediate hits, result is "<<the_ans<<endmsg;
  
  return the_ans;
}

void PrDebugUTTruthTool::printTSeedInfoShort(LHCb::Track* tr){
  int nMatchedHits = 0;
  int ntot = 0;
  for(auto id: tr->lhcbIDs()){
    if(!id.isFT())continue;//don't consider other hits on the track
    ntot++;
    if(isIDOnMCParticle(id, *tr))nMatchedHits++;
  }
  
  info()<<"Seed xyz("
        <<tr->closestState(StateParameters::ZBegT).x()<<","
        <<tr->closestState(StateParameters::ZBegT).y()<<","
        <<tr->closestState(StateParameters::ZBegT).z()<<") ,tx,ty("
        <<tr->closestState(StateParameters::ZBegT).tx()<<","
        <<tr->closestState(StateParameters::ZBegT).ty()
        <<"), nMatchedHits "<<nMatchedHits<<" / "
        <<ntot<<endmsg;
  return;
}


void PrDebugUTTruthTool::printTSeedInfo(LHCb::Track* tr){
    std::cout<<"***************************************************"<<std::endl;
    std::cout<<"                Start of Tseed Loop"<<std::endl;
    std::cout<<"***************************************************"<<std::endl;
    std::cout<<"Got the TSeed, parameters are"<<std::endl;
    std::cout<<std::setw(15)<<"Begin T"<<std::setw(15)<<"End T"<<std::endl;
    std::cout<<"-----------------------------------------------------------"<<std::endl;
    std::cout<<"x"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).x()<<std::setw(15)
             <<tr->closestState(StateParameters::ZEndT).x()<<std::endl;
    std::cout<<"y"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).y()<<std::setw(15)
             <<tr->closestState(StateParameters::ZEndT).y()<<std::endl;
    std::cout<<"z"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).z()<<std::setw(15)
             <<tr->closestState(StateParameters::ZEndT).z()<<std::endl;
    std::cout<<"tx"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).tx()<<std::setw(15)
             <<tr->closestState(StateParameters::ZEndT).tx()<<std::endl;
    std::cout<<"ty"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).ty()<<std::setw(15)
             <<tr->closestState(StateParameters::ZEndT).ty()<<std::endl;
    //std::cout<<"p_of_state"<<std::setw(15)<<tr->closestState(StateParameters::ZBegT).p()<<std::setw(15)<<std::endl;
    std::cout<<std::endl;
}

bool PrDebugUTTruthTool::isPreselGood(std::vector<LHCb::LHCbID> x1ids,
                                      std::vector<LHCb::LHCbID> uids,
                                      std::vector<LHCb::LHCbID> vids,
                                      std::vector<LHCb::LHCbID> x2ids,
                                      LHCb::Track * track){
  //just a wrapper to make the preselection prettier.
  //also sets the maps for every other check, since every other check
  //will be a subset of the preselection hits.
  int nx1_good(0),nx2_good(0),nu_good(0),nv_good(0);
  for(auto x1id: x1ids){
    m_x1_matching_map[x1id.lhcbID()] = isIDOnMCParticle(x1id,*track);
    if(m_x1_matching_map[x1id.lhcbID()]==1) nx1_good++;
  }
  for(auto uid: uids){ 
    m_u_matching_map[uid.lhcbID()]=isIDOnMCParticle(uid,*track);
    if(m_u_matching_map[uid.lhcbID()]==1) nu_good++;
  }
  for(auto vid: vids){ 
    m_v_matching_map[vid.lhcbID()]=isIDOnMCParticle(vid,*track);
    if(m_v_matching_map[vid.lhcbID()]==1) nv_good++;
  }
  for(auto x2id: x2ids){ 
    m_x2_matching_map[x2id.lhcbID()]=isIDOnMCParticle(x2id,*track);
    if(m_x2_matching_map[x2id.lhcbID()]==1) nx2_good++;
  }
  std::cout<<"In Checking Preselection"<<std::endl;
  std::cout<<"number of good hits: x1 "<<nx1_good<<", u "<<nu_good<<", v "<<nv_good<<", x2 "<<nx2_good
           <<std::endl;
  bool the_ans = (nx1_good+nx2_good+nu_good >2 //x1,u,x2 fired
                  || nx1_good+nx2_good+nv_good>2//x1,v, x2 fired
                  || nx2_good+nu_good+nv_good >2//u,v,x2 fired
                  || nx1_good+nu_good+nv_good >2)?true:false;//x1,u,v fired
  std::cout<<"result from is presel good = "<<the_ans<<std::endl;
  
  return the_ans;
}

bool PrDebugUTTruthTool::CheckMCinTrackHits(std::vector<int> ids)
{
  int nx1_good(0),nx2_good(0),nu_good(0),nv_good(0);
  for(auto id: ids){
    if(m_x1_matching_map[id] ==1) nx1_good++;
    if(m_x2_matching_map[id] ==1) nx2_good++;
    if(m_u_matching_map[id] ==1) nu_good++;
    if(m_v_matching_map[id] ==1) nv_good++;
  }
  //require 3 planes
  bool the_ans = (nx1_good+nx2_good+nu_good >2 //x1,u,x2 fired
                  || nx1_good+nx2_good+nv_good>2//x1,v, x2 fired
                  || nx2_good+nu_good+nv_good >2//u,v,x2 fired
                  || nx1_good+nu_good+nv_good >2)?true:false;//x1,u,v fired
  return the_ans;
}

void PrDebugUTTruthTool::PrintHitTable(PrUTHit* hit)
{
  
  //info()<<"Calling Print Hit Table"<<endmsg;
  std::cout<<"................................."<<std::endl;
  
  //print all the information associated to a hit.
  std::cout<<"Printing Table of Hit Information"<<std::endl;           
    /*<<"xTrue"<<std::setw(15)
        <<"yTrue"<<std::setw(15)
        <<"zTrue"<<std::endl;*/
  std::vector<double> dummy = TrueHitXYZ(hit->hit()->lhcbID(),true);
  std::cout<<std::endl;
  //  std::cout<<std::setfill('-')<<std::setw(80)<<"-"<<std::endl;
  std::cout<<"Plane"<<std::setw(15)<<hit->planeCode()<<std::endl
           <<"length"<<std::setw(15)<<hit->hit()->length()<<std::endl
           <<"zAtYEq0"<<std::setw(15)<<hit->hit()->zAtYEq0()<<std::endl
           <<"z(hit)"<<std::setw(15)<<hit->hit()->z()<<std::endl
           <<"zmid"<<std::setw(15)<<hit->hit()->zMid()<<std::endl
           <<"yBegin"<<std::setw(15)<<hit->hit()->yBegin()<<std::endl
           <<"yEnd"<<std::setw(15)<<hit->hit()->yEnd()<<std::endl
           <<"y"<<std::setw(15)<<hit->hit()->y()<<std::endl
           <<"yMin"<<std::setw(15)<<hit->hit()->yMin()<<std::endl
           <<"yMid"<<std::setw(15)<<hit->hit()->yMid()<<std::endl
           <<"yMax"<<std::setw(15)<<hit->hit()->yMax()<<std::endl
           <<"isX"<<std::setw(15)<<hit->hit()->isX()<<std::endl
           <<"xAtYEq0"<<std::setw(15)<<hit->hit()->xAtYEq0()<<std::endl
           <<"x" <<std::setw(15)<<hit->hit()->x()<<std::endl
           <<"xT"<<std::setw(15)<<hit->hit()->xT()<<std::endl
           <<"xMin"<<std::setw(15)<<hit->hit()->xMin()<<std::endl
           <<"xMid"<<std::setw(15)<<hit->hit()->xMid()<<std::endl
           <<"xMax"<<std::setw(15)<<hit->hit()->xMax()<<std::endl
           <<"dxDy"<<std::setw(15)<<hit->hit()->dxDy()<<std::endl
  
    //           <<dummy[0]<<std::setw(15)
    //           <<dummy[1]<<std::setw(15)
    //           <<dummy[2]
           <<std::endl;

  std::cout<<"True Hit info"<<std::endl
           <<"xmid"<<std::setw(15)<<dummy[0]<<std::endl
           <<"ymid"<<std::setw(15)<<dummy[1]<<std::endl
           <<"zmid"<<std::setw(15)<<dummy[2]<<std::endl
           <<"tx"<<std::setw(15)<<dummy[3]<<std::endl
           <<"ty"<<std::setw(15)<<dummy[4]<<std::endl
           <<"qop"<<std::setw(15)<<dummy[5]<<std::endl
           <<"entryx"<<std::setw(15)<<dummy[6]<<std::endl
           <<"entryy"<<std::setw(15)<<dummy[7]<<std::endl
           <<"entryz"<<std::setw(15)<<dummy[8]<<std::endl
           <<"true_p"<<std::setw(15)<<dummy[9]<<std::endl;
  std::cout<<"angles are"<<std::endl;
  std::cout<<"cosTheta"<<std::setw(15)<<hit->hit()->cosT()<<std::endl;
  std::cout<<"sinTheta"<<std::setw(15)<<hit->hit()->sinT()<<std::endl;
  std::cout<<"errors"<<std::setw(15)<<hit->hit()->weight()<<std::endl;
  std::cout<<"Matched"<<std::setw(15)<<(m_x1_matching_map[hit->hit()->lhcbID().lhcbID()] ||
                                        m_u_matching_map[hit->hit()->lhcbID().lhcbID()] ||
                                        m_v_matching_map[hit->hit()->lhcbID().lhcbID()] ||
                                        m_x2_matching_map[hit->hit()->lhcbID().lhcbID()] )<<std::endl;
  
  std::cout<<"................................................"<<std::endl;


  return;
}
void PrDebugUTTruthTool::PrintHitTableShort(const PrUTHit* hit){
  //int id = hit->hit()->lhcbID().lhcbID();
  //  info()<<"printing stream to check"<<endmsg;
  //  hit->hit()->lhcbID().fillStream(std::cout);
  //  info()<<"done"<<endmsg;
  //info()<<"testing id in print hit table short = "<<id<<endmsg;
  
  bool truth_matched = ((m_x1_matching_map[hit->hit()->lhcbID().lhcbID()] ==1)||
                        (m_x2_matching_map[hit->hit()->lhcbID().lhcbID()] ==1)||
                        (m_u_matching_map[hit->hit()->lhcbID().lhcbID()] ==1)||
                        (m_v_matching_map[hit->hit()->lhcbID().lhcbID()] ==1))?1:0;
  
  info()<<"Hit ("<< hit->hit()->lhcbID().lhcbID()<<") Plane "<<hit->planeCode()
        <<" x,y,z("<<hit->hit()->x(hit->hit()->y())<<","<<hit->hit()->y()<<","<<hit->hit()->z(hit->hit()->y())
        <<"), truthMatched "<<truth_matched<<endmsg;
  
  return;
}

StatusCode PrDebugUTTruthTool::writeExtraInfoToDownstreamTrack(LHCb::Track& dsTrack){
  //write the extra info of the algorithm to the downstream track
  //for use with PrDownstreamChecker.
  //first write all the info from the flags
  if(!dsTrack.addInfo(1801,m_flags["getPreSelection"]))Warning("Could not add info 1801 to DS track").ignore();
  if(!dsTrack.addInfo(1802,m_flags["findMatchingHits"]))Warning("Could not add info 1802 to DS track").ignore();
  if(!dsTrack.addInfo(1803,m_flags["fitXProjection"]))Warning("Could not add info 1803 to DS track").ignore();
  if(!dsTrack.addInfo(1804,m_flags["addUVHits"]))Warning("Could not add info 1804 to DS track").ignore();
  if(!dsTrack.addInfo(1805,m_flags["fitAndRemove"]))Warning("Could not add info 1805 to DS track").ignore();
  if(!dsTrack.addInfo(1806,m_flags["acceptCandidate"]))Warning("Could not add info 1806 to DS track").ignore();
  if(!dsTrack.addInfo(1807,m_flags["beforeStore"]))Warning("Could not add info 1807 to DS track").ignore();
  std::cout<<"writing extra info to downstream track. Test info(1801) = "<<dsTrack.info(1801,-999)<<std::endl;
  //next, categorize all the hits used in the downstream track.
  int n_good_x1(0),n_good_x2(0),n_good_u(0),n_good_v(0);
  int n_x1(0),n_x2(0),n_u(0),n_v(0);
  for(auto id: dsTrack.lhcbIDs()){
    //only categorize UT hits
    if(!id.isUT())continue;
    if(id.stID().station() == id.stID().layer()){//xhits
      //goes as 
      //x1 = station(1) layer(1)
      //u  = station(1) layer(2)
      //v  = station(2) layer(1)
      //x2 = station(2) layer(2)
      if(id.stID().station()==1){//x1
        n_x1++;
        if(m_x1_matching_map[id.lhcbID()])n_good_x1++;
      }
      if(id.stID().station()==2){//x2
        n_x2++;
        if(m_x2_matching_map[id.lhcbID()])n_good_x2++;
      }
    }
    else{
      if(id.stID().station()==1){//x1
        n_u++;
        if(m_u_matching_map[id.lhcbID()])n_good_u++;
      }
      if(id.stID().station()==2){//x2
        n_v++;
        if(m_v_matching_map[id.lhcbID()])n_good_v++;
      }
    }
  }
  if(!dsTrack.addInfo(1900, n_good_x1))Warning("Couldn't add good x1 hits on track").ignore();
  if(!dsTrack.addInfo(1901, n_good_u))Warning("Couldn't add good u hits on track").ignore();
  if(!dsTrack.addInfo(1902, n_good_v))Warning("Couldn't add good v hits on track").ignore();
  if(!dsTrack.addInfo(1903, n_good_x2))Warning("Couldn't add good x2 hits on track").ignore();
  if(!dsTrack.addInfo(1904, n_x1))Warning("Couldn't add total x1 hits on track").ignore();
  if(!dsTrack.addInfo(1905, n_u))Warning("Couldn't add total u hits on track").ignore();
  if(!dsTrack.addInfo(1906, n_v))Warning("Couldn't add total v hits on track").ignore();
  if(!dsTrack.addInfo(1907, n_x2))Warning("Couldn't add total x2 hits on track").ignore();
  //add chi2 information from fit and remove?
  return StatusCode::SUCCESS; 
}

StatusCode PrDebugUTTruthTool::writeExtraInfoToSeed(LHCb::Track& dTr, LHCb::Track& seed){
  //first, get the information we wrote to the downstream track and write it to the seed.
  if(!seed.addInfo(1801,dTr.info(1801,-9999)))Warning("Could not add info 1801 to DS track").ignore();
  if(!seed.addInfo(1802,dTr.info(1802,-9999)))Warning("Could not add info 1802 to DS track").ignore();
  if(!seed.addInfo(1803,dTr.info(1803,-9999)))Warning("Could not add info 1803 to DS track").ignore();
  if(!seed.addInfo(1804,dTr.info(1804,-9999)))Warning("Could not add info 1804 to DS track").ignore();
  if(!seed.addInfo(1805,dTr.info(1805,-9999)))Warning("Could not add info 1805 to DS track").ignore();
  if(!seed.addInfo(1806,dTr.info(1806,-9999)))Warning("Could not add info 1806 to DS track").ignore();
  if(!seed.addInfo(1807,dTr.info(1807,-9999)))Warning("Could not add info 1807 to DS track").ignore();

  if(!seed.addInfo(1900,dTr.info(1900,-9999)))Warning("Could not add info 1901 to DS track").ignore();
  if(!seed.addInfo(1901,dTr.info(1901,-9999)))Warning("Could not add info 1901 to DS track").ignore();
  if(!seed.addInfo(1902,dTr.info(1902,-9999)))Warning("Could not add info 1902 to DS track").ignore();
  if(!seed.addInfo(1903,dTr.info(1903,-9999)))Warning("Could not add info 1903 to DS track").ignore();
  if(!seed.addInfo(1904,dTr.info(1904,-9999)))Warning("Could not add info 1904 to DS track").ignore();
  if(!seed.addInfo(1905,dTr.info(1905,-9999)))Warning("Could not add info 1905 to DS track").ignore();
  if(!seed.addInfo(1906,dTr.info(1906,-9999)))Warning("Could not add info 1906 to DS track").ignore();
  if(!seed.addInfo(1907,dTr.info(1907,-9999)))Warning("Could not add info 1907 to DS track").ignore();
  std::cout<<"writing extra info to seed. Test info(1801) = "<<seed.info(1801,-999)<<std::endl;
  // Now that all of that info is there, we can write the rest of the extra information.
  return StatusCode::SUCCESS;
}
