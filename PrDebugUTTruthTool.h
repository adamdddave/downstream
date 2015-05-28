// $Id: PrDebugUTTruthTool.h,v 1.4 2008-12-04 09:05:07 cattanem Exp $
#ifndef PRDEBUGUTTRUTHTOOL_H 
#define PRDEBUGUTTRUTHTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "PrKernel/IPrDebugUTTool.h"            // Interface
#include "PrKernel/PrUTHit.h"//pr ut hits, to see if it helps
#include "Event/MCParticle.h"
#include "Linker/LinkerTool.h"
#include "TrackInterfaces/ITrackExtrapolator.h"

#include <map>
#include <vector>
class DeSTDetector;
class ILHCbIDsToMCParticles;
class ILHCbIDsToMCHits;
//struct for PrDebugUTTruth
struct TSeed_Analyzer 
{
  //PrDownTrack initial_candidate;
  int tot_candidates;
  int presel_count;
  PrUTHits presel_hits_x;
  PrUTHits presel_hits_uv;
  std::vector<bool> mc_candidate;        
};// added 1-29-14. AD  

/** @class PrDebugUTTruthTool PrDebugUTTruthTool.h
 *  
 *
 *  @author Olivier Callot
 *  @date   2007-10-22
 */
class PrDebugUTTruthTool : public GaudiTool, virtual public IPrDebugUTTool {
public: 
  /// Standard constructor
  PrDebugUTTruthTool( const std::string& type, 
                       const std::string& name,
                       const IInterface* parent);

  virtual ~PrDebugUTTruthTool( ); ///< Destructor

  virtual StatusCode initialize(); /// initialize

  virtual void debugUTClusterOnTrack( const LHCb::Track* track, 
                                      const PrUTHits::const_iterator beginCoord,
                                      const PrUTHits::const_iterator endCoord   );
  
  virtual void debugUTCluster( MsgStream& msg, const PrUTHit* hit );

  //added by AD 4/15/15, for use with mika's tool.
  
  virtual void recordStepInProcess(std::string step, bool result);
  
  //virtual void write_info_to_TSeed(LHCb::Track& seed, std::string OutputLocation);
  
  virtual StatusCode writeExtraInfoToDownstreamTrack(LHCb::Track& dsTrack);

  virtual StatusCode writeExtraInfoToSeed(LHCb::Track& dTr, LHCb::Track& seed);

  virtual bool isIDOnMCParticle(LHCb::LHCbID id, LHCb::Track& track);//ad
  
  virtual bool CheckMCinIntermediateHits(std::vector<int> containerx,
                                         std::vector<int> containeruv);

  virtual bool isTrackReconstructible(LHCb::Track& track);
  
  virtual std::vector<double>TrueHitXYZ(LHCb::LHCbID id,bool extended = false);
  
  virtual bool CheckMCinTrackHits(std::vector<int>ids);
  virtual void printTSeedInfo(LHCb::Track* tr);
  virtual void printTSeedInfoShort(LHCb::Track* tr);
  virtual bool isPreselGood(std::vector<LHCb::LHCbID> x1ids,
                            std::vector<LHCb::LHCbID> uids,
                            std::vector<LHCb::LHCbID> vids,
                            std::vector<LHCb::LHCbID> x2ids,
                            LHCb::Track * track);
  virtual void PrintHitTable(PrUTHit* hit);
  virtual void PrintHitTableShort(const PrUTHit* hit);
  
  
protected:

private:

  DeSTDetector* m_tracker;
  std::map<std::string,bool> m_flags;
  std::map<int,bool>m_x1_matching_map, m_x2_matching_map;//maps to tell whether the LHCbID
  std::map<int,bool>m_u_matching_map, m_v_matching_map;// is matched to the tseed. save some calc time.
  
  typedef LinkerTool<LHCb::Track, LHCb::MCParticle> MyAsct;
  typedef MyAsct::DirectType            Table;
  typedef MyAsct::DirectType::Range     Range;
  typedef Table::iterator               iterator;

  typedef MyAsct::InverseType           InvTable;
  typedef InvTable::Range               InvRange;
  typedef InvTable::iterator            InvIterator;

  MyAsct*         m_link;
  const InvTable* m_invTable;

  ILHCbIDsToMCParticles * m_lhcbid2mcparticles;
  ILHCbIDsToMCHits * m_lhcbid2mchits;

};
#endif // PRDEBUGUTTRUTHTOOL_H
