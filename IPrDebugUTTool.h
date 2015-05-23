// $Id: IPrDebugUTTool.h,v 1.1 2008-12-04 09:07:35 cattanem Exp $
#ifndef PATKERNEL_IPATDEBUGUTTOOL_H 
#define PATKERNEL_IPATDEBUGUTTOOL_H 1

// Include files
#include "PrKernel/PrUTHit.h"

// Forward declarations
class MsgStream;

namespace LHCb{
  class Track;
}

// from Gaudi
#include "GaudiKernel/IAlgTool.h"
//std
#include <map>
#include <vector>
static const InterfaceID IID_IPrDebugUTTool ( "IPrDebugUTTool", 1, 0 );

/** @class IPrDebugUTTool IPrDebugUTTool.h PrKernel/IPrDebugUTTool.h
 *  
 *
 *  @author Olivier Callot
 *  @date   2007-10-22
 */
class IPrDebugUTTool : virtual public IAlgTool {
public: 

  // Return the interface ID
  static const InterfaceID& interfaceID() { return IID_IPrDebugUTTool; }

  virtual void debugUTClusterOnTrack( const LHCb::Track* track, 
                                      const PrUTHits::const_iterator beginCoord,
                                      const PrUTHits::const_iterator endCoord   ) = 0;
  
  virtual void debugUTCluster( MsgStream& msg, const PrUTHit* hit ) = 0;


  //added by AD 4/28/15, for use with mika's tool.
    
  virtual void recordStepInProcess(std::string step, bool result) = 0;
  

  virtual StatusCode writeExtraInfoToDownstreamTrack(LHCb::Track& dsTrack) = 0;

  virtual StatusCode writeExtraInfoToSeed(LHCb::Track& dTr, LHCb::Track& seed) = 0;

  virtual bool isIDOnMCParticle(LHCb::LHCbID id, LHCb::Track& track) = 0;//ad
  
  virtual bool CheckMCinIntermediateHits(std::vector<int> containerx,
                                         std::vector<int> containeruv) = 0;

  virtual bool isTrackReconstructible(LHCb::Track& track) = 0;
  
  virtual std::vector<double>TrueHitXYZ(LHCb::LHCbID id,bool extended = false) = 0;
  
  virtual bool CheckMCinTrackHits(std::vector<int>ids) = 0;

  virtual void printTSeedInfo(LHCb::Track* tr) = 0;
  virtual void printTSeedInfoShort(LHCb::Track* tr) = 0;

  virtual bool isPreselGood(std::vector<LHCb::LHCbID> x1ids,
                            std::vector<LHCb::LHCbID> uids,
                            std::vector<LHCb::LHCbID> vids,
                            std::vector<LHCb::LHCbID> x2ids,
                            LHCb::Track * track) = 0;
  virtual void PrintHitTable(PrUTHit* hit) = 0;
  virtual void PrintHitTableShort(PrUTHit* hit) = 0;


protected:

private:

};
#endif // PATKERNEL_IPATDEBUGUTTOOL_H
