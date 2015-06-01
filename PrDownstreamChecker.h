// Include files 

#ifndef PrDownstreamChecker_H
#define PrDownstreamChecker_H
// Include files
#include "Event/Track.h"
#include "Event/MCParticle.h"

// from Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"

// interfaces
#include "MCInterfaces/IMCReconstructible.h"

// linkers
#include "Linker/LinkerTool.h"
#include "Event/MCTrackInfo.h"

/** @class PrDownstreamChecker PrDownstreamChecker.h "PrCheckers/PrDownstreamChecker"
 *
 * Class for upgrade downstream tracking studies
 *  @author M. Vesterinen
 *  @date   16-10-2013
 */

class PrDownstreamChecker : public GaudiTupleAlg
{  
 public:
  /** Standard construtor */
  PrDownstreamChecker( const std::string& name, ISvcLocator* pSvcLocator );
  
  /** Destructor */
  virtual ~PrDownstreamChecker();
  
  /** Algorithm execute */
  virtual StatusCode execute();
  
  /** Algorithm initialize */
  virtual StatusCode initialize();
  
  /** Algorithm finalize */
  virtual StatusCode finalize();
  
protected:
  bool bAncestor(const LHCb::MCParticle* ) const;
  StatusCode loopOverMCParticles();
  StatusCode loopOverTracks();
  StatusCode loopOverSeeds();
  StatusCode loopOverUpdatedSeeds();
  
  void fillMCParticle(const LHCb::MCParticle* p, Tuples::Tuple ntuple);
  void fillTrack(const LHCb::Track* t, Tuples::Tuple ntuple);
  int getUT_layers(const LHCb::MCParticle* p) const;
  int getUT_hits(const LHCb::MCParticle* p) const;

  
private:
  /** Data members */
  std::string m_TSeedContainer;          ///< Input Tracks container location
  std::string m_Test_Container;          ///< Input Tracks container location after modification//ad
  std::string m_DownstreamContainer;          ///< Input Tracks container location
  IMCReconstructible* m_selector;           ///< Pointer to selector
  IMCReconstructible::RecCategory m_recCat; ///<



  typedef LinkerTool<LHCb::Track, LHCb::MCParticle> AsctTool;
  AsctTool m_associator;
  AsctTool::DirectType* m_directTable;
  AsctTool::InverseType* m_inverseTable;
  AsctTool m_associator_seed;
  AsctTool::DirectType* m_directTable_seed;
  AsctTool::InverseType* m_inverseTable_seed;
  //AD 1-10-14
  AsctTool m_associator_updatedseed;
  AsctTool::DirectType* m_directTable_updatedseed;
  AsctTool::InverseType* m_inverseTable_updatedseed;
  
  int nSeedTracks;
  int nDownstreamTracks;
  int nDownstreamTestTracks;//ad
  int evt;
  int n_in_evt;

};

#endif// PrDownstreamChecker_H



