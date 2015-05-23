// $Id: PrDownTrack2.h,v 1.3 2009-10-08 10:09:46 sstahl Exp $
#ifndef PRDOWNTRACK2_H 
#define PRDOWNTRACK2_H 1

// Include files
#include "GaudiKernel/Point3DTypes.h"

#include "Event/Track.h"
#include "Event/State.h"

#include "PrKernel/PrUTHit.h"
#include "TfKernel/RecoFuncs.h"

/** @class PrDownTrack2 PrDownTrack2.h
 *  Track helper for Downstream track search
 *  Adapted from PrDownTrack to accomodate the new chi2 calculation better
 *  
 *  There are two parts to this helper. First, all the information which pertains
 *  to the T Seed. This has m_seed_ as the prefix.
 *
 *  The second part only pertains to the UT. This only has m_ as the prefix.
 * 
 *  @author Adam Davis
 *  @date   2015-05-18
 */

class PrDownTrack2 {
public:
  /// Standard constructor
  PrDownTrack2( LHCb::Track* tr, //the seed
		//the following are commented out for now, but will be changed later
                const std::vector<double>& magnetParams,
                const std::vector<double>& momentumParams,
                //const std::vector<double>& yParams,
                //double errZMag,
                double magnetScale );
  virtual ~PrDownTrack2( ) {} ///< Destructor
  LHCb::Track* track() const {return m_seed_track;}//tseed things
  LHCb::State* state() const {return m_seed_state;}//same
  

  //getters
  double xMagnet() const {return m_seed_magnet.x();}//these look the same as before
  double yMagnet() const {return m_seed_magnet.y();}//but are drastically different due
  double zMagnet() const {return m_seed_magnet.z();}//to straight line projections.
  double slopeX() const {return m_slopeX;}
  double slopeY() const {return m_slopeY;}
  PrUTHits& hits(){return m_hits;}
  const PrUTHits & hits() const {return m_hits;}
  double xAtZ(double z) const { return m_x0 + m_slopeX*z; }//just straight lines
  double yAtZ(double z) const { return m_y0 + m_slopeY*z; }
  //old momentum for now
  double moment()  const {return ( (*m_momPar_old)[0] +    
                                        (*m_momPar_old)[1] * m_seed_state->tx() * m_seed_state->tx() +
                                        (*m_momPar_old)[2] * m_seed_state->ty() * m_seed_state->ty() ) / 
      ( m_seed_state->tx() - m_slopeX ) * m_magnetScale;
  }
  double chisq(){return m_chi2;}  
  double distance( const PrUTHit* hit ) const {return hit->x() - xAtZ( hit->z() );}

  //setters
  void setXat0(double x0){m_x0 = x0;}
  void setYat0(double y0){m_y0 = y0;}
  void setSlopeX(double slx){m_slopeX = slx;}
  void setSlopeY(double sly){m_slopeY = sly;}
  void updateTrack(double dx,double dslx, double dy, double dsly){
    m_x0+=dx;//after the fit to the candidate, we will adjust the guesses by the difference
    m_y0+=dy;//between the guess and the truth. That's the reason that we use the dx,etc
    m_slopeX+=dslx;//here
    m_slopeY+=dsly;
  }
  void setChisq(double chi2){m_chi2 = chi2;}
  
  //member function
  void startNewCandidate(PrUTHit* firstHit){
    m_hits.clear();
    m_hits.push_back(firstHit);
  }
  void startNewCandidate(){
    m_hits.clear();
    m_x0 = 0;
    m_y0 = 0;
  }

  void sortFinalHits() {
    std::sort( m_hits.begin(), m_hits.end(), Tf::increasingByZ<PrUTHit>() );
  }

  
  
    
private:
  LHCb::Track* m_seed_track;
  LHCb::State* m_seed_state;
  Gaudi::XYZPoint m_seed_magnet;
  double m_slopeX;
  double m_slopeY;

  double m_x0;
  double m_y0;
  PrUTHits m_hits;
  double m_magnetScale;
  const std::vector<double>* m_magPar_old;
  const std::vector<double>* m_momPar_old;
  const std::vector<double>* m_magPar;
  double m_chi2;
  
};

#endif // PRDOWNTRACK2_H
