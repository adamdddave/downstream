// $Id: PrDownTrack2.cpp,v 1.0 2015-05-19 11:07:00 adavis Exp $
// Include files 

// local
#include "PrDownTrack2.h"
#include "PrDownTrack.h"
//-----------------------------------------------------------------------------
// Implementation file for class : PrDownTrack2, from Pat/PatKShort package and
// from PrDownTrack.
//
// 2007-10-18 : Olivier Callot
// 2015-05-19 : Adam Davis
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
PrDownTrack2::PrDownTrack2( LHCb::Track* tr,
                          const std::vector<double>& magnetParams,
                          const std::vector<double>& momentumParams,
                          double magnetScale) : 
  m_magPar_old(&magnetParams), m_momPar_old(&momentumParams)
{
  m_hits.reserve(8);
  m_seed_track = tr;
  m_seed_state = &tr->closestState( 10000. );
  m_slopeY = m_seed_state->ty();
  m_magnetScale = magnetScale;
  m_y0 = m_seed_state->y() - m_slopeY*m_seed_state->z();
  //now hard code a bunch of new magnet parameters right now.
  std::vector<double> params;
  params.push_back(5393);//p0
  params.push_back(-1.363e5);//p1
  params.push_back(1.102e9);//p2
  //so this can be replaced later.
  m_magPar = &params;
  //now use the old momentum tuning to get ther right momentum estimate.
  double zmagnet_old = (*m_magPar_old)[0] + 
    (*m_magPar_old)[1] * m_seed_state->ty() * m_seed_state->ty() +
    (*m_magPar_old)[2] * m_seed_state->tx() * m_seed_state->tx() + (*m_magPar)[3] / m_seed_state->p();
  //this takes the tseed momentum and uses it. We refine using the other momentum.
  double dz = zmagnet_old - m_seed_state->z();
  double xmagnet_old = m_seed_state->x() + dz * m_seed_state->tx();
  m_slopeX = xmagnet_old/zmagnet_old;//starting guess.
  double moment_old = ( (*m_momPar_old)[0] +    
			(*m_momPar_old)[1] * m_seed_state->tx() * m_seed_state->tx() +
			(*m_momPar_old)[2] * m_seed_state->ty() * m_seed_state->ty() ) / 
    ( m_seed_state->tx() - m_slopeX ) * m_magnetScale;
  //now we have the momentum guess from the tseed, which is generally ok. now let's construct the new magnet point.
  double zmagnet_new =(*m_magPar)[0] + (*m_magPar)[1]/fabs(moment_old)+(*m_magPar)[2]/moment_old/moment_old;
  //p0 + p1/momentum + p2/momentum^2
  double xmagnet_new = (m_seed_state->x()) + m_seed_state->tx()*(zmagnet_new - m_seed_state->z());
  double ymagnet_new = (m_seed_state->y()) + m_seed_state->ty()*(zmagnet_new - m_seed_state->z());
  m_seed_magnet = Gaudi::XYZPoint(xmagnet_new, ymagnet_new,zmagnet_new);
}

//alternative constructor
PrDownTrack2::PrDownTrack2(PrDownTrack &tr){
  m_hits = tr.hits();
  m_seed_track = tr.track();
  m_seed_state = tr.state();
  m_slopeY = tr.slopeY();
  m_momPar_old = tr.MomPars();
  m_magPar_old = tr.MagPars();
  m_magnetScale = tr.MagScale();
  m_y0 = tr.yAtZ(0);
  //now hard code a bunch of new magnet parameters right now.
  std::vector<double> params;
  params.push_back(5393);//p0
  params.push_back(-1.363e5);//p1
  params.push_back(1.102e9);//p2
  //so this can be replaced later.
  m_magPar = &params;
  m_slopeX = tr.slopeX();
  
  m_seed_magnet = Gaudi::XYZPoint(tr.xMagnet(), tr.yMagnet(),tr.zMagnet());
  
}
