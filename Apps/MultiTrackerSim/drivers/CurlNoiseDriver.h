//
//  CurlNoiseDriver.h
//
//  adapted from code by Tyson Brochu 2008
//  Christopher Batty, Fang Da 2014
//
//

#ifndef CURLNOISEDRIVER_H
#define CURLNOISEDRIVER_H

#include "curlnoise/curlnoise.h"
#include "curlnoise/noise.h"

#include "surftrack.h"

class CurlNoiseDriver: public CurlNoise3
{
public:
    static bool step(LosTopos::SurfTrack * st, double dt, bool bbwall = false);
    
public:
    CurlNoiseDriver(void);
    
    virtual bool seed_particles(std::vector<LosTopos::Vec3d>& /*x*/, double /*dt*/ ) const { return true; }
    void advance_time(double /*dt*/) {}

    // z-component of noise vector (x and y are 0)
    double noise2(double x, double y, double z) const { return noise(z-203.994, x+169.47, y-205.31); }
    
    // Take the curl of this function to get velocity
    LosTopos::Vec3d potential(double x, double y, double z) const;
    
    // Set velocity on each mesh vertex
    void set_predicted_vertex_positions(const LosTopos::SurfTrack & surf, std::vector<LosTopos::Vec3d> & predicted_positions, double current_t, double & adaptive_dt);
    
    // Parameters for each scale of noise
    std::vector<double> noise_lengthscale, noise_gain;
    
    // Noise generator
    FlowNoise3 noise;
    
protected:
    static CurlNoiseDriver * s_singleton;
    
};

#endif
