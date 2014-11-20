//
//  MCFDriver.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MCFDriver__
#define __MCFDriver__

#include "surftrack.h"

class MCFDriver
{
public:
    static bool step(LosTopos::SurfTrack * st, double dt, bool bbwall = false);
    
    static double determineMaxDt(LosTopos::SurfTrack * st, std::vector<LosTopos::Vec3d> & v);
    
    static void evaluateV(LosTopos::SurfTrack * st, std::vector<LosTopos::Vec3d> & v);
    
};

#endif
