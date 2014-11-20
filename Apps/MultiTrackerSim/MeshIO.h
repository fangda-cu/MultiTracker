//
//  MeshIO.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MeshIO__
#define __MeshIO__

#include <iostream>
#include <vector>
#include "surftrack.h"

class MeshIO
{
public:
    static bool save(LosTopos::SurfTrack & st, const std::string & filename, bool binary = true);
    static bool load(LosTopos::SurfTrack & st, const std::string & filename, bool binary = true);
    
    static bool loadIntoRaw(std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls, const std::string & filename, bool binary = true);
    
};



#endif /* defined(__MeshIO__) */
