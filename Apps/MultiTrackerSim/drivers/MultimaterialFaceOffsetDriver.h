//
//  MultimaterialFaceOffsetDriver.h
//
//  adapted from code by Tyson Brochu 2008, with extensions for multimaterial
//  Christopher Batty, Fang Da 2014
//
//  Mesh driver for motion in the normal direction using the faceoff method (entropy solution).
//

#ifndef FACEOFF_MULTI_H
#define FACEOFF_MULTI_H

#include <mat.h>
#include <vec.h> 
#include <vector>

#include "surftrack.h"

// Motion in the normal direction using Face Offsetting [Jiao 2007].
class MultimaterialFaceOffsetDriver
{
public:
    static bool step(LosTopos::SurfTrack * st, double dt, std::vector<std::vector<double> > & speeds, int expanding_surf, bool nm_stationary, bool bbwall = false);
    
public:
    MultimaterialFaceOffsetDriver() : speed_matrix(), expanding_surface(-1), nonmanifold_stationary(false) { }

    // Get the quadric metric tensor at a vertex from the given incident triangles
    void compute_quadric_metric_tensor(const std::vector<LosTopos::Vec3d>& triangle_normals, const std::vector<double>& triangle_areas, const std::vector<size_t>& incident_triangles, LosTopos::Mat33d & quadric_metric_tensor);
    
    // Return intersection point between a set of planes in the least-squares sense
    void intersection_point(const std::vector<LosTopos::Vec3d>& triangle_normals, const std::vector<double>& triangle_plane_distances, const std::vector<double>& triangle_areas, const std::vector<LosTopos::Vec2i>& triangle_labels, const std::vector<size_t>& incident_triangles, LosTopos::Vec3d& out);
    
    // Assign a velocity vector to each mesh vertex
    void set_predicted_vertex_positions(const LosTopos::SurfTrack & surf, std::vector<LosTopos::Vec3d> & predicted_positions, double current_t, double & adaptive_dt);
    
    // Speed of normal motion
    std::vector<std::vector<double> > speed_matrix; //dictate speeds based on pairwise labels
    
    // (outer) surface region to use for offsetting and null-space smoothing (i.e. for handling non-manifold cases)
    int expanding_surface;

    // turn on if the non-manifold vertices are not to move (e.g. in curling sphere example)
    bool nonmanifold_stationary;

protected:
    static MultimaterialFaceOffsetDriver * s_singleton;
    
};

#endif


