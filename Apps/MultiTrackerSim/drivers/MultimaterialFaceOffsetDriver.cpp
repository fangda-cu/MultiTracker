//
//  MultimaterialFaceOffsetDriver.h
//
//  adapted from code by Tyson Brochu 2008, with extensions for multimaterial
//  Christopher Batty, Fang Da 2014
//
//  Mesh driver for motion in the normal direction using the faceoff method (entropy solution).
//

#include "MultimaterialFaceOffsetDriver.h"

#include <cassert>
#include <lapack_wrapper.h>
#include <meshsmoother.h>
#include <nondestructivetrimesh.h>
#include <runstats.h>
#include <surftrack.h>

using namespace LosTopos;

MultimaterialFaceOffsetDriver * MultimaterialFaceOffsetDriver::s_singleton = new MultimaterialFaceOffsetDriver;

bool MultimaterialFaceOffsetDriver::step(LosTopos::SurfTrack * st, double dt, std::vector<std::vector<double> > & speeds, int expanding_surf, bool nm_stationary, bool bbwall)
{
    assert(!bbwall);

    s_singleton->speed_matrix = speeds;
    s_singleton->expanding_surface = expanding_surf;
    s_singleton->nonmanifold_stationary = nm_stationary;
    
    std::vector<LosTopos::Vec3d> new_positions = st->pm_positions;
    double curr_dt = dt;
    s_singleton->set_predicted_vertex_positions(*st, new_positions, dt, curr_dt);
    st->set_all_newpositions(new_positions);
    
    return true;
}

// Get the quadric metric tensor at a vertex from the given incident triangles
void MultimaterialFaceOffsetDriver::compute_quadric_metric_tensor(const std::vector<Vec3d> & triangle_normals,
                                                                  const std::vector<double> & triangle_areas,
                                                                  const std::vector<size_t> & incident_triangles,
                                                                  Mat33d & quadric_metric_tensor)
{
    std::vector< Vec3d > N(incident_triangles.size());
    std::vector< double > W(incident_triangles.size());
    
    for ( size_t i = 0; i < incident_triangles.size(); ++i )
    {
        size_t triangle_index = incident_triangles[i];
        N[i] = triangle_normals[triangle_index];
        W[i] = triangle_areas[triangle_index];
        
        assert( triangle_areas[triangle_index] > 0.0 );
    }
    
    zero( quadric_metric_tensor );
    
    // Ax = b from N^TWNx = N^TWd
    for ( size_t i = 0; i < N.size(); ++i )
    {
        quadric_metric_tensor(0,0) += N[i][0] * W[i] * N[i][0];
        quadric_metric_tensor(1,0) += N[i][1] * W[i] * N[i][0];
        quadric_metric_tensor(2,0) += N[i][2] * W[i] * N[i][0];
        
        quadric_metric_tensor(0,1) += N[i][0] * W[i] * N[i][1];
        quadric_metric_tensor(1,1) += N[i][1] * W[i] * N[i][1];
        quadric_metric_tensor(2,1) += N[i][2] * W[i] * N[i][1];
        
        quadric_metric_tensor(0,2) += N[i][0] * W[i] * N[i][2];
        quadric_metric_tensor(1,2) += N[i][1] * W[i] * N[i][2];
        quadric_metric_tensor(2,2) += N[i][2] * W[i] * N[i][2];
    }
}

// Return intersection point between a set of planes in the least-squares sense
void MultimaterialFaceOffsetDriver::intersection_point(const std::vector<Vec3d> & triangle_normals,
                                                       const std::vector<double> & triangle_plane_distances,
                                                       const std::vector<double> & triangle_areas,
                                                       const std::vector<Vec2i> & triangles_labels,
                                                       const std::vector<size_t> & incident_triangles,
                                                       Vec3d & out)
{
    
    std::vector< Vec3d > N;
    std::vector< double > W;
    std::vector< double > d;
    std::vector<size_t> reduced_tris;
    
    reduced_tris.reserve(incident_triangles.size());
    N.reserve(incident_triangles.size());
    W.reserve(incident_triangles.size());
    d.reserve(incident_triangles.size());

    for ( size_t i = 0; i < incident_triangles.size(); ++i )
    {
        size_t triangle_index = incident_triangles[i];
        
        Vec2i label = triangles_labels[i]; //don't process if it's not a growing face
        if(expanding_surface != -1) {
           if(label[0] != expanding_surface && label[1] != expanding_surface)
              continue;
        }
        
        reduced_tris.push_back(triangle_index);
        N.push_back( triangle_normals[triangle_index] );
        W.push_back( triangle_areas[triangle_index] );
        d.push_back( triangle_plane_distances[triangle_index] );
    }
    
    Mat33d A(0,0,0,0,0,0,0,0,0);
    Vec3d b(0,0,0);
    
    compute_quadric_metric_tensor( triangle_normals, triangle_areas, reduced_tris, A );
    
    for ( size_t i = 0; i < N.size(); ++i )
    {
        
        b[0] += N[i][0] * W[i] * d[i];
        b[1] += N[i][1] * W[i] * d[i];
        b[2] += N[i][2] * W[i] * d[i];      
    }
    
    double eigenvalues[3];
    double work[9];
    int info = ~0, n = 3, lwork = 9;
    LAPACK::get_eigen_decomposition( &n, A.a, &n, eigenvalues, work, &lwork, &info );
    
    if ( info != 0 )
    {
        assert(0);
    }
    
    Vec3d displacement_vector(0,0,0);
    
    for ( unsigned int i = 0; i < 3; ++i )
    {
        if ( eigenvalues[i] > G_EIGENVALUE_RANK_RATIO * eigenvalues[2] )
        {
            Vec3d eigenvector( A(0,i), A(1,i), A(2,i) );
            displacement_vector += dot(eigenvector, b) * eigenvector / eigenvalues[i];
        }
    }
    
    out = displacement_vector;
}


// Assign a velocity vector to each vertex in the input mesh, and adjust the time step size
void MultimaterialFaceOffsetDriver::set_predicted_vertex_positions(const SurfTrack & surf,
                                                                   std::vector<Vec3d> & new_positions,
                                                                   double current_t,
                                                                   double & adaptive_dt)
{
    const NonDestructiveTriMesh& mesh = surf.m_mesh;
    std::vector<double> triangle_areas;
    triangle_areas.reserve(mesh.num_triangles());
    std::vector<Vec3d> triangle_normals;
    triangle_normals.reserve(mesh.num_triangles());
    std::vector<Vec3d> triangle_centroids;
    triangle_centroids.reserve(mesh.num_triangles());
    std::vector<double> triangle_plane_distances;
    triangle_plane_distances.reserve(mesh.num_triangles());
    std::vector<Vec2i> triangle_labels;
    triangle_labels.reserve(mesh.num_triangles());

    const std::vector<Vec3st>& tris = mesh.get_triangles();
    
    for ( size_t i = 0; i < tris.size(); ++i )
    {
        Vec2i label;
        if ( tris[i][0] == tris[i][1] )
        {
            triangle_areas.push_back( 0 );
            triangle_normals.push_back( Vec3d(0,0,0) );
            triangle_centroids.push_back( Vec3d(0,0,0) );
            label = Vec2i(-1,-1);
            triangle_labels.push_back(label);
        }
        else
        {
            triangle_areas.push_back( surf.get_triangle_area( i ) );
            triangle_normals.push_back( surf.get_triangle_normal( i ) );
            triangle_centroids.push_back( (surf.get_position(tris[i][0]) + surf.get_position(tris[i][1]) + surf.get_position(tris[i][2])) / 3 );
            label = surf.m_mesh.get_triangle_label(i);

            triangle_labels.push_back( label );
        }
        
        double switch_speed = 0;
        if(label[0] >= 0 && label[1] >= 0)
           switch_speed = speed_matrix[label[0]][label[1]];

        triangle_plane_distances.push_back( adaptive_dt * switch_speed );
    }
    
    std::vector<Vec3d> displacements;
    displacements.resize( surf.get_num_vertices(), Vec3d(0,0,0) );
    
    //
    // Null space smoothing - to get tangential component
    //
    
    {
        for ( size_t i = 0; i < surf.get_num_vertices(); ++i )
        {
            const MeshSmoother& smoother = surf.m_smoother;
            smoother.null_space_smooth_vertex( i, triangle_areas, triangle_normals, triangle_centroids, displacements[i] );
        }
    }
    

    //
    // Primary space displacement - to get normal component
    //
    //#pragma omp parallel for
    for ( int p = 0; p < (int)surf.get_num_vertices(); ++p )
    {

       if(nonmanifold_stationary) {
          //if it's a non-manifold vertex, and we're holding those stationary
          //don't move it. e.g. for the curling sphere example.
          std::set<int> labelset;
          for(size_t i = 0; i < mesh.m_vertex_to_triangle_map[p].size(); ++i) {
             size_t tri = mesh.m_vertex_to_triangle_map[p][i];
             Vec2i label = mesh.get_triangle_label(tri);
             labelset.insert(label[0]);
             labelset.insert(label[1]);
          }
          if(labelset.size() > 2) continue;
       }
        Vec3d normal_displacement;
        intersection_point( triangle_normals, triangle_plane_distances, triangle_areas, triangle_labels, mesh.m_vertex_to_triangle_map[p], normal_displacement );
        
        //
        // Entropy solution
        //
        
        if(mag(normal_displacement) <= 1e-8) //the surface is not moving
        {
           continue;
        }

        if ( surf.m_mesh.m_vertex_to_triangle_map[p].empty() )
        {
            continue;
        }
        
        double sum_mu_l = 0, sum_mu = 0;
        
        const std::vector<size_t>& incident_triangles = mesh.m_vertex_to_triangle_map[p];
        std::vector<size_t> manifold_tris;
        for(size_t j = 0; j < incident_triangles.size(); ++j) {
           Vec2i labels = triangle_labels[incident_triangles[j]];

           if(expanding_surface != -1) {
              if(labels[0] != expanding_surface && labels[1] != expanding_surface) 
                 continue;
           }
           manifold_tris.push_back(incident_triangles[j]);
        }

        for ( size_t j = 0; j < manifold_tris.size(); ++j )
        {
            size_t triangle_index = manifold_tris[j];
            
            const Vec3st& tri = surf.m_mesh.get_triangle( triangle_index );
            
            Vec3d edge_vector;
            if ( tri[0] == p )
            {
                edge_vector = surf.get_position(tri[1]) - surf.get_position(tri[2]);
            }
            else if ( tri[1] == p )
            {
                edge_vector = surf.get_position(tri[2]) - surf.get_position(tri[0]);
            }
            else
            {
                edge_vector = surf.get_position(tri[0]) - surf.get_position(tri[1]);
            }
            
            Vec3d tri_normal = triangle_normals[triangle_index];

            Vec3d s = cross( tri_normal, edge_vector );   // orthogonal to normal and edge opposite vertex

            Vec3d total_displacement_estimate = displacements[p] + normal_displacement;
            bool contracting = dot( s, total_displacement_estimate) >= 0.0;
            
            double cos_theta = dot( tri_normal, normal_displacement ) / mag(normal_displacement);
            
            assert(!std::isnan(cos_theta));

            double mu = triangle_areas[triangle_index];
            
            double li = fabs( triangle_plane_distances[triangle_index] ); 
            
            if ( contracting )
            {
                li /= fabs( cos_theta );
            }
            
            sum_mu_l += mu * li;
            sum_mu += mu;
        }
        
        double length = sum_mu_l / sum_mu;

        //add the normal displacement to the tangential one we already have
        displacements[p] += length * normal_displacement / mag(normal_displacement);

    }
    
    //compute timestep that will avoid self intersection and mesh folding
    double beta = MeshSmoother::compute_max_timestep_quadratic_solve( surf.m_mesh.get_triangles(), surf.get_positions(), displacements, false );
    
    adaptive_dt *= beta;

    for(size_t i = 0; i < surf.get_num_vertices(); i++)
    {
        new_positions[i] = surf.get_position(i) + beta * displacements[i];
        assert(!std::isnan(new_positions[i][0]) && !std::isnan(new_positions[i][1]) && !std::isnan(new_positions[i][2]));
    }
}

