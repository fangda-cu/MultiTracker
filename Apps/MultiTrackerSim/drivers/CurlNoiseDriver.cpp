//
//  CurlNoiseDriver.h
//
//  adapted from code by Tyson Brochu 2008
//  Christopher Batty, Fang Da 2014
//
//

#include "CurlNoiseDriver.h"

using namespace LosTopos;

CurlNoiseDriver * CurlNoiseDriver::s_singleton = new CurlNoiseDriver;

bool CurlNoiseDriver::step(LosTopos::SurfTrack * st, double dt, bool bbwall)
{
    assert(!bbwall);

    std::vector<LosTopos::Vec3d> new_positions = st->pm_positions;
    double curr_dt = dt;
    s_singleton->set_predicted_vertex_positions(*st, new_positions, dt, curr_dt);
    st->set_all_newpositions(new_positions);

    return true;
}

CurlNoiseDriver::CurlNoiseDriver() :
   CurlNoise3(),
   noise_lengthscale(1),
   noise_gain(1),
   noise()
{
    noise_lengthscale[0]=1.5;
    noise_gain[0]=1.3;  
}

// 3D vector field which defines the velocity via its curl.
Vec3d CurlNoiseDriver::potential(double x, double y, double z) const
{
    Vec3d psi(0,0,0);
    double height_factor=0.5;
    
    static const Vec3d centre( 0.0, 1.0, 0.0 );
    static double radius = 4.0;
    
    for(unsigned int i=0; i<noise_lengthscale.size(); ++i)
    {
        double sx=x/noise_lengthscale[i];
        double sy=y/noise_lengthscale[i];
        double sz=z/noise_lengthscale[i];
        
        Vec3d psi_i( 0.f, 0.f, noise2(sx,sy,sz));
        
        double dist = mag( Vec3d(x,y,z) - centre );      
        double scale = max( (radius - dist)/radius, 0.0 );
        psi_i *= scale;
        
        psi+=height_factor*noise_gain[i]*psi_i;
    }
    
    return psi;
}

// For each vertex on the mesh, take the curl of potential to get velocity vector at that vertex location.
// Uses RK4 to get x(t+1) from x(t), then returns linear trajectory v = (x(t+1) - x(t)) / dt.
void CurlNoiseDriver::set_predicted_vertex_positions(const SurfTrack & surf, std::vector<Vec3d> & predicted_positions, double /*current_t*/, double & adaptive_dt)
{
    const std::vector<Vec3d>& positions = surf.get_positions();
    Vec3d v, midx;
    
    // Update velocities
    for ( unsigned int i = 0; i < positions.size(); ++i )
    {    
        
        const Vec3d& p = positions[i];
        
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() )
        {
            predicted_positions[i] = p;
            continue;
        }
        
        // RK4
        // -----------
        // k1 = dt * f( t, x );
        get_velocity( p, v );
        Vec3d k1 = adaptive_dt * v;
        
        // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
        advance_time( 0.5*adaptive_dt );
        get_velocity( p + 0.5 * k1, v );
        Vec3d k2 = adaptive_dt * v;
        
        //k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
        get_velocity( p + 0.5 * k2, v );
        Vec3d k3 = adaptive_dt * v;
        
        //k4 = dt * f( t + dt, x + k3 );
        advance_time( 0.5*adaptive_dt );
        get_velocity( p + 0.5 * k3, v );
        Vec3d k4 = adaptive_dt * v;
        
        predicted_positions[i] = p + 1.0/6.0 * ( k1 + k4 ) + 1.0/3.0 * ( k2 + k3 );
        
        advance_time( -adaptive_dt );
    }
    
    advance_time( adaptive_dt );
}


