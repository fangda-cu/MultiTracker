//
//  MCFDriver.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include "MCFDriver.h"

using namespace LosTopos;

bool MCFDriver::step(SurfTrack * st, double dt, bool bbwall)
{
    std::vector<Vec3d> v;
    
    // adaptive time stepping
    double dt_left = dt;
    double dt_largest_possible = 0;
    
    int substep_count = 0;
    while (dt_left > 0 && substep_count < 100)
    {
        evaluateV(st, v);
        
        // determine max dt
        dt_largest_possible = determineMaxDt(st, v);
        std::cout << "dt_left = " << dt_left << " max dt = " << dt_largest_possible << std::endl;
        if (dt_largest_possible >= dt_left)
            dt_largest_possible = dt_left;

        for (size_t i = 0; i < st->m_mesh.nv(); i++)
        {
            Vec3d m = st->m_masses[i];
            Vec3d newpos = st->get_position(i) + Vec3d(v[i][0] / m[0], v[i][1] / m[1], v[i][2] / m[2]) * dt_largest_possible;
            assert(newpos == newpos);
            if (bbwall) newpos = Vec3d(std::min(1.0, std::max(0.0, newpos[0])), std::min(1.0, std::max(0.0, newpos[1])), std::min(1.0, std::max(0.0, newpos[2])));
            st->set_newposition(i, newpos);
        }
        
        dt_left -= dt_largest_possible;
        substep_count++;
    }
    
    if (dt_left > 0)
    {
        std::cout << "[MCFDriver::step] Warning: time step incomplete after maximum number of substeps." << std::endl;
        return false;
    } else
    {
        std::cout << "[MCFDriver::step] Notice: time step finished after " << substep_count << " substeps." << std::endl;
        return true;
    }
}

double MCFDriver::determineMaxDt(SurfTrack * st, std::vector<Vec3d> & v)
{
    double global_max_dt = 1.0;
    for (size_t i = 0; i < st->m_mesh.ne(); i++)
    {
        size_t v0 = st->m_mesh.m_edges[i][0];
        size_t v1 = st->m_mesh.m_edges[i][1];
        Vec3d edge = st->get_position(v1) - st->get_position(v0);
        Vec3d rv = v[v1] - v[v0];
        double approaching_velocity = -dot(rv, edge) / dot(edge, edge);
        
        double max_dt = 0;
        if (approaching_velocity <= 0)
            max_dt = std::numeric_limits<double>::infinity();
        else
            max_dt = 1 / approaching_velocity * 0.8;    // allow the edge to be shrunk by 80% in one time step at most
        
        if (max_dt < global_max_dt)
            global_max_dt = max_dt;
    }
    
    return global_max_dt;
}

void MCFDriver::evaluateV(LosTopos::SurfTrack * st, std::vector<LosTopos::Vec3d> & v)
{
    const double COEF = 1;
    
    v.resize(st->m_mesh.nv(), Vec3d(0, 0, 0));
    for (size_t i = 0; i < st->m_mesh.nt(); i++)
    {
        Vec3d p0 = st->get_position(st->m_mesh.m_tris[i][0]);
        Vec3d p1 = st->get_position(st->m_mesh.m_tris[i][1]);
        Vec3d p2 = st->get_position(st->m_mesh.m_tris[i][2]);
        
        Vec3d v01 = p1 - p0;
        Vec3d v20 = p0 - p2;
        
        Vec3d A = cross(v01, -v20);
        double Anorm = mag(A);
        Vec3d mul = A / Anorm;
        
        Vec3d p2part = COEF * cross(v01, mul);
        Vec3d p1part = COEF * cross(mul, -v20);
        Vec3d p0part = -(p1part + p2part);

        v[st->m_mesh.m_tris[i][0]] += p0part;
        v[st->m_mesh.m_tris[i][1]] += p1part;
        v[st->m_mesh.m_tris[i][2]] += p2part;
    }
}
