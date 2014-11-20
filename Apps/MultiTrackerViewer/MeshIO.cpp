//
//  MeshIO.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include <fstream>
#include "MeshIO.h"

bool MeshIO::save(LosTopos::SurfTrack & st, const std::string & filename, bool binary)
{
    if (binary)
    {
        std::ofstream os(filename.c_str(), std::ios::binary);
        size_t n;
        
        n = st.m_mesh.nv();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x = st.get_position(i);
            os.write((char *)&(x[0]), sizeof (x[0]));
            os.write((char *)&(x[1]), sizeof (x[1]));
            os.write((char *)&(x[2]), sizeof (x[2]));
        }
        
        n = st.m_mesh.nt();
        os.write((char *)&n, sizeof (size_t));
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os.write((char *)&(t[0]), sizeof (t[0]));
            os.write((char *)&(t[1]), sizeof (t[1]));
            os.write((char *)&(t[2]), sizeof (t[2]));
            
            os.write((char *)&(l[0]), sizeof (l[0]));
            os.write((char *)&(l[1]), sizeof (l[1]));
        }
        
        return os.good();
    } else
    {
        std::ofstream os(filename.c_str());
        os << st.m_mesh.nv() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nv(); i++)
        {
            os << st.get_position(i) << std::endl;
        }
        
        os << st.m_mesh.nt() << std::endl;
        for (size_t i = 0; i < st.m_mesh.nt(); i++)
        {
            LosTopos::Vec3st t = st.m_mesh.get_triangle(i);
            if (t[0] == t[1])
                continue;
            LosTopos::Vec2i l = st.m_mesh.get_triangle_label(i);
            
            os << t << " " << l << std::endl;
        }
        
        os.close();
        
        return os.good();
    }
}

bool MeshIO::load(LosTopos::SurfTrack & st, const std::string & filename, bool binary)
{
    std::ifstream test(filename.c_str());
    if (!test.is_open())
    {
        std::cout << "[MeshIO::load] Error: file " << filename << " not found." << std::endl;
        return false;
    }
    
    for (size_t i = 0; i < st.m_mesh.nt(); i++)
    {
        if (st.m_mesh.get_triangle(i)[0] == st.m_mesh.get_triangle(i)[1])
            continue;
        st.remove_triangle(i);
    }
    
    for (size_t i = 0; i < st.m_mesh.nv(); i++)
        st.remove_vertex(i);
    
    if (binary)
    {
        std::ifstream is(filename.c_str(), std::ios::binary);
        
        size_t n;
        n = st.m_mesh.nv();
        is.read((char *)&n, sizeof (size_t));
        
        st.m_mesh.set_num_vertices(n);
        std::vector<LosTopos::Vec3d> pos(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is.read((char *)&(x[0]), sizeof (x[0]));
            is.read((char *)&(x[1]), sizeof (x[1]));
            is.read((char *)&(x[2]), sizeof (x[2]));
            pos[i] = x;
        }
        
        st.m_masses.resize(n);
        for (size_t i = 0; i < n; i++)
            st.m_masses[i] = LosTopos::Vec3d(1, 1, 1);
        
        st.set_all_positions(pos);
        st.set_all_newpositions(pos);
        st.set_all_remesh_velocities(std::vector<LosTopos::Vec3d>(n, LosTopos::Vec3d(0)));
        
        n = st.m_mesh.nt();
        is.read((char *)&n, sizeof (size_t));
        
        std::vector<LosTopos::Vec3st> tris;
        std::vector<LosTopos::Vec2i> labels;
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is.read((char *)&(t[0]), sizeof (t[0]));
            is.read((char *)&(t[1]), sizeof (t[1]));
            is.read((char *)&(t[2]), sizeof (t[2]));
            tris.push_back(t);
            
            LosTopos::Vec2i l;
            is.read((char *)&(l[0]), sizeof (l[0]));
            is.read((char *)&(l[1]), sizeof (l[1]));
            labels.push_back(l);
        }
        
        st.m_mesh.replace_all_triangles(tris, labels);
        
        size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
        st.pm_positions.resize(nv);
        st.pm_newpositions.resize(nv);
        st.pm_velocities.resize(nv);
        st.m_velocities.resize(nv);
        
        is.close();
        
        return is.good();
    } else
    {
        std::ifstream is(filename.c_str());
        
        size_t n;
        is >> n;
        
        st.m_mesh.set_num_vertices(n);
        std::vector<LosTopos::Vec3d> pos(n);
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3d x;
            is >> x[0] >> x[1] >> x[2];
            pos[i] = x;
        }
        
        st.m_masses.resize(n);
        for (size_t i = 0; i < n; i++)
            st.m_masses[i] = LosTopos::Vec3d(1, 1, 1);
        
        st.set_all_positions(pos);
        st.set_all_newpositions(pos);
        st.set_all_remesh_velocities(std::vector<LosTopos::Vec3d>(n, LosTopos::Vec3d(0)));
        
        n = st.m_mesh.nt();
        is >> n;
        
        std::vector<LosTopos::Vec3st> tris;
        std::vector<LosTopos::Vec2i> labels;
        for (size_t i = 0; i < n; i++)
        {
            LosTopos::Vec3st t;
            is >> t[0] >> t[1] >> t[2];
            tris.push_back(t);
            
            LosTopos::Vec2i l;
            is >> l[0] >> l[1];
            labels.push_back(l);
        }
        
        st.m_mesh.replace_all_triangles(tris, labels);
        
        size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
        st.pm_positions.resize(nv);
        st.pm_newpositions.resize(nv);
        st.pm_velocities.resize(nv);
        st.m_velocities.resize(nv);
        
        is.close();
        
        return is.good();
    }
}


