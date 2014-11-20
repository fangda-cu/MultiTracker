//
//  MTRenderer.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include "MTRenderer.h"
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <vector>
#include <algorithm>
#include <Eigen/Core>
#include <Eigen/Dense>

typedef Eigen::Matrix<double, 4, 4> Mat4d;
typedef Eigen::Matrix<double, 3, 3> Mat3d;
typedef Eigen::Matrix<double, 4, 1> Vec4d;
typedef Eigen::Matrix<double, 3, 1> Vec3d;
typedef Eigen::Matrix<double, 2, 1> Vec2d;

MTRenderer::MTRenderer() :
    m_st(NULL),
    m_nregions(0),
    m_current_region(0),
    m_mode(RM_REGIONS),
    m_mouse_x(0),
    m_mouse_y(0)
{
    
}

namespace
{
    //Color routine to use HSV colors instead of RGB taken
    // from http://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
    // HSV \in (0-360, 0-1, 0-1)
    void hsvToRgb(double _h, double _s, double _v, double & r, double & g, double & b)
    {
        double c, m, x, hp;
        double h = std::min(std::max(_h, 0.), 360.);
        double s = std::min(std::max(_s, 0.), 1.);
        double v = std::min(std::max(_v, 0.), 1.);
        c = v * s;
        m = v - c;
        hp = h / 60.;
        x = c * (double)(1 - abs((int)hp % 2 - 1));
        assert (hp < 6. && hp >= 0);
        if (hp < 1)
        {
            r = c;
            g = x;
            b = 0;
        } else if (hp < 2)
        {
            r = x;
            g = c;
            b = 0;
        } else if (hp < 3)
        {
            r = 0;
            g = c;
            b = x;
        } else if (hp < 4 )
        {
            r = 0;
            g = x;
            b = c;
        } else if (hp < 5)
        {
            r = x;
            g = 0;
            b = c;
        } else if (hp < 6)
        {
            r = c;
            g = 0;
            b = x;
        }
        r += m;
        g += m;
        c += m;
    }

    void glColorHSV3d(double h, double s, double v)
    {
        double r, g, b;
        hsvToRgb(h, s, v, r, g, b);
        glColor3d(r, g, b);
    }

    void glVertexVec3d(const Vec3d& v)
    {
        glVertex3d((GLdouble)v[0], (GLdouble)v[1], (GLdouble)v[2]);
    }

    void glNormalVec3d(const Vec3d& v)
    {
        glNormal3d((GLdouble)v[0], (GLdouble)v[1], (GLdouble)v[2]);
    }

    class FaceComp
    {
    public:
        bool operator() (const std::pair<size_t, double> & f1, const std::pair<size_t, double> & f2) const
        {
            return f1.second < f2.second;
        }
    };

    int junction(const LosTopos::NonDestructiveTriMesh & mesh, size_t e)
    {
        return mesh.m_edge_to_triangle_map[e].size();
    }

    bool junctionNeighborV(const LosTopos::NonDestructiveTriMesh & mesh, size_t v)
    {
        for (size_t i = 0; i < mesh.m_vertex_to_edge_map[v].size(); i++)
            if (junction(mesh, mesh.m_vertex_to_edge_map[v][i]) > 2)
                return true;
        
        return false;
    }

    bool junctionNeighborE(const LosTopos::NonDestructiveTriMesh & mesh, size_t e)
    {
        if (junctionNeighborV(mesh, mesh.m_edges[e][0]))
            return true;
        
        if (junctionNeighborV(mesh, mesh.m_edges[e][1]))
            return true;
        
        return false;
    }

    bool junctionNeighborF(const LosTopos::NonDestructiveTriMesh & mesh, size_t f)
    {
        for (size_t i = 0; i < 3; i++)
            if (junctionNeighborV(mesh, mesh.m_tris[f][i]))
                return true;
        
        return false;
    }

    int onBBWall(const Vec3d & pos)
    {
        static const double WALL_THRESHOLD = 1e-6;
        
        int walls = 0;
        if (pos[0] < 0 + WALL_THRESHOLD) walls |= (1 << 0);
        if (pos[1] < 0 + WALL_THRESHOLD) walls |= (1 << 1);
        if (pos[2] < 0 + WALL_THRESHOLD) walls |= (1 << 2);
        if (pos[0] > 1 - WALL_THRESHOLD) walls |= (1 << 3);
        if (pos[1] > 1 - WALL_THRESHOLD) walls |= (1 << 4);
        if (pos[2] > 1 - WALL_THRESHOLD) walls |= (1 << 5);
        
        return walls;
    }

    LosTopos::Vec3d EVec(const Vec3d & v) { return LosTopos::Vec3d(v[0], v[1], v[2]); }
    Vec3d          Vec(const LosTopos::Vec3d & v) { return Vec3d(v[0], v[1], v[2]); }
 
}

bool MTRenderer::visibleV(size_t v) const
{
    bool visible = false;
    for (size_t i = 0; i < m_st->m_mesh.m_vertex_to_triangle_map[v].size(); i++)
    {
        const LosTopos::Vec2i & l = m_st->m_mesh.get_triangle_label(m_st->m_mesh.m_vertex_to_triangle_map[v][i]);
        if (l[0] >= 0 && m_region_visible[l[0]]) visible = true;
        if (l[1] >= 0 && m_region_visible[l[1]]) visible = true;
    }
    
    return visible;
}

bool MTRenderer::visibleE(size_t e) const
{
    bool visible = false;
    for (size_t i = 0; i < m_st->m_mesh.m_edge_to_triangle_map[e].size(); i++)
    {
        const LosTopos::Vec2i & l = m_st->m_mesh.get_triangle_label(m_st->m_mesh.m_edge_to_triangle_map[e][i]);
        if (l[0] >= 0 && m_region_visible[l[0]]) visible = true;
        if (l[1] >= 0 && m_region_visible[l[1]]) visible = true;
    }
    
    return visible;
}

bool MTRenderer::visibleF(size_t f) const
{
    bool visible = false;
    const LosTopos::Vec2i & l = m_st->m_mesh.get_triangle_label(f);
    if (l[0] >= 0 && m_region_visible[l[0]]) visible = true;
    if (l[1] >= 0 && m_region_visible[l[1]]) visible = true;

    return visible;
}

void MTRenderer::cycleMode(int inc)
{
    m_mode = (MTRenderer::Mode)((m_mode + RM_COUNT + inc) % RM_COUNT);
}

void MTRenderer::keyboard(unsigned char key, int x, int y)
{
    if (key == 'm' || key == 'M')
    {
        cycleMode(key == 'm' ? 1 : -1);
        glutPostRedisplay();
    } else if (key == '=' || key == '+')
    {
        m_current_region = (m_current_region + (key == '+' ? 10 : 1)) % m_nregions;
        glutPostRedisplay();
        std::cout << "current region = " << m_current_region << "/" << m_nregions << std::endl;
    } else if (key == '-' || key == '_')
    {
        m_current_region = (m_current_region + m_nregions - (key == '_' ? 10 : 1)) % m_nregions;
        glutPostRedisplay();
        std::cout << "current region = " << m_current_region << "/" << m_nregions << std::endl;
    } else if (key == 'l' || key == 'L')
    {
        while ((int)m_region_visible.size() < m_nregions)
            m_region_visible.push_back(true);
        m_region_visible[m_current_region] = !m_region_visible[m_current_region];
        glutPostRedisplay();
        for (size_t i = 0; i < m_region_visible.size(); i++)
            std::cout << "region " << i << " " << (m_region_visible[i] ? "visible" : "invisible") << std::endl;
    } else if (key == 'h' || key == 'H')
    {
        while ((int)m_region_visible.size() < m_nregions)
            m_region_visible.push_back(true);
        bool allinvisible = true;
        for (size_t i = 0; i < m_region_visible.size(); i++)
            if (m_region_visible[i])
                allinvisible = false;
        for (size_t i = 0; i < m_region_visible.size(); i++)
            m_region_visible[i] = allinvisible;
        glutPostRedisplay();
    }
    
}

void MTRenderer::render()
{
    assert(m_st);
    
    // find the vertex/edge/face the mouse cursor is nearest to
    Vec2d mousepos = Vec2d(m_mouse_x, m_mouse_y);
    
    Mat4d MV;
    Mat4d PJ;
    {
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
        float pj[16];
        glGetFloatv(GL_PROJECTION_MATRIX, pj);
        MV << mv[0], mv[4], mv[8], mv[12], mv[1], mv[5], mv[9], mv[13], mv[2], mv[6], mv[10], mv[14], mv[3], mv[7], mv[11], mv[15];
        PJ << pj[0], pj[4], pj[8], pj[12], pj[1], pj[5], pj[9], pj[13], pj[2], pj[6], pj[10], pj[14], pj[3], pj[7], pj[11], pj[15];
    }
    Mat4d MVP = PJ * MV;
    
    double mind = -1;
    int mind_vertex = -1;
    int mind_edge = -1;
    int mind_face = -1;
    
    for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
    {
        if (!visibleV(i))
            continue;
        
        Vec3d pos = Vec(m_st->get_position(i));
        Vec4d scrpos_h = MVP * Vec4d(pos.x(), pos.y(), pos.z(), 1.0);
        Vec2d scrpos = Vec2d(scrpos_h.x(), scrpos_h.y()) / scrpos_h.w();
        
        double distance = (scrpos - mousepos).norm();
        if (distance < mind || mind < 0)
        {
            mind = distance;
            mind_vertex = i;
        }
    }
    
    for (size_t i = 0; i < m_st->m_mesh.ne(); i++)
    {
        if (!visibleE(i))
            continue;
        
        Vec3d v0 = Vec(m_st->get_position(m_st->m_mesh.m_edges[i][0]));
        Vec3d v1 = Vec(m_st->get_position(m_st->m_mesh.m_edges[i][1]));
        
        Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
        Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
        Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
        Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
        
        double distance = (mousepos - (scrv0 + scrv1) / 2).norm();
//        double distance = (mousepos - scrv0 - (mousepos - scrv0).dot(scrv1 - scrv0) * (scrv1 - scrv0) / (scrv1 - scrv0).squaredNorm()).norm();
        if (distance < mind || mind < 0)
        {
            mind = distance;
            mind_vertex = -1;
            mind_edge = i;
        }
    }
    
    for (size_t i = 0; i < m_st->m_mesh.nt(); i++)
    {
        if (!visibleF(i))
            continue;
        
        Vec3d v0 = Vec(m_st->get_position(m_st->m_mesh.m_tris[i][0]));
        Vec3d v1 = Vec(m_st->get_position(m_st->m_mesh.m_tris[i][1]));
        Vec3d v2 = Vec(m_st->get_position(m_st->m_mesh.m_tris[i][2]));
        
        Vec4d scrv0_h = MVP * Vec4d(v0.x(), v0.y(), v0.z(), 1.0);
        Vec2d scrv0 = Vec2d(scrv0_h.x(), scrv0_h.y()) / scrv0_h.w();
        Vec4d scrv1_h = MVP * Vec4d(v1.x(), v1.y(), v1.z(), 1.0);
        Vec2d scrv1 = Vec2d(scrv1_h.x(), scrv1_h.y()) / scrv1_h.w();
        Vec4d scrv2_h = MVP * Vec4d(v2.x(), v2.y(), v2.z(), 1.0);
        Vec2d scrv2 = Vec2d(scrv2_h.x(), scrv2_h.y()) / scrv2_h.w();
        
        double distance = (mousepos - (scrv0 + scrv1 + scrv2) / 3).norm();
        if (distance < mind || mind < 0)
        {
            mind = distance;
            mind_vertex = -1;
            mind_edge = -1;
            mind_face = i;
        }
    }
    
    if (mind_vertex < 0 && mind_edge < 0 && mind_face < 0)  // this can only happen when no region is visible, thus no primitives to pick from
    {
        mind_vertex = 0;
        mind = 0;
    }
    
    m_nearest_vertex = mind_vertex;
    m_nearest_edge = mind_edge;
    m_nearest_face = mind_face;
    
    glPushMatrix();
    
    const LosTopos::NonDestructiveTriMesh & mesh = m_st->m_mesh;
    
    if (m_mode == RM_INTERFACES || m_mode == RM_JUNCTIONS)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        std::vector<Vec3d> faceNormals(m_st->m_mesh.nt());
        for (size_t i = 0; i < m_st->m_mesh.nt(); i++) faceNormals[i] = Vec(m_st->get_triangle_normal(i));
        
        // stats on total number of labels
        int maxlabel = -1;
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            LosTopos::Vec2i regions = mesh.get_triangle_label(i);
            if (regions[0] >= maxlabel) maxlabel = regions[0];
            if (regions[1] >= maxlabel) maxlabel = regions[1];
        }
        m_nregions = maxlabel + 1;
        while (m_region_visible.size() < m_nregions)
            m_region_visible.push_back(true);
        
        // render bounding box
        glLineWidth(1);
        glColor3f(0.5, 0.5, 0.5);
        glBegin(GL_LINES);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 1.0);
        
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 0.0, 1.0);
        glEnd();
        
        // Render all edges
        glLineWidth(2);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBegin(GL_LINES);
        glColor3f(0, 0, 0);
        
        for (size_t i = 0; i < mesh.ne(); i++)
        {
            if (m_mode == RM_JUNCTIONS)
                if (!junctionNeighborE(mesh, i))
                    continue;
            
            Vec3d p0 = Vec(m_st->get_position(mesh.m_edges[i][0]));
            Vec3d p1 = Vec(m_st->get_position(mesh.m_edges[i][1]));
            //Vec3d dir = (p1-p0);
            //p0 = p0 + 0.05*dir;
            //p1 = p1 - 0.05*dir;
            glColor4f(0.0, 0.0, 0.0, 0.2);
            
            if (m_mode == RM_JUNCTIONS)
            {
                int ne = junction(mesh, i);
                if (ne == 3)
                    glColor4f(1.0, 0.0, 1.0, 0.2);
                else if (ne == 4)
                    glColor4f(0.3, 0.8, 0.9, 0.2);
                else if (ne > 4)
                    glColor4f(0.2, 0.3, 1.0, 0.2);
            }
            
            if ((p0 - p1).norm() <= 0.1)
                glColor4f(0.0, 0.5, 1.0, 1.0);
            
            if (!visibleE(i))
                glColor4f(0.0, 0.0, 0.0, 0.1);
            
            glVertexVec3d(p0);
            glVertexVec3d(p1);
        }
        glEnd();
        
        if (mind_edge >= 0)
        {
            glColor4f(0.0, 0.0, 0.0, 1.0);
            glLineWidth(4);
            glBegin(GL_LINES);
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][0])));
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][1])));
            glEnd();
            glLineWidth(2);
        }
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
        
        // render face labels
        if (true)
        {
            // generate a list of colors for all the labels present
            std::vector<Vec3d> labelcolors;
            labelcolors.push_back(Vec3d(0, 0, 0));
            for (int i = 0; i < maxlabel + 1; i++)
            {
                float r, g, b;
                float t = (float)i / (maxlabel == 0 ? 1 : maxlabel);
                if (t < 0.5)
                {
                    r = 1 - t * 1.5f;
                    g = t * 1.5f;
                    b = 0;
                } else
                {
                    r = 0;
                    g = 1 - (t - 0.5f) * 1.5f;
                    b = (t - 0.5f) * 1.5f;
                }
                labelcolors.push_back(Vec3d(r, g, b));
            }
            
            glBegin(GL_LINES);
            for (size_t i = 0; i < mesh.nt(); i++)
            {
                if (m_mode == RM_JUNCTIONS)
                    if (!junctionNeighborF(mesh, i))
                        continue;
                
                LosTopos::Vec2i regions = mesh.get_triangle_label(i);
                if ((regions[0] < 0 || !m_region_visible[regions[0]]) &&
                    (regions[1] < 0 || !m_region_visible[regions[1]]))
                    continue;
                
                Vec3d p0 = Vec(m_st->get_position(mesh.m_tris[i][0]));
                Vec3d p1 = Vec(m_st->get_position(mesh.m_tris[i][1]));
                Vec3d p2 = Vec(m_st->get_position(mesh.m_tris[i][2]));
                
                Vec3d c = (p0 + p1 + p2) / 3;
                Vec3d n = (p1 - p0).cross(p2 - p0).normalized();
                
                double mean_edge_length = ((p0 - p1).norm() + (p1 - p2).norm() + (p2 - p0).norm()) / 3;
                
                Vec3d color0 = labelcolors[std::max(0, regions[0] + 1)];
                glColor3d(color0.x(), color0.y(), color0.z());
                
                if (regions[0] >= 0)
                {
                    glVertexVec3d(c);
                    glVertexVec3d(c - n * mean_edge_length * 0.05);
                }
                
                Vec3d color1 = labelcolors[std::max(0, regions[1] + 1)];
                glColor3d(color1.x(), color1.y(), color1.z());
                
                if (regions[1] >= 0)
                {
                    glVertexVec3d(c);
                    glVertexVec3d(c + n * mean_edge_length * 0.05);
                }
            }
            glEnd();
        }
        
        // render back to front semi transparent
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
        Vec3d view_vec(mv[2], mv[6], mv[10]);  // assuming ModelView matrix contains only translation, rotation and uniform scaling
        
        std::vector<std::pair<size_t, double> > sorted_faces;
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            if (m_mode == RM_JUNCTIONS)
                if (!junctionNeighborF(mesh, i))
                    continue;
            
            Vec3d barycentre(0, 0, 0);
            for (int j = 0; j < 3; j++)
                barycentre += Vec(m_st->get_position(mesh.m_tris[i][j]));
            barycentre /= 3.0;
            
            double depth = barycentre.dot(view_vec);
            sorted_faces.push_back(std::pair<size_t, double>(i, depth));
        }
        
        FaceComp fc;
        std::sort(sorted_faces.begin(), sorted_faces.end(), fc);
        
        // Render all faces
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBegin(GL_TRIANGLES);
        
        for (size_t i = 0; i < sorted_faces.size(); i++)
        {
            Vec3d barycentre(0, 0, 0);
            for (int j = 0; j < 3; j++)
                barycentre += Vec(m_st->get_position(mesh.m_tris[sorted_faces[i].first][j]));
            barycentre /= 3.0;
            
            double alpha = 0.1;
            glColor4f(1.0, 0.0, 0.0, alpha);
            
            if (!visibleF(sorted_faces[i].first))
                glColor4f(0.0, 0.0, 0.0, 0.02);
            
            double edge_shrink = 0.1;
            
            if (sorted_faces[i].first == mind_face)
            {
                glColor4f(1.0, 0.0, 0.0, 0.5);
                edge_shrink = 0.0;
            }
            
            for (int j = 0; j < 3; j++)
            {
                Vec3d pos = Vec(m_st->get_position(mesh.m_tris[sorted_faces[i].first][j]));
                pos += (barycentre - pos) * edge_shrink;
                glVertexVec3d(pos);
            }
        }
        
        glEnd();
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        
        glColor3f(0.0, 0.0, 0.0);
        
        //Render all vertices
        glPointSize(4);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBegin(GL_POINTS);
        
        for (size_t i = 0; i < mesh.nv(); i++)
        {
            Vec3d vertPos = Vec(m_st->get_position(i));
            
            double minangle = -1;
            double maxangle = -1;
            for (size_t j = 0; j < mesh.m_vertex_to_triangle_map[i].size(); j++)
            {
                size_t v0 = mesh.m_tris[mesh.m_vertex_to_triangle_map[i][j]][0];
                size_t v1 = mesh.m_tris[mesh.m_vertex_to_triangle_map[i][j]][1];
                size_t v2 = mesh.m_tris[mesh.m_vertex_to_triangle_map[i][j]][2];
                
                double angle = 0;
                Vec3d a, b;
                if (v0 == i)
                {
                    a = Vec(m_st->get_position(v1));
                    b = Vec(m_st->get_position(v2));
                } else if (v1 == i)
                {
                    a = Vec(m_st->get_position(v2));
                    b = Vec(m_st->get_position(v0));
                } else
                {
                    a = Vec(m_st->get_position(v0));
                    b = Vec(m_st->get_position(v1));
                }
                angle = acos(((a - vertPos).squaredNorm() + (b - vertPos).squaredNorm() - (a - b).squaredNorm()) / (2 * (a - vertPos).norm() * (b - vertPos).norm()));
                if (minangle < 0 || angle < minangle)
                    minangle = angle;
                if (maxangle < 0 || angle > maxangle)
                    maxangle = angle;
            }
            
            double minedge = -1;
            for (size_t j = 0; j < mesh.m_vertex_to_edge_map[i].size(); j++)
            {
                double edgelength = m_st->get_edge_length(mesh.m_vertex_to_edge_map[i][j]);
                if (minedge < 0 || edgelength < minedge)
                    minedge = edgelength;
            }
            
            if (m_mode == RM_JUNCTIONS)
                if (!junctionNeighborV(mesh, i))
                    continue;
            
            if (minedge < 0.01)
                glColor4f(1.0, 1.0, 0.0, 1.0);
            else if (minangle < 3 * M_PI / 180)
                glColor4f(0.0, 1.0, 1.0, 1.0);
            else if (maxangle > 177 * M_PI / 180)
                glColor4f(1.0, 0.0, 1.0, 1.0);
            else if(onBBWall(vertPos))
                glColor4f(0.0, 1.0, 0.0, 1.0);
            else
                glColor4f(0.0, 0.0, 0.0, 0.1);
            
            glVertexVec3d(vertPos);
        }
        glEnd();
        
        if (mind_vertex >= 0)
        {
            glPointSize(8);
            glColor4f(0, 0, 0, 1);
            glBegin(GL_POINTS);
            glVertexVec3d(Vec(m_st->get_position(mind_vertex)));
            glEnd();
        }
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        
    } else if (m_mode == RM_REGIONS)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        std::vector<Vec3d> faceNormals(m_st->m_mesh.nt());
        for (size_t i = 0; i < m_st->m_mesh.nt(); i++) faceNormals[i] = Vec(m_st->get_triangle_normal(i));
        
        // stats on total number of labels
        int maxlabel = -1;
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            LosTopos::Vec2i regions = mesh.get_triangle_label(i);
            if (regions[0] >= maxlabel) maxlabel = regions[0];
            if (regions[1] >= maxlabel) maxlabel = regions[1];
        }
        m_nregions = maxlabel + 1;
        while (m_region_visible.size() < m_nregions)
            m_region_visible.push_back(true);
        
        // render bounding box
        glLineWidth(1);
        glColor3f(0.5, 0.5, 0.5);
        glBegin(GL_LINES);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 1.0);
        
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 0.0, 1.0);
        glEnd();
        
        // Render all edges
        glLineWidth(1);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBegin(GL_LINES);

        for (size_t i = 0; i < mesh.ne(); i++)
        {
            Vec3d p0 = Vec(m_st->get_position(mesh.m_edges[i][0]));
            Vec3d p1 = Vec(m_st->get_position(mesh.m_edges[i][1]));
            
            bool visible = visibleE(i);
            
            bool xb = ((p0.x() < 1e-4 && p1.x() < 1e-4) || (p0.x() > 1 - 1e-4 && p1.x() > 1 - 1e-4));
            bool yb = ((p0.y() < 1e-4 && p1.y() < 1e-4) || (p0.y() > 1 - 1e-4 && p1.y() > 1 - 1e-4));
            bool zb = ((p0.z() < 1e-4 && p1.z() < 1e-4) || (p0.z() > 1 - 1e-4 && p1.z() > 1 - 1e-4));
            
            bool cubeedge = ((xb && yb) || (xb && zb) || (yb && zb));
            
            if (!visible && !cubeedge)
                continue;
            
            glColor4f(0.0f, 0.0f, 0.0f, 0.3f);
            
            glVertexVec3d(p0);
            glVertexVec3d(p1);
        }
        glEnd();
        
        if (mind_edge >= 0)
        {
            glColor4f(0.0, 0.0, 0.0, 1.0);
            glLineWidth(2);
            glBegin(GL_LINES);
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][0])));
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][1])));
            glEnd();
            glLineWidth(1);
        }
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        
        // generate a list of colors for all the labels present
        std::vector<Vec3d> labelcolors;
        labelcolors.push_back(Vec3d(0, 0, 0));
        for (int i = 0; i < maxlabel + 1; i++)
        {
            float r, g, b;
            float t = (float)i / (maxlabel == 0 ? 1 : maxlabel);
            if (t < 0.5)
            {
                r = 1 - t * 1.5f;
                g = t * 1.5f;
                b = 0;
            } else
            {
                r = 0;
                g = 1 - (t - 0.5f) * 1.5f;
                b = (t - 0.5f) * 1.5f;
            }
            labelcolors.push_back(Vec3d(r, g, b));
        }
        
        // render back to front semi transparent
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
        Vec3d view_vec(mv[2], mv[6], mv[10]);  // assuming ModelView matrix contains only translation, rotation and uniform scaling
        
        std::vector<std::pair<size_t, double> > sorted_faces;
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            Vec3d barycentre(0, 0, 0);
            for (int j = 0; j < 3; j++)
                barycentre += Vec(m_st->get_position(mesh.m_tris[i][j]));
            barycentre /= 3.0;
            
            double depth = barycentre.dot(view_vec);
            sorted_faces.push_back(std::pair<size_t, double>(i, depth));
        }
        
        FaceComp fc;
        std::sort(sorted_faces.begin(), sorted_faces.end(), fc);
        
        // Render all faces
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glBegin(GL_TRIANGLES);
        
        for (size_t i = 0; i < sorted_faces.size(); i++)
        {
            Vec3d barycentre(0, 0, 0);
            for (int j = 0; j < 3; j++)
                barycentre += Vec(m_st->get_position(mesh.m_tris[sorted_faces[i].first][j]));
            barycentre /= 3.0;
            
            LosTopos::Vec2i regions = mesh.get_triangle_label(sorted_faces[i].first);
            
            Vec3d color0 = labelcolors[regions[0] + 1];
            Vec3d color1 = labelcolors[regions[1] + 1];
            Vec3d color_combined = Vec3d::Zero();
            double alpha = 0.02;
            int visible_count = 0;
            if (regions[0] >= 0 && m_region_visible[regions[0]])
            {
                color_combined += color0;
                visible_count++;
                alpha += 0.3;
            }
            if (regions[1] >= 0 && m_region_visible[regions[1]])
            {
                color_combined += color1;
                visible_count++;
                alpha += 0.3;
            }
            if (visible_count == 0)
                color_combined = Vec3d(0.0,0.0,0.0);
            else
                color_combined /= visible_count;
            
            glColor4f(color_combined.x(), color_combined.y(), color_combined.z(), alpha);
            
            if (sorted_faces[i].first == mind_face)
                glColor4f(color_combined.x(), color_combined.y(), color_combined.z(), 1.0);
            
            for (int j = 0; j < 3; j++)
            {
                Vec3d pos = Vec(m_st->get_position(mesh.m_tris[sorted_faces[i].first][j]));
                glVertexVec3d(pos);
            }
        }
        
        glEnd();
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        
        //Render vertex of interest
        if (mind_vertex >= 0)
        {
            glPointSize(4);
            glColor4f(0, 0, 0, 1);
            glBegin(GL_POINTS);
            glVertexVec3d(Vec(m_st->get_position(mind_vertex)));
            glEnd();
        }
        
    } else if (m_mode == RM_OPAQUE)
    {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        
        std::vector<Vec3d> faceNormals(m_st->m_mesh.nt());
        for (size_t i = 0; i < m_st->m_mesh.nt(); i++) faceNormals[i] = Vec(m_st->get_triangle_normal(i));
        
        // stats on total number of labels
        int maxlabel = -1;
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            LosTopos::Vec2i regions = mesh.get_triangle_label(i);
            if (regions[0] >= maxlabel) maxlabel = regions[0];
            if (regions[1] >= maxlabel) maxlabel = regions[1];
        }
        m_nregions = maxlabel + 1;
        while (m_region_visible.size() < m_nregions)
            m_region_visible.push_back(true);
        
        // render bounding box
        glLineWidth(1);
        glColor3f(0.5, 0.5, 0.5);
        glBegin(GL_LINES);
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 1.0);
        
        glVertex3f(0.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 0.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(1.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 1.0, 0.0);
        glVertex3f(0.0, 0.0, 0.0);
        
        glVertex3f(0.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 0.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(1.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 1.0, 1.0);
        glVertex3f(0.0, 0.0, 1.0);
        glEnd();
        
        // Render all non-manifold junctions
        glLineWidth(1);
        glDisable(GL_BLEND);
        glDepthMask(GL_TRUE);
        glPolygonOffset(1, 2);
        glEnable(GL_POLYGON_OFFSET_FILL);
        glBegin(GL_LINES);
        
        for (size_t i = 0; i < mesh.ne(); i++)
        {
            if (junction(mesh, i) == 2)
                continue;
            
            Vec3d p0 = Vec(m_st->get_position(mesh.m_edges[i][0]));
            Vec3d p1 = Vec(m_st->get_position(mesh.m_edges[i][1]));
            
            bool visible = visibleE(i);
            
            bool xb = ((p0.x() < 1e-4 && p1.x() < 1e-4) || (p0.x() > 1 - 1e-4 && p1.x() > 1 - 1e-4));
            bool yb = ((p0.y() < 1e-4 && p1.y() < 1e-4) || (p0.y() > 1 - 1e-4 && p1.y() > 1 - 1e-4));
            bool zb = ((p0.z() < 1e-4 && p1.z() < 1e-4) || (p0.z() > 1 - 1e-4 && p1.z() > 1 - 1e-4));
            
            bool cubeedge = ((xb && yb) || (xb && zb) || (yb && zb));
            
            if (!visible && !cubeedge)
                continue;
            
            glColor4f(0.0f, 0.0f, 0.0f, 0.3f);
            
            glVertexVec3d(p0);
            glVertexVec3d(p1);
        }
        glEnd();
        
        if (mind_edge >= 0)
        {
            glColor4f(0.0, 0.0, 0.0, 1.0);
            glLineWidth(2);
            glBegin(GL_LINES);
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][0])));
            glVertexVec3d(Vec(m_st->get_position(mesh.m_edges[mind_edge][1])));
            glEnd();
            glLineWidth(1);
        }
        
        // generate a list of colors for all the labels present
        std::vector<Vec3d> labelcolors;
        labelcolors.push_back(Vec3d(0, 0, 0));
        for (int i = 0; i < maxlabel + 1; i++)
        {
            float r, g, b;
            float t = (float)i / (maxlabel == 0 ? 1 : maxlabel);
            if (t < 0.5)
            {
                r = 1 - t * 1.5f;
                g = t * 1.5f;
                b = 0;
            } else
            {
                r = 0;
                g = 1 - (t - 0.5f) * 1.5f;
                b = (t - 0.5f) * 1.5f;
            }
            labelcolors.push_back(Vec3d(r, g, b));
        }
        
        // render all faces opaque
        float mv[16];
        glGetFloatv(GL_MODELVIEW_MATRIX, mv);
//        Vec3d view_vec(mv[2], mv[6], mv[10]);  // assuming ModelView matrix contains only translation, rotation and uniform scaling
        // since the face between two regions would switch between two distinct colors when it flips from facing the viewer with one side to another, simply
        // using the z axis in the model view matrix does not suffice in the perspective projection mode.
        Vec3d eye = (MV.inverse() * Vec4d(0, 0, 0, 1)).segment<3>(0);
        
        glBegin(GL_TRIANGLES);
        
        for (size_t i = 0; i < mesh.nt(); i++)
        {
            const LosTopos::Vec3st & f = mesh.m_tris[i];
            
            Vec3d x0 = Vec(m_st->get_position(f[0]));
            Vec3d x1 = Vec(m_st->get_position(f[1]));
            Vec3d x2 = Vec(m_st->get_position(f[2]));
            
            LosTopos::Vec2i regions = mesh.get_triangle_label(i);
            
            Vec3d color0 = labelcolors[regions[0] + 1];
            Vec3d color1 = labelcolors[regions[1] + 1];
            Vec3d color_final = Vec3d::Zero();
            if ((regions[0] >= 0 && m_region_visible[regions[0]]) && (regions[1] >= 0 && m_region_visible[regions[1]]))
            {
                color_final = ((x1 - x0).cross(x2 - x0).dot(eye - x0) < 0 ? color0 : color1);
            } else if (regions[0] >= 0 && m_region_visible[regions[0]])
            {
                color_final = color0;
            } else if (regions[1] >= 0 && m_region_visible[regions[1]])
            {
                color_final = color1;
            } else
            {
                continue; // neither region visible: don't render this face at all
            }
            
            glColor3f(0.4 + 0.6 * color_final.x(), 0.4 + 0.6 * color_final.y(), 0.4 + 0.6 * color_final.z());
            
            if (i == mind_face)
                glColor3f(color_final.x() * 0.6, color_final.y() * 0.6, color_final.z() * 0.6);
            
            glVertexVec3d(x0);
            glVertexVec3d(x1);
            glVertexVec3d(x2);
        }
        
        glEnd();
        glDisable(GL_POLYGON_OFFSET_FILL);

        //Render vertex of interest
        if (mind_vertex >= 0)
        {
            glPointSize(4);
            glColor4f(0, 0, 0, 1);
            glBegin(GL_POINTS);
            glVertexVec3d(Vec(m_st->get_position(mind_vertex)));
            glEnd();
        }
        
    }

    glPopMatrix();
    
}
