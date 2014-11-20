//
//  MTRenderer.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __MTRenderer__
#define __MTRenderer__

#include "surftrack.h"

class MTRenderer
{
public:
    MTRenderer();
    
public:
    void render();
    
    void setMousePos(double x, double y) { m_mouse_x = x; m_mouse_y = y; }
    
    void setSurfTrack(LosTopos::SurfTrack * st, int nregions) { m_st = st; m_nregions = nregions; m_region_visible = std::vector<bool>(nregions, true); }
    
    enum Mode
    {
        RM_INTERFACES,
        RM_JUNCTIONS,
        RM_REGIONS,
        RM_OPAQUE,
        
        RM_COUNT
    };
    
    void cycleMode(int inc = 1);
    
    void keyboard(unsigned char key, int x, int y);

public:
    int nearestVertex() const { return m_nearest_vertex; }
    int nearestEdge() const   { return m_nearest_edge; }
    int nearestFace() const   { return m_nearest_face; }
    
protected:
    bool visibleV(size_t v) const;
    bool visibleE(size_t e) const;
    bool visibleF(size_t f) const;
    
protected:
    LosTopos::SurfTrack * m_st;
    
    int m_nregions;
    int m_current_region;
    std::vector<bool> m_region_visible;
    
    Mode m_mode;
    
    double m_mouse_x;
    double m_mouse_y;
    
    int m_nearest_vertex;
    int m_nearest_edge;
    int m_nearest_face;
    
};

#endif
