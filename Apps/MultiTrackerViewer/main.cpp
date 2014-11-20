#include <iostream>
#include <sstream>
#include <iomanip>
#include <map>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#else
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#endif
#include <math.h>
#include "MTRenderer.h"
#include "MeshIO.h"

struct View
{
    static double LDRAG_GAIN;
    static double RDRAG_GAIN;
    
    double alpha;
    double theta;
    double dist;
    
    bool ldrag;
    int x_ldrag;
    int y_ldrag;
    
    bool rdrag;
    int x_rdrag;
    int y_rdrag;

    int win_w;
    int win_h;
    
    MTRenderer renderer;
};

double View::LDRAG_GAIN = 0.005;
double View::RDRAG_GAIN = 0.002;

View g_view;

struct Playback
{
    std::string rec_directory;
    
    int current_frame;
    int current_subframe;
    
    LosTopos::SurfTrack * st;
};

Playback g_playback;

void display()
{
    glClearColor(1, 1, 1, 1);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -g_view.dist);
    glRotatef(g_view.alpha * 180 / M_PI, 1, 0, 0);
    glRotatef(g_view.theta * 180 / M_PI, 0, 1, 0);
    glTranslatef(-0.5, -0.5, -0.5);
    glPushMatrix();
    
//    // coordinate frame
//    glBegin(GL_LINES);
//    glColor3f(1, 0, 0);
//    glVertex3f(0, 0, 0);
//    glVertex3f(1, 0, 0);
//    glColor3f(0, 1, 0);
//    glVertex3f(0, 0, 0);
//    glVertex3f(0, 1, 0);
//    glColor3f(0, 0, 1);
//    glVertex3f(0, 0, 0);
//    glVertex3f(0, 0, 1);
//    glEnd();
    
    // render the mesh
    g_view.renderer.render();    
    
    glPopMatrix();
    
    glutSwapBuffers();
}

void reshape(int w, int h)
{
    glViewport(0, 0, w, h);
    
    g_view.win_w = w;
    g_view.win_h = h;
    double ar = (double)w / h;
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-0.1 * ar, 0.1 * ar, -0.1, 0.1, 0.3, 10);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef(0, 0, -g_view.dist);
    glRotatef(g_view.alpha * 180 / M_PI, 1, 0, 0);
    glRotatef(g_view.theta * 180 / M_PI, 0, 1, 0);
    glTranslatef(-0.5, -0.5, -0.5);
    
    glEnable(GL_DEPTH_TEST);
}

void idle()
{

}

namespace
{
    class Vec2iLess
    {
    public:
        bool operator () (const LosTopos::Vec2i & x, const LosTopos::Vec2i & y) const { return (x[0] < y[0] || (x[0] == y[0] && x[1] < y[1])); }
    };
}

void keyboard(unsigned char k, int x, int y)
{
    if (k == 27 || k == 'q' || k == 'Q')
    {
        exit(0);
    } else if (k == '[' || k == '{')
    {
        // prev frame
        g_playback.current_frame--;
        if (g_playback.current_frame < 0) g_playback.current_frame = 0;
        g_playback.current_subframe = 0;

        if (k == '[')
        {
            std::stringstream ss;
            ss << g_playback.rec_directory << "/mesh" << std::setfill('0') << std::setw(6) << g_playback.current_frame << ".rec";
            bool success = MeshIO::load(*g_playback.st, ss.str());
            if (success) std::cout << "Frame " << ss.str() << " loaded: nv = " << g_playback.st->m_mesh.nv() << " ne = " << g_playback.st->m_mesh.ne() << " nf = " << g_playback.st->m_mesh.nt() << std::endl;
            glutPostRedisplay();
        }
        
        std::cout << "Current frame: " << g_playback.current_frame << std::endl;
        
    } else if (k == ']' || k == '}')
    {
        // next frame
        g_playback.current_frame++;
        g_playback.current_subframe = 0;
        
        if (k == ']')
        {
            std::stringstream ss;
            ss << g_playback.rec_directory << "/mesh" << std::setfill('0') << std::setw(6) << g_playback.current_frame << ".rec";
            bool success = MeshIO::load(*g_playback.st, ss.str());
            if (success) std::cout << "Frame " << ss.str() << " loaded: nv = " << g_playback.st->m_mesh.nv() << " ne = " << g_playback.st->m_mesh.ne() << " nf = " << g_playback.st->m_mesh.nt() << std::endl;
            glutPostRedisplay();
        }
        
        std::cout << "Current frame: " << g_playback.current_frame << std::endl;
        
    } else if (k == ',' || k == '<')
    {
        // prev subframe
        g_playback.current_subframe--;
        
        std::stringstream ss;
        ss << g_playback.rec_directory << "/mesh" << std::setfill('0') << std::setw(6) << g_playback.current_frame << "." << g_playback.current_subframe << ".rec";
        bool success = MeshIO::load(*g_playback.st, ss.str());
        if (success) std::cout << "Subframe " << ss.str() << " loaded: nv = " << g_playback.st->m_mesh.nv() << " ne = " << g_playback.st->m_mesh.ne() << " nf = " << g_playback.st->m_mesh.nt() << std::endl;
        glutPostRedisplay();
        
        std::cout << "Current frame: " << g_playback.current_frame << " subframe: " << g_playback.current_subframe << std::endl;

    } else if (k == '.' || k == '>')
    {
        // next subframe
        g_playback.current_subframe++;
        
        std::stringstream ss;
        ss << g_playback.rec_directory << "/mesh" << std::setfill('0') << std::setw(6) << g_playback.current_frame << "." << g_playback.current_subframe << ".rec";
        bool success = MeshIO::load(*g_playback.st, ss.str());
        if (success) std::cout << "Subframe " << ss.str() << " loaded: nv = " << g_playback.st->m_mesh.nv() << " ne = " << g_playback.st->m_mesh.ne() << " nf = " << g_playback.st->m_mesh.nt() << std::endl;
        glutPostRedisplay();
        
        std::cout << "Current frame: " << g_playback.current_frame << " subframe: " << g_playback.current_subframe << std::endl;
        
    } else if (k == 'n' || k == 'N')
    {
        int v = g_view.renderer.nearestVertex();
        int e = g_view.renderer.nearestEdge();
        int f = g_view.renderer.nearestFace();
        
        LosTopos::NonDestructiveTriMesh & mesh = g_playback.st->m_mesh;
        
        if (v >= 0)
        {
            std::cout << "VOI: " << v << " (" << g_playback.st->get_position(v) << ")" << std::endl;
            std::cout << "  number of incident edges = " << mesh.m_vertex_to_edge_map[v].size() << " number of incident triangles = " << mesh.m_vertex_to_triangle_map[v].size() << std::endl;
            std::cout << "  Region graph: " << std::endl;
            
            std::set<int> regionset;
            std::map<LosTopos::Vec2i, bool, Vec2iLess> graph;
            for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[v].size(); i++)
            {
                LosTopos::Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[v][i]);
                if (label[1] < label[0]) std::swap(label[0], label[1]);
                regionset.insert(label[0]);
                regionset.insert(label[1]);
                graph[label] = true;
            }
            
            std::vector<int> regions;
            regions.assign(regionset.begin(), regionset.end());
            std::sort(regions.begin(), regions.end());
            std::cout << "   "; for (size_t i = 0; i < regions.size(); i++) std::cout << std::setfill(' ') << std::setw(2) << regions[i] << " "; std::cout << std::endl;
            for (size_t i = 0; i < regions.size(); i++)
            {
                std::cout << std::setfill(' ') << std::setw(2) << regions[i] << " ";
                for (size_t j = 0; j < regions.size(); j++)
                {
                    LosTopos::Vec2i regionpair = (regions[i] < regions[j] ? LosTopos::Vec2i(regions[i], regions[j]) : LosTopos::Vec2i(regions[j], regions[i]));
                    std::cout << " " << (graph[regionpair] ? "*" : " ") << " ";
                }
                std::cout << std::endl;
            }
        }
        if (e >= 0)
        {
            std::cout << "EOI: " << e << ": " << mesh.m_edges[e][0] << " (" << g_playback.st->get_position(mesh.m_edges[e][0]) << ") - " << mesh.m_edges[e][1] << " (" << g_playback.st->get_position(mesh.m_edges[e][1]) << ")" << std::endl;
            std::cout << "  length = " << g_playback.st->get_edge_length(e) << " number of incident triangles = " << mesh.m_edge_to_triangle_map[e].size() << std::endl;
        }
        if (f >= 0)
        {
            size_t v0 = mesh.m_tris[f][0];
            size_t v1 = mesh.m_tris[f][1];
            size_t v2 = mesh.m_tris[f][2];
            std::cout << "FOI: " << f << ": " << v0 << " (" << g_playback.st->get_position(v0) << "), " << v1 << " (" << g_playback.st->get_position(v1) << "), " << v2 << " (" << g_playback.st->get_position(v2) << ")" << std::endl;
            std::cout << "  area = " << g_playback.st->get_triangle_area(f) << " label = " << mesh.get_triangle_label(f) << std::endl;
        }

    }
    
    g_view.renderer.keyboard(k, x, y);
}

void mouse(int button, int state, int x, int y)
{
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
    {
        g_view.ldrag = true;
        g_view.x_ldrag = x;
        g_view.y_ldrag = y;
        
    } else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP)
    {
        g_view.ldrag = false;
    } else if (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
    {
        g_view.rdrag = true;
        g_view.x_rdrag = x;
        g_view.y_rdrag = y;
        
    } else if (button == GLUT_RIGHT_BUTTON && state == GLUT_UP)
    {
        g_view.rdrag = false;
    }
    
    glutPostRedisplay();
}

void motion(int x, int y)
{
    if (g_view.ldrag)
    {
        g_view.alpha += (y - g_view.y_ldrag) * View::LDRAG_GAIN;
        g_view.theta += (x - g_view.x_ldrag) * View::LDRAG_GAIN;
        
        g_view.x_ldrag = x;
        g_view.y_ldrag = y;
    }
    
    if (g_view.rdrag)
    {
        g_view.dist *= exp((y - g_view.y_rdrag) * View::RDRAG_GAIN);
        
        g_view.x_rdrag = x;
        g_view.y_rdrag = y;
    }

    g_view.renderer.setMousePos(-1 + 2 * (double)x / g_view.win_w, 1 - 2 * (double)y / g_view.win_h);
    glutPostRedisplay();
}

void passiveMotion(int x, int y)
{
    g_view.renderer.setMousePos(-1 + 2 * (double)x / g_view.win_w, 1 - 2 * (double)y / g_view.win_h);
    glutPostRedisplay();
}

int main(int argc, char * argv[])
{
    if (argc != 2)
    {
        std::cout << "Usage:\n\tMultiTrackerViewer rec_directory" << std::endl;
        exit(1);
    }
    
    g_playback.rec_directory = argv[1];
    
    glutInit(&argc, argv);
    
    glutInitWindowPosition(0, 0);
    glutInitWindowSize(800, 800);
    glutCreateWindow("MultiTrackerViewer");
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutIdleFunc(idle);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    glutPassiveMotionFunc(passiveMotion);
    
    g_view.alpha = 0;
    g_view.theta = 0;
    g_view.dist = 3;
    g_view.ldrag = false;
    g_view.x_ldrag = g_view.y_ldrag = 0;
    g_view.rdrag = false;
    g_view.x_rdrag = g_view.y_rdrag = 0;
    g_view.renderer.setMousePos(0.5, 0.5);

    g_playback.current_frame = 0;
    
    LosTopos::SurfTrackInitializationParameters stip;
    // set the required parameters just to be able to construct the SurfTrack object; their values don't matter because the SurfTrack object is only used as a mesh data container.
    stip.m_min_edge_length = 1;
    stip.m_max_edge_length = 1;
    stip.m_max_volume_change = std::numeric_limits<double>::max();
    std::vector<LosTopos::Vec3d> vs(3); vs[0] = LosTopos::Vec3d(0, 0, 0); vs[1] = LosTopos::Vec3d(1, 0, 0); vs[2] = LosTopos::Vec3d(0, 1, 0);
    std::vector<LosTopos::Vec3d> ms(3); ms[0] = LosTopos::Vec3d(1, 1, 1); ms[1] = LosTopos::Vec3d(1, 1, 1); ms[2] = LosTopos::Vec3d(1, 1, 1);
    std::vector<LosTopos::Vec3st> fs(1); fs[0] = LosTopos::Vec3st(0, 1, 2);
    std::vector<LosTopos::Vec2i> ls(1); ls[0] = LosTopos::Vec2i(0, 1);
    g_playback.st = new LosTopos::SurfTrack(vs, fs, ls, ms, stip);
    bool success = MeshIO::load(*g_playback.st, g_playback.rec_directory + "/mesh000000.rec");
    if (success) std::cout << "Frame " << g_playback.rec_directory << "/mesh000000.rec loaded: nv = " << g_playback.st->m_mesh.nv() << " ne = " << g_playback.st->m_mesh.ne() << " nf = " << g_playback.st->m_mesh.nt() << std::endl;
    
    int nregions = 0;
    for (size_t i = 0; i < g_playback.st->m_mesh.nt(); i++)
    {
        const LosTopos::Vec2i & l = g_playback.st->m_mesh.get_triangle_label(i);
        if (l[0] >= nregions) nregions = l[0] + 1;
        if (l[1] >= nregions) nregions = l[1] + 1;
    }
    g_view.renderer.setSurfTrack(g_playback.st, nregions);
    
    glutMainLoop();
    
	return 0;
}

