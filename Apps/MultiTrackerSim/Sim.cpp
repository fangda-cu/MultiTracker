//
//  Sim.cpp
//
//  Christopher Batty, Fang Da 2014
//
//

#include <sstream>
#include <iomanip>
#include "Sim.h"
#include "SimOptions.h"
#include "subdivisionscheme.h"
#include "MeshIO.h"
#include "drivers/MCFDriver.h"
#include "drivers/MultimaterialFaceOffsetDriver.h"
#include "drivers/CurlNoiseDriver.h"

using namespace LosTopos;

Sim::Sim(bool verbose) :
    m_verbose(verbose),
    m_scene("unspecified"),
    m_output_directory(""),
    m_assets_directory(""),
    m_bbwall(false),
    m_t1vel(false),
    m_dt(0),
    m_time(0),
    m_frameid(0),
    m_finished(false),
    m_st(NULL)
{
    
}

Sim::~Sim()
{
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::step()
{
    assert(m_scene != "unspecified");
    assert(m_st);
    if (m_verbose)
        std::cout << "Time stepping: t = " << m_time << ", dt = " << m_dt << std::endl;
    
    
    // scene-specific time stepping
    if (m_scene == "T1")
        stepT1(m_dt);
    else if (m_scene == "T2")
        stepT2(m_dt);
    else if (m_scene == "Merge")
        stepMerge(m_dt);
    else if (m_scene == "Zalesak")
        stepZalesak(m_dt);
    else if (m_scene == "Enright")
        stepEnright(m_dt);
    else if (m_scene == "CyclicFlow")
        stepCyclicFlow(m_dt);
    else if (m_scene == "Noise")
        stepNoise(m_dt);
    else if (m_scene == "MCF")
        stepMCF(m_dt);
    
    
    // handle collisions during the dynamics
    double actual_dt;
    m_st->integrate(m_dt, actual_dt);
    if (actual_dt != m_dt)
        std::cout << "Warning: SurfTrack::integrate() failed to step the full length of the time step!" << std::endl;
    
    
    // BB wall constraint: update the infinite masses
    updateBBWallConstraints();
    
    
    // remove faces completely inside the BB wall
    removeBBWallFaces();

    
    // mesh improvement
    for(int i = 0; i < Options::intValue("remeshing-iterations"); i++)
    {
        m_st->topology_changes();
        m_st->improve_mesh();
    }
    
    
    // remove faces completely inside the BB wall
    removeBBWallFaces();
    
    
    // defrag the mesh in the end, to ensure the next step starts with a clean mesh
    m_st->defrag_mesh();
    
    
    // advance time
    m_frameid++;
    m_time += m_dt;
    if (m_time >= Options::doubleValue("simulation-time"))
        m_finished = true;
   
   
    // output the mesh
    std::stringstream filename;
    filename << m_output_directory << "/mesh" << std::setfill('0') << std::setw(6) << m_frameid << ".rec";
    MeshIO::save(*m_st, filename.str());
   
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  General initialization of a simulation
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::init(const std::string & option_file, const std::string & output_directory, const std::string & assets_directory)
{
    // declare and load the options
    Options::addStringOption("scene", "T1");
    Options::addDoubleOption("time-step", 0.01);
    Options::addDoubleOption("simulation-time", 1.0);
    Options::addDoubleOption("remeshing-resolution", 0.1);
    Options::addIntegerOption("remeshing-iterations", 1);
    
    Options::addDoubleOption ("lostopos-collision-epsilon-fraction", 1e-4);       // lostopos collision epsilon (fraction of mean edge length)
    Options::addDoubleOption ("lostopos-merge-proximity-epsilon-fraction", 0.02); // lostopos merge proximity epsilon (fraction of mean edge length)
    Options::addBooleanOption("lostopos-perform-smoothing", false);               // whether or not to perform smoothing
    Options::addDoubleOption ("lostopos-max-volume-change-fraction", 1e-4);       // maximum allowed volume change during a remeshing operation (fraction of mean edge length cubed)
    Options::addDoubleOption ("lostopos-min-triangle-angle", 3.0);                // min triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-max-triangle-angle", 177.0);              // max triangle angle (in degrees)
    Options::addDoubleOption ("lostopos-large-triangle-angle-to-split", 160.0);   // threshold for large angles to be split
    Options::addDoubleOption ("lostopos-min-triangle-area-fraction", 0.02);       // minimum allowed triangle area (fraction of mean edge length squared)
    Options::addBooleanOption("lostopos-t1-transition-enabled", true);            // whether t1 is enabled
    Options::addDoubleOption ("lostopos-t1-pull-apart-distance-fraction", 0.1);   // t1 pull apart distance (fraction of mean edge legnth)
    Options::addBooleanOption("lostopos-smooth-subdivision", false);              // whether to use smooth subdivision during remeshing
    Options::addBooleanOption("lostopos-allow-non-manifold", true);               // whether to allow non-manifold geometry in the mesh
    Options::addBooleanOption("lostopos-allow-topology-changes", true);           // whether to allow topology changes
    
    Options::parseOptionFile(option_file, m_verbose);
    

    // select the scene
    m_scene = Options::strValue("scene");
    m_output_directory = output_directory;
    m_assets_directory = assets_directory;
    
    std::vector<Vec3d> vertices;
    std::vector<Vec3st> faces;
    std::vector<Vec2i> face_labels;
    std::vector<Vec3d> masses;
    
    if (m_scene == "T1")
        sceneT1(vertices, faces, face_labels);
    else if (m_scene == "T2")
        sceneT2(vertices, faces, face_labels);
    else if (m_scene == "Merge")
        sceneMerge(vertices, faces, face_labels);
    else if (m_scene == "Zalesak")
        sceneZalesak(vertices, faces, face_labels);
    else if (m_scene == "Enright")
        sceneEnright(vertices, faces, face_labels);
    else if (m_scene == "CyclicFlow")
        sceneCyclicFlow(vertices, faces, face_labels);
    else if (m_scene == "Noise")
        sceneNoise(vertices, faces, face_labels);
    else if (m_scene == "MCF")
        sceneMCF(vertices, faces, face_labels);

    
    // construct the surface tracker
    double mean_edge_len = Options::doubleValue("remeshing-resolution");
    double min_edge_len = mean_edge_len * 0.5;
    double max_edge_len = mean_edge_len * 1.5;
    
    SurfTrackInitializationParameters params;
    params.m_proximity_epsilon = Options::doubleValue("lostopos-collision-epsilon-fraction") * mean_edge_len;
    params.m_merge_proximity_epsilon = Options::doubleValue("lostopos-merge-proximity-epsilon-fraction") * mean_edge_len;
    params.m_allow_vertex_movement_during_collapse = true;
    params.m_perform_smoothing = Options::boolValue("lostopos-perform-smoothing");
    params.m_min_edge_length = min_edge_len;
    params.m_max_edge_length = max_edge_len;
    params.m_max_volume_change = Options::doubleValue("lostopos-max-volume-change-fraction") * pow(mean_edge_len, 3);
    params.m_min_triangle_angle = Options::doubleValue("lostopos-min-triangle-angle");
    params.m_max_triangle_angle = Options::doubleValue("lostopos-max-triangle-angle");
    params.m_large_triangle_angle_to_split = Options::doubleValue("lostopos-large-triangle-angle-to-split");
    params.m_min_triangle_area = Options::doubleValue("lostopos-min-triangle-area-fraction") * pow(mean_edge_len, 2);
    params.m_verbose = false;
    params.m_allow_non_manifold = Options::boolValue("lostopos-allow-non-manifold");
    params.m_allow_topology_changes = Options::boolValue("lostopos-allow-topology-changes");
    params.m_collision_safety = true;
    params.m_remesh_boundaries = true;
    params.m_t1_transition_enabled = Options::boolValue("lostopos-t1-transition-enabled");
    params.m_pull_apart_distance = Options::doubleValue("lostopos-t1-pull-apart-distance-fraction") * mean_edge_len;
    
    params.m_velocity_field_callback = NULL;
    if (m_t1vel)
        params.m_velocity_field_callback = this; // this is only turned on for specific scenes
    
    if (Options::boolValue("lostopos-smooth-subdivision"))
        params.m_subdivision_scheme = new ModifiedButterflyScheme();
    else
        params.m_subdivision_scheme = new MidpointScheme();
    
    params.m_use_curvature_when_collapsing = false;
    params.m_use_curvature_when_splitting = false;

    masses.resize(vertices.size(), Vec3d(1, 1, 1));
    m_st = new SurfTrack(vertices, faces, face_labels, masses, params);
    if (m_bbwall) m_st->m_solid_vertices_callback = this;
    
    
    // BB wall constraint: update the infinite masses
    updateBBWallConstraints();

    
    // prepare to start the simulation
    m_time = 0;
    m_dt = Options::doubleValue("time-step");
    m_finished = false;
    
    
    // output the initial frame
    MeshIO::save(*m_st, m_output_directory + "/mesh000000.rec");
    
    return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific initialization
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::sceneT1(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = true;

    vs.resize(12);
    vs[ 0] = Vec3d(0.3, 0.0, 0.0);
    vs[ 1] = Vec3d(0.3, 0.5, 0.0);
    vs[ 2] = Vec3d(0.3, 1.0, 0.0);
    vs[ 3] = Vec3d(0.3, 0.0, 1.0);
    vs[ 4] = Vec3d(0.3, 0.5, 1.0);
    vs[ 5] = Vec3d(0.3, 1.0, 1.0);
    vs[ 6] = Vec3d(0.7, 0.0, 0.0);
    vs[ 7] = Vec3d(0.7, 0.5, 0.0);
    vs[ 8] = Vec3d(0.7, 1.0, 0.0);
    vs[ 9] = Vec3d(0.7, 0.0, 1.0);
    vs[10] = Vec3d(0.7, 0.5, 1.0);
    vs[11] = Vec3d(0.7, 1.0, 1.0);

    fs.resize(10); ls.resize(10);
    fs[ 0] = Vec3st( 0,  1,  3);    ls[ 0] = Vec2i(0, 2);
    fs[ 1] = Vec3st( 4,  3,  1);    ls[ 1] = Vec2i(0, 2);
    fs[ 2] = Vec3st( 1,  2,  4);    ls[ 2] = Vec2i(0, 3);
    fs[ 3] = Vec3st( 5,  4,  2);    ls[ 3] = Vec2i(0, 3);
    fs[ 4] = Vec3st( 6,  7,  9);    ls[ 4] = Vec2i(2, 1);
    fs[ 5] = Vec3st(10,  9,  7);    ls[ 5] = Vec2i(2, 1);
    fs[ 6] = Vec3st( 7,  8, 10);    ls[ 6] = Vec2i(3, 1);
    fs[ 7] = Vec3st(11, 10,  8);    ls[ 7] = Vec2i(3, 1);
    fs[ 8] = Vec3st( 1,  4,  7);    ls[ 8] = Vec2i(2, 3);
    fs[ 9] = Vec3st(10,  7,  4);    ls[ 9] = Vec2i(2, 3);
    
}

void Sim::sceneT2(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = true;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/T2.arec", false);
}

void Sim::sceneMerge(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = false;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/merge.arec", false);
}

void Sim::sceneZalesak(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = false;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/zalesak.arec", false);
}

void Sim::sceneEnright(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = false;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/enright.arec", false);
}

void Sim::sceneCyclicFlow(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = false;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/cyclicflow.arec", false);
}

void Sim::sceneNoise(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = false;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/noise.arec", false);
}

void Sim::sceneMCF(std::vector<Vec3d> & vs, std::vector<Vec3st> & fs, std::vector<Vec2i> & ls)
{
    m_bbwall = true;
    MeshIO::loadIntoRaw(vs, fs, ls, m_assets_directory + "/mcf.arec", false);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Scene-specific time stepping
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::stepT1(double dt)
{
    MCFDriver::step(m_st, dt, m_bbwall);
}

void Sim::stepT2(double dt)
{
    MCFDriver::step(m_st, dt, m_bbwall);
    
    // constraint the lateral motion of wall vertices
    for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
        if (onBBWall(m_st->get_position(i)))
            m_st->set_newposition(i, m_st->get_position(i));
}

void Sim::stepMerge(double dt)
{
    std::vector<std::vector<double> > speeds(3, std::vector<double>(3));
    if (m_time < 0.15)
    {
        speeds[0][0] = 0;   speeds[0][1] = -1;  speeds[0][2] = -1;
        speeds[1][0] = 1;   speeds[1][1] = 0;   speeds[1][2] = 0;
        speeds[2][0] = 1;   speeds[2][1] = 0;   speeds[2][2] = 0;
    } else
    {
        speeds[0][0] = 0;   speeds[0][1] = 1;   speeds[0][2] = 1;
        speeds[1][0] = -1;  speeds[1][1] = 0;   speeds[1][2] = 0;
        speeds[2][0] = -1;  speeds[2][1] = 0;   speeds[2][2] = 0;
    }
    
    MultimaterialFaceOffsetDriver::step(m_st, dt, speeds, 0, false);
}

namespace
{
    void zalesak_velocity(double t, const Vec3d & pos, Vec3d & out)
    {
        double x = pos[0];
        double y = pos[1];
        double z = pos[2];
        
        out = Vec3d((z - 0.5), 0, -(x - 0.5));
    }
}

void Sim::stepZalesak(double dt)
{
    for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
    {
        // RK4 integration
        Vec3d v;
        Vec3d x = m_st->get_position(i);
        zalesak_velocity(m_time, x, v);
        Vec3d k1 = v;
        zalesak_velocity(m_time + 0.5 * dt, x + 0.5 * dt * k1, v);
        Vec3d k2 = v;
        zalesak_velocity(m_time + 0.5 * dt, x + 0.5 * dt * k2, v);
        Vec3d k3 = v;
        zalesak_velocity(m_time + dt, x + dt * k3, v);
        Vec3d k4 = v;
        v = (1./6. * (k1 + k4) + 1./3. * (k2 + k3));
        m_st->set_newposition(i, x + v * dt);
    }
}

namespace
{
    void enright_velocity(double t, const Vec3d & pos, Vec3d & out)
    {
        //
        //  code adapted from El Topo's Enright driver:
        //    https://github.com/tysonbrochu/eltopo/blob/master/talpa/drivers/enrightdriver.h
        //
        double x = pos[0];
        double y = pos[1];
        double z = pos[2];
        
        out = Vec3d(2.0 * std::sin(    M_PI*x) * std::sin(    M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(2.0*M_PI*z),
                         -std::sin(2.0*M_PI*x) * std::sin(    M_PI*y) * std::sin(    M_PI*y) * std::sin(2.0*M_PI*z),
                         -std::sin(2.0*M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(    M_PI*z) * std::sin(    M_PI*z));
        
        out *= sin(M_PI * t * 2 / 3);    // modulate with a period of 3
    }
}

void Sim::stepEnright(double dt)
{
    //
    //  RK4 integration of the Enright velocity field
    //  code adapted from El Topo's Enright driver:
    //    https://github.com/tysonbrochu/eltopo/blob/master/talpa/drivers/enrightdriver.h
    //
    for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
    {
        Vec3d v;
        Vec3d x = m_st->get_position(i);
        enright_velocity(m_time, x, v);
        Vec3d k1 = v;
        enright_velocity(m_time + 0.5 * dt, x + 0.5 * dt * k1, v);
        Vec3d k2 = v;
        enright_velocity(m_time + 0.5 * dt, x + 0.5 * dt * k2, v);
        Vec3d k3 = v;
        enright_velocity(m_time + dt, x + dt * k3, v);
        Vec3d k4 = v;
        v = (1./6. * (k1 + k4) + 1./3. * (k2 + k3));
        m_st->set_newposition(i, x + v * dt);
    }
}

void Sim::stepCyclicFlow(double dt)
{
    std::vector<std::vector<double> > speeds(3, std::vector<double>(3));
    speeds[0][0] = 0;   speeds[0][1] = 1;   speeds[0][2] = -1;
    speeds[1][0] = -1;  speeds[1][1] = 0;   speeds[1][2] = 1;
    speeds[2][0] = 1;   speeds[2][1] = -1;  speeds[2][2] = 0;
    
    MultimaterialFaceOffsetDriver::step(m_st, dt, speeds, -1, true);
}

void Sim::stepNoise(double dt)
{
    static bool stepfaceoff = true;
    
    if (stepfaceoff)
    {
        std::vector<std::vector<double> > speeds(5, std::vector<double>(5));
        speeds[0][0] = 0;   speeds[0][1] = 1;   speeds[0][2] = 1;   speeds[0][3] = 1;   speeds[0][4] = 1;
        speeds[1][0] = -1;  speeds[1][1] = 0;   speeds[1][2] = 0;   speeds[1][3] = 0;   speeds[1][4] = 0;
        speeds[2][0] = -1;  speeds[2][1] = 0;   speeds[2][2] = 0;   speeds[2][3] = 0;   speeds[2][4] = 0;
        speeds[3][0] = -1;  speeds[3][1] = 0;   speeds[3][2] = 0;   speeds[3][3] = 0;   speeds[3][4] = 0;
        speeds[4][0] = -1;  speeds[4][1] = 0;   speeds[4][2] = 0;   speeds[4][3] = 0;   speeds[4][4] = 0;
        
        MultimaterialFaceOffsetDriver::step(m_st, dt, speeds, 0, false);
    } else
    {
        CurlNoiseDriver::step(m_st, dt);
    }
    
    stepfaceoff = !stepfaceoff;
}

void Sim::stepMCF(double dt)
{
    MCFDriver::step(m_st, dt, m_bbwall);
   
   static bool s_speedup_1 = false;
   static bool s_speedup_2 = false;
   if (m_time > 6 && !s_speedup_1)
   {
      s_speedup_1 = true;
      m_dt *= 4;
   }
   if (m_time > 30 && !s_speedup_2)
   {
      s_speedup_2 = true;
      m_dt *= 4;
   }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Utilities for maintaining a [0, 1]^3 bounding box
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Sim::updateBBWallConstraints()
{
    if (m_bbwall)
    {
        for (size_t i = 0; i < m_st->m_mesh.nv(); i++)
        {
            Vec3d & x = m_st->pm_positions[i];
            if (x[0] > 1) x[0] = 1;
            if (x[0] < 0) x[0] = 0;
            if (x[1] > 1) x[1] = 1;
            if (x[1] < 0) x[1] = 0;
            if (x[2] > 1) x[2] = 1;
            if (x[2] < 0) x[2] = 0;
            
            int onwall = onBBWall(x);
            Vec3d & mass = m_st->m_masses[i];
            mass[0] = ((onwall & (1 << 0)) || (onwall & (1 << 3)) ? std::numeric_limits<double>::infinity() : 1);
            mass[1] = ((onwall & (1 << 1)) || (onwall & (1 << 4)) ? std::numeric_limits<double>::infinity() : 1);
            mass[2] = ((onwall & (1 << 2)) || (onwall & (1 << 5)) ? std::numeric_limits<double>::infinity() : 1);
        }
    }
}

void Sim::removeBBWallFaces()
{
    if (!m_bbwall) // scene doesn't use the BB
        return;

    // remove all faces completely inside wall
    for (size_t i = 0; i < m_st->m_mesh.nt(); i++)
    {
        Vec3st & f = m_st->m_mesh.m_tris[i];
        int w0 = onBBWall(m_st->get_position(f[0]));
        int w1 = onBBWall(m_st->get_position(f[1]));
        int w2 = onBBWall(m_st->get_position(f[2]));
        
        if (w0 & w1 & w2)
            m_st->remove_triangle(i);
    }
}

int Sim::onBBWall(const LosTopos::Vec3d & pos) const
{
    if (!m_bbwall) // scene doesn't use the BB
        return 0;
    
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

Vec3d Sim::enforceBBWallConstraint(const Vec3d & input, int constraints) const
{
    if (!m_bbwall) // scene doesn't use the BB
        return input;
    
    Vec3d output = input;
    if (constraints & (1 << 0)) output[0] = 0;
    if (constraints & (1 << 1)) output[1] = 0;
    if (constraints & (1 << 2)) output[2] = 0;
    if (constraints & (1 << 3)) output[0] = 1;
    if (constraints & (1 << 4)) output[1] = 1;
    if (constraints & (1 << 5)) output[2] = 1;
    
    return output;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  Callbacks
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Sim::generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    LosTopos::Vec3d x0 = st.get_position(v0);
    LosTopos::Vec3d x1 = st.get_position(v1);
    
    int label0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int label1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    if (label0 == label1)
    {
        // on the same wall(s), prefer the one with higher max edge valence
        size_t maxedgevalence0 = 0;
        size_t maxedgevalence1 = 0;
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v0].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size() > maxedgevalence0)
                maxedgevalence0 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size();
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v1].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size() > maxedgevalence1)
                maxedgevalence1 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size();
        
        if (maxedgevalence0 == maxedgevalence1) // same max edge valence, use their midpoint
            pos = (x0 + x1) / 2;
        else if (maxedgevalence0 < maxedgevalence1)
            pos = x1;
        else
            pos = x0;
        
        return true;
        
    } else if ((label0 & ~label1) == 0)
    {
        // label0 is a proper subset of label1 (since label0 != label1)
        pos = x1;
        
        return true;
        
    } else if ((label1 & ~label0) == 0)
    {
        // label1 is a proper subset of label0
        pos = x0;
        
        return true;
        
    } else
    {
        // label0 and label1 are not subset of each other
        int newlabel = label0 | label1;
        assert(label0 != newlabel); // not subset of each other
        assert(label1 != newlabel);
        
        assert(!((label0 & (1 << 0)) != 0 && (label0 & (1 << 3)) != 0)); // can't have conflicting constraints in label0 and label1 already
        assert(!((label0 & (1 << 1)) != 0 && (label0 & (1 << 4)) != 0));
        assert(!((label0 & (1 << 2)) != 0 && (label0 & (1 << 5)) != 0));
        assert(!((label1 & (1 << 0)) != 0 && (label1 & (1 << 3)) != 0));
        assert(!((label1 & (1 << 1)) != 0 && (label1 & (1 << 4)) != 0));
        assert(!((label1 & (1 << 2)) != 0 && (label1 & (1 << 5)) != 0));
        
        bool conflict = false;
        if ((newlabel & (1 << 0)) != 0 && (newlabel & (1 << 3)) != 0) conflict = true;
        if ((newlabel & (1 << 1)) != 0 && (newlabel & (1 << 4)) != 0) conflict = true;
        if ((newlabel & (1 << 2)) != 0 && (newlabel & (1 << 5)) != 0) conflict = true;
        
        if (conflict)
        {
            // the two vertices are on opposite walls (conflicting constraints). Can't collapse this edge (which shouldn't have become a collapse candidate in the first place)
            return false;
        }
        
        pos = (x0 + x1) / 2;
        if (newlabel & (1 << 0))  pos[0] = 0;    // project the midpoint onto the constraint manifold (BB walls)
        if (newlabel & (1 << 1))  pos[1] = 0;
        if (newlabel & (1 << 2))  pos[2] = 0;
        if (newlabel & (1 << 3))  pos[0] = 1;
        if (newlabel & (1 << 4))  pos[1] = 1;
        if (newlabel & (1 << 5))  pos[2] = 1;
        
        return true;
    }
    
    return false;
}

bool Sim::generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos)
{
    pos = (st.get_position(v0) + st.get_position(v1)) / 2;
    
    return true;
}

LosTopos::Vec3c Sim::generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    LosTopos::Vec3d x0 = st.get_position(v0);
    LosTopos::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    LosTopos::Vec3c result;  // if either endpoint is constrained, the collapsed point shold be constrained. more specifically it should be on all the walls any of the two endpoints is on (implemented in generate_collapsed_position())
    int result_constraint = (constraint0 | constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
}

LosTopos::Vec3c Sim::generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1)
{
    LosTopos::Vec3d x0 = st.get_position(v0);
    LosTopos::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
  
    LosTopos::Vec3c result;  // the splitting midpoint has a positive constraint label only if the two endpoints are on a same wall (sharing a bit in their constraint bitfield representation)
    int result_constraint = (constraint0 & constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
    
}

bool Sim::generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower)
{
    LosTopos::Vec3d original_pos = st.get_position(oldv);
    int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
    
    Vec3d new_pos_upper = enforceBBWallConstraint(Vec3d(pos_upper[0], pos_upper[1], pos_upper[2]), original_constraint);
    Vec3d new_pos_lower = enforceBBWallConstraint(Vec3d(pos_lower[0], pos_lower[1], pos_lower[2]), original_constraint);
    
    pos_upper = new_pos_upper;
    pos_lower = new_pos_lower;
    
    return true;
}

bool Sim::generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b)
{
    LosTopos::Vec3d original_pos = st.get_position(oldv);
    int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
    
    Vec3d new_pos_a = enforceBBWallConstraint(Vec3d(pos_a[0], pos_a[1], pos_a[2]), original_constraint);
    Vec3d new_pos_b = enforceBBWallConstraint(Vec3d(pos_b[0], pos_b[1], pos_b[2]), original_constraint);
    
    pos_a = new_pos_a;
    pos_b = new_pos_b;
    
    return true;
}

bool Sim::solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e)
{
    LosTopos::Vec3d x0 = st.get_position(st.m_mesh.m_edges[e][0]);
    LosTopos::Vec3d x1 = st.get_position(st.m_mesh.m_edges[e][1]);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    if (constraint0 & constraint1)  // edge is completely inside a wall
        return true;
    else
        return false;
}

LosTopos::Vec3d Sim::sampleVelocity(LosTopos::Vec3d & pos)
{
    return Vec3d(0, 0, 0);
}

bool Sim::sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output)
{
    return false;
}




