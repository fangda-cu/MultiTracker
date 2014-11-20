//
//  Sim.h
//
//  Christopher Batty, Fang Da 2014
//
//

#ifndef __Sim__
#define __Sim__

#include <iostream>
#include <string>
#include "surftrack.h"

class Sim : public LosTopos::SurfTrack::SolidVerticesCallback, LosTopos::T1Transition::VelocityFieldCallback
{
public:
    Sim(bool verbose);
    ~Sim();
    
public:
    bool init(const std::string & option_file, const std::string & output_directory, const std::string & assets_directory);
    
public:
    void step();
    bool isFinished() const { return m_finished; }
    
protected:
    // scene-specific initialization
    void sceneT1        (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneT2        (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneMerge     (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneZalesak   (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneEnright   (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneCyclicFlow(std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneNoise     (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    void sceneMCF       (std::vector<LosTopos::Vec3d> & vs, std::vector<LosTopos::Vec3st> & fs, std::vector<LosTopos::Vec2i> & ls);
    
    // scene-specific time evolution
    void stepT1         (double dt);
    void stepT2         (double dt);
    void stepMerge      (double dt);
    void stepZalesak    (double dt);
    void stepEnright    (double dt);
    void stepCyclicFlow (double dt);
    void stepNoise      (double dt);
    void stepMCF        (double dt);
    
protected:
    void updateBBWallConstraints();
    void removeBBWallFaces();
    
    int onBBWall(const LosTopos::Vec3d & pos) const;
    LosTopos::Vec3d enforceBBWallConstraint(const LosTopos::Vec3d & input, int constraints) const;

protected:
    // SurfTrack::SolidVerticesCallback method
    bool            generate_collapsed_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    bool            generate_split_position(LosTopos::SurfTrack & st, size_t v0, size_t v1, LosTopos::Vec3d & pos);
    LosTopos::Vec3c generate_collapsed_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    LosTopos::Vec3c generate_split_solid_label(LosTopos::SurfTrack & st, size_t v0, size_t v1, const LosTopos::Vec3c & label0, const LosTopos::Vec3c & label1);
    bool            generate_edge_popped_positions(LosTopos::SurfTrack & st, size_t oldv, const LosTopos::Vec2i & cut, LosTopos::Vec3d & pos_upper, LosTopos::Vec3d & pos_lower);
    bool            generate_vertex_popped_positions(LosTopos::SurfTrack & st, size_t oldv, int A, int B, LosTopos::Vec3d & pos_a, LosTopos::Vec3d & pos_b);
    bool            solid_edge_is_feature(const LosTopos::SurfTrack & st, size_t e);
    
    // T1Transition::VelocityFieldCallback methods
    LosTopos::Vec3d sampleVelocity(LosTopos::Vec3d & pos);
    bool sampleDirectionalDivergence(const LosTopos::Vec3d & pos, const LosTopos::Vec3d & dir, double & output);

    
protected:
    bool m_verbose;
    
    std::string m_scene;
    std::string m_output_directory;
    std::string m_assets_directory;
    bool m_bbwall;
    bool m_t1vel;
    
    double m_dt;
    double m_time;
    int m_frameid;
    bool m_finished;
   
    LosTopos::SurfTrack * m_st;
    
};

#endif
