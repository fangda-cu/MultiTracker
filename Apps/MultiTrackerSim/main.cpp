#include <iostream>

#include "Sim.h"

Sim g_Sim(true);

int main(int argc, const char * argv[])
{
    if (argc != 3 && argc != 4)
    {
        std::cout << "Usage:\n\tMultiTrackerSim option_file output_directory [assets_directory]\n\nThe assets directory is optional and defaults to \"./assets\".\n" << std::endl;
        return 0;
    } else
    {
        bool success = g_Sim.init(argv[1], argv[2], argc == 4 ? argv[3] : "./assets");
        if (!success)
            return 1;
    }
    
    std::cout << "Initialization complete. Starting the simulation..." << std::endl;
    
    while (!g_Sim.isFinished())
        g_Sim.step();
    
}

