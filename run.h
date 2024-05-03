//
// Created by boris on 22.02.24.
//

#ifndef COVECTOR_FLUIDS_CPP_RUN_H
#define COVECTOR_FLUIDS_CPP_RUN_H

#include "StGridv2.h"

namespace run {

    class run {
    public:
        static void Simulation(StGridv2::StGridv2 grid, double sim_time, int interval, std::string Name, std::string Mode);
        static void test1(std::string Name);
        static void test2(std::string Name);
    };

} // run

#endif //COVECTOR_FLUIDS_CPP_RUN_H
