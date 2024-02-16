//
// Created by boris on 19.01.24.
//

#ifndef COVECTOR_FLUIDS_CPP_FLOWCPPV2_H
#define COVECTOR_FLUIDS_CPP_FLOWCPPV2_H
#include <eigen3/Eigen/Dense>
#include "StGridv2.h"


namespace FlowCppv2{
    // Functions for Covector Fluid Solver
    Eigen::ArrayXd GetVelocityU(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys);
    Eigen::ArrayXd GetVelocityV(StGridv2::StGridv2 grid, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys);
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd> GetVelocity(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys);
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd> FlowMapPsi(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys, double dt);
    void BasicFlow(StGridv2::StGridv2& grid, double dt);
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd> CovectorAdvectionSL(StGridv2::StGridv2 grid, Eigen::ArrayXd u, Eigen::ArrayXd v, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, double dt);
    void CovectorFluids1(StGridv2::StGridv2& grid, double dt);
}


#endif //COVECTOR_FLUIDS_CPP_FLOWCPPV2_H
