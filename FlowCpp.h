//
// Created by boris on 20.06.23.
//

#ifndef COVECTOR_FLUIDS_CPP_FLOWCPP_H
#define COVECTOR_FLUIDS_CPP_FLOWCPP_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;

namespace FlowCpp{
    // Functions for Covector Fluid Solver
    Eigen::Vector2d GetVelocity(StGrid::StGrid grid, MatrixXd u_f, MatrixXd v_f, double x, double y);
    std::tuple<MatrixXd, MatrixXd> GetVelocityVec(const StGrid::StGrid& grid, MatrixXd u_f, MatrixXd v_f, Eigen::MatrixXd xs, Eigen::MatrixXd ys);
    Eigen::Vector2d FlowMapPsi(StGrid::StGrid grid, const MatrixXd& u_f, const MatrixXd& v_f, double x, double y, double dt);
    std::tuple<MatrixXd, MatrixXd> FlowMapPsiVec(StGrid::StGrid grid, const MatrixXd& u_f, const MatrixXd& v_f, Eigen::MatrixXd xs, Eigen::MatrixXd ys, double dt);
    std::tuple<MatrixXd , MatrixXd> CovectorAdvectionSL(StGrid::StGrid grid, const MatrixXd& u, const MatrixXd& v, const MatrixXd& u_f, const MatrixXd& v_f, double dt);
    std::tuple<MatrixXd, MatrixXd> CovectorAdvectionSLVec(StGrid::StGrid grid, const MatrixXd& u, const MatrixXd& v, const MatrixXd &u_f, const MatrixXd &v_f, double dt);
    std::tuple<MatrixXd , MatrixXd> CovectorAdvectionBFECC(StGrid::StGrid grid, const MatrixXd& u, const MatrixXd& v, const MatrixXd& u_f, const MatrixXd& v_f, double dt);

    void CovectorFluid1(StGrid::StGrid grid);
    void CovectorFluid2(StGrid::StGrid grid);
}

#endif //COVECTOR_FLUIDS_CPP_FLOWCPP_H
