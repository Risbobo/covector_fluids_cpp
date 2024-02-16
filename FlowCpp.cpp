//
// Created by boris on 20.06.23.
//

#include "StGrid.h"
#include "FlowCpp.h"
#include <iostream>
#include <utility>
#include "cmath"
#include <tuple>
//#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;


// /!\ Obsolete
namespace FlowCpp{
    Eigen::Vector2d GetVelocity(StGrid::StGrid grid, MatrixXd u_f, MatrixXd v_f, double x, double y) {
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getCols();
        // for u
        double upos_x = x;
        double upos_y = y - 0.5*dy;
        int i_u = floor(upos_x / dx);
        i_u = std::min(std::max(i_u, 0), rows - 2);
        int j_u = floor(upos_y / dy);
        j_u = std::min(std::max(j_u, 0), cols - 1);
        double u00 = u_f(i_u, j_u);
        double u01 = u_f(i_u, j_u + 1);
        double u10 = u_f(i_u + 1, j_u);
        double u11 = u_f(i_u + 1, j_u + 1);
        double a_u = upos_x / dx - i_u;
        double b_u = upos_y / dy - j_u;
        double u_xy = (1 - a_u) * (1 - b_u) * u00 + (1 - a_u) * b_u * u01 + a_u * (1 - b_u) * u10 + a_u * b_u * u11;

        // for v
        double vpos_x = x - 0.5*dx;
        double vpos_y = y;
        int i_v = floor(vpos_x / dx);
        i_v = std::min(std::max(i_v, 0), rows - 1);
        int j_v = floor(vpos_y / dy);
        j_v = std::min(std::max(j_v, 0), cols - 2);
        double v00 = v_f(i_v, j_v);
        double v01 = v_f(i_v, j_v + 1);
        double v10 = v_f(i_v + 1, j_v);
        double v11 = v_f(i_v + 1, j_v + 1);
        double a_v = vpos_x / dx - i_v;
        double b_v = vpos_y / dy - j_v;

        double v_xy = (1 - a_v) * (1 - b_v) * v00 + (1 - a_v) * b_v * v01 + a_v * (1 - b_v) * v10 + a_v * b_v * v11;

        return (Eigen::Vector2d(2) <<  u_xy, v_xy).finished();
    }

    std::tuple<MatrixXd, MatrixXd> GetVelocityVec(const StGrid::StGrid& grid, MatrixXd u_f, MatrixXd v_f, Eigen::MatrixXd xs, Eigen::MatrixXd ys){
        struct index{
            Eigen::Index size() const {return X.size();}
            Eigen::Index operator[] (Eigen::Index i) const {return (X(i) + X.rows() * Y(i))%X.size();}
            Eigen::ArrayXXi X, Y;
        };

        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getCols();

        // for u
        MatrixXd uxs = xs;
        MatrixXd uys = ys.array() - 0.5 * dy;
        // for every x in xs, we compute the corresponding index in the grid between 0 and rows - 2
        Eigen::ArrayXXi I_u = (uxs.array() / dx).floor().max(0).min(rows - 2).cast<int>();
        Eigen::ArrayXXi J_u = (uys.array() / dy).floor().max(0).min(cols - 1).cast<int>();
        Eigen::ArrayXXd u00 = u_f.reshaped()(index{I_u,     J_u    }).reshaped(u_f.rows(), u_f.cols());
        Eigen::ArrayXXd u01 = u_f.reshaped()(index{I_u,     J_u + 1}).reshaped(u_f.rows(), u_f.cols());
        Eigen::ArrayXXd u10 = u_f.reshaped()(index{I_u + 1, J_u    }).reshaped(u_f.rows(), u_f.cols());
        Eigen::ArrayXXd u11 = u_f.reshaped()(index{I_u + 1, J_u + 1}).reshaped(u_f.rows(), u_f.cols());
        Eigen::ArrayXXd A_u = (uxs / dx) - I_u.matrix().cast<double>();
        Eigen::ArrayXXd B_u = (uys / dy) - J_u.matrix().cast<double>();
        MatrixXd velocity_u = (1 - A_u) * (1 - B_u) * u00 + (1 - A_u) * B_u * u01 + A_u * (1 - B_u) * u10 + A_u * B_u * u11;

        // for v
        MatrixXd vxs = xs.array() - 0.5 * dy;
        MatrixXd vys = ys;
        Eigen::ArrayXXi I_v = (vxs.array() / dx).floor().max(0).min(rows - 1).cast<int>();
        Eigen::ArrayXXi J_v = (vys.array() / dy).floor().max(0).min(cols - 2).cast<int>();
        Eigen::ArrayXXd v00 = v_f.reshaped()(index{I_v,     J_v    }).reshaped(v_f.rows(), v_f.cols());
        Eigen::ArrayXXd v01 = v_f.reshaped()(index{I_v,     J_v + 1}).reshaped(v_f.rows(), v_f.cols());
        Eigen::ArrayXXd v10 = v_f.reshaped()(index{I_v + 1, J_v    }).reshaped(v_f.rows(), v_f.cols());
        Eigen::ArrayXXd v11 = v_f.reshaped()(index{I_v + 1, J_v + 1}).reshaped(v_f.rows(), v_f.cols());
        Eigen::ArrayXXd A_v = (vxs / dx) - I_v.matrix().cast<double>();
        Eigen::ArrayXXd B_v = (vys / dy) - J_v.matrix().cast<double>();
        MatrixXd velocity_v = (1 - A_v) * (1 - B_v) * v00 + (1 - A_v) * B_v * v01 + A_v * (1 - B_v) * v10 + A_v * B_v * v11;

        return {velocity_u, velocity_v};
    }

    Eigen::Vector2d FlowMapPsi(StGrid::StGrid grid, const MatrixXd& u_f, const MatrixXd& v_f, double x, double y, double dt) {
        double c1 = 1.0 / 6.0 * dt;
        double c2 = 1.0 / 3.0 * dt;
        double c3 = c2;
        double c4 = c1;
        Eigen::Vector2d pos(x, y);
        Eigen::Vector2d k1 = GetVelocity(grid, u_f, v_f, x, y);
        Eigen::Vector2d mid1 = pos + dt * k1 * 0.5;
        Eigen::Vector2d k2 = GetVelocity(grid, u_f, v_f, mid1(0), mid1(1));
        Eigen::Vector2d mid2 = pos + dt * k2 * 0.5;
        Eigen::Vector2d k3 = GetVelocity(grid, u_f, v_f, mid2(0), mid2(1));
        Eigen::Vector2d mid3 = pos + dt * k3;
        Eigen::Vector2d k4 = GetVelocity(grid, u_f, v_f, mid3(0), mid3(1));
        return c1 * k1 + c2 * k2 + c3 * k3 + c4 * k4;
    }

    std::tuple<MatrixXd, MatrixXd> FlowMapPsiVec(StGrid::StGrid grid, const MatrixXd& u_f, const MatrixXd& v_f, Eigen::MatrixXd xs, Eigen::MatrixXd ys, double dt){
        double c1 = 1.0 / 6.0 * dt;
        double c2 = 1.0 / 3.0 * dt;
        double c3 = c2;
        double c4 = c1;

        auto [k1u, k1v] = GetVelocityVec(grid, u_f, v_f, xs, ys);
        MatrixXd mid1x = xs + 0.5 * dt * k1u;
        MatrixXd mid1y = ys + 0.5 * dt * k1v;
        auto [k2u, k2v] = GetVelocityVec(grid, u_f, v_f, mid1x, mid1y);
        MatrixXd mid2x = xs + 0.5 * dt * k2u;
        MatrixXd mid2y = ys + 0.5 * dt * k2v;
        auto [k3u, k3v] = GetVelocityVec(grid, u_f, v_f, mid2x, mid2y);
        MatrixXd mid3x = xs + dt * k3u;
        MatrixXd mid3y = ys + dt * k3v;
        auto [k4u, k4v] = GetVelocityVec(grid, u_f, v_f, mid3x, mid3y);

        return {c1 * k1u + c2 * k2u + c3 * k3u + c4 * k4u, c1 * k1v + c2 * k2v + c3 * k3v + c4 * k4v};
    }

    std::tuple<MatrixXd, MatrixXd> CovectorAdvectionSL(StGrid::StGrid grid, const MatrixXd& u, const MatrixXd& v, const MatrixXd &u_f, const MatrixXd &v_f, double dt) {
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getRows();
        MatrixXd u_next = u; //MatrixXd::Zero(rows, cols + 1);
        MatrixXd v_next = v; //MatrixXd::Zero(rows + 1, cols);
        for(int i = 1; i < rows; ++i){
            for(int j = 1; j < cols; ++j){
                //std::cout << " coord : " << i << " " << j << std::endl;
                double x = i * dx;
                double y = j * dy;
                Eigen::Vector2d psi_center = FlowMapPsi(grid, u_f, v_f, x + 0.5 * dx, y + 0.5 * dy, dt);
                grid.setPsiX_ij(i, j, psi_center(0));
                grid.setPsiY_ij(i, j, psi_center(1));
                if (i > 0 && j > 0){
                    // Compute dPsi for u
                    double dPsi_ux = 1.0 / dx * (grid.getPsiX_ij(i, j) - grid.getPsiX_ij(i - 1, j));
                    double dPsi_uy = 1.0 / dx * (grid.getPsiY_ij(i, j) - grid.getPsiY_ij(i - 1, j));
                    // Compute u(Psi(x)) where x is the center of the face of u
                    Eigen::Vector2d psiFu = FlowMapPsi(grid, u_f, v_f, x, y -0.5 * dy, dt);
                    Eigen::Vector2d uF = GetVelocity(grid, u_f, v_f, psiFu(0), psiFu(1));
                    // Compute u <- dPsi^t * u(Psi(x))
                    u_next(i, j) = dPsi_ux * uF(0) + dPsi_uy * uF(1);

                    // Compute dPsi for v
                    double dPsi_vx = 1.0 / dy * (grid.getPsiX_ij(i, j) - grid.getPsiX_ij(i, j - 1));
                    double dPsi_vy = 1.0 / dy * (grid.getPsiY_ij(i, j) - grid.getPsiY_ij(i, j - 1));
                    // Compute u(Psi(x)) where x is the center of the face of u
                    Eigen::Vector2d psiFv = FlowMapPsi(grid, u_f, v_f, x - 0.5 * dx, y, dt);
                    Eigen::Vector2d vF = GetVelocity(grid, u_f, v_f, psiFv(0), psiFv(1));
                    // Compute u <- dPsi^t * u(Psi(x))
                    v_next(i, j) = dPsi_vx * vF(0) + dPsi_vy * vF(1);
                }
            }
        }
        return {u_next, v_next};
    }

    std::tuple<MatrixXd, MatrixXd> CovectorAdvectionSLVec(StGrid::StGrid grid, const MatrixXd& u, const MatrixXd& v, const MatrixXd &u_f, const MatrixXd &v_f, double dt){
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getRows();
        MatrixXd u_next = u;
        MatrixXd v_next = v;

        // TODO : We can generate xs, ys, xms and yms in constructor rather than in every time step
        // Generate Xs coordinates and Ys coordinates
        Eigen::ArrayXXd xs = Eigen::ArrayXd::LinSpaced(rows, 0, rows - 1).rowwise().replicate(cols) * dx;
        Eigen::ArrayXXd ys = Eigen::ArrayXd::LinSpaced(cols, 0, cols - 1).transpose().colwise().replicate(rows) * dy;
        // Generate middle coordinates (for face coord and center of cell coord)
        Eigen::ArrayXXd xms = xs + 0.5 * dx;
        Eigen::ArrayXXd yms = ys + 0.5 * dy;
        // Compute psi at center of cells
        auto [psi_mx, psi_my] = FlowMapPsiVec(grid, u_f, v_f, xms, yms, dt);
        //std::cout << "What is psi_mx here : \n" << psi_mx << std::endl;
        // Compute psi on the vertical faces
        auto [psi_fux, psi_fuy] = FlowMapPsiVec(grid, u_f, v_f, xs, yms, dt);
        // Compute psi on the horizontal faces
        auto [psi_fvx, psi_fvy] = FlowMapPsiVec(grid, u_f, v_f, xms, ys, dt);

        // Compute u <- dPsi^t * u(Psi(x))
        // dPsi_u = 1/dx * (psi_m[1:last, all] - psi_m[0:last-1, all]
        Eigen::ArrayXXd dPsi_ux = 1.0 / dx * (psi_mx(Eigen::all, Eigen::seq(1, Eigen::last)) - psi_mx(Eigen::all, Eigen::seq(0, Eigen::last - 1))).array();
        Eigen::ArrayXXd dPsi_uy = 1.0 / dx * (psi_my(Eigen::all, Eigen::seq(1, Eigen::last)) - psi_my(Eigen::all,Eigen::seq(0, Eigen::last - 1))).array();
        auto [u_psix, u_psiy] = GetVelocityVec(grid, u_f, v_f, psi_fux, psi_fuy);
        //std::cout << "What is u_next here : \n" << u_next << std::endl;
        //TODO : Vérifier les dimension ici, j'ai peur qu'il y ait toujours un problème
        u_next(Eigen::seq(0, Eigen::last - 1), Eigen::seq(1, Eigen::last))
                = (dPsi_ux * u_psix(Eigen::seq(1, Eigen::last), Eigen::all).array()).matrix()
                + (dPsi_uy * u_psiy(Eigen::seq(1, Eigen::last), Eigen::all).array()).matrix();
        //std::cout << "What is u_next here : \n" << u_next << std::endl;

        // Compute v <- dPsi^t * v(Psi(x))
        Eigen::ArrayXXd dPsi_vx = 1.0 / dy * (psi_mx(Eigen::seq(1, Eigen::last), Eigen::all) - psi_mx(Eigen::seq(0, Eigen::last - 1), Eigen::all)).array();
        Eigen::ArrayXXd dPsi_vy = 1.0 / dy * (psi_my(Eigen::seq(1, Eigen::last), Eigen::all) - psi_my(Eigen::seq(0, Eigen::last - 1), Eigen::all)).array();
        auto [v_psix, v_psiy] = GetVelocityVec(grid, u_f, v_f, psi_fvx, psi_fvy);
        v_next(Eigen::seq(1, Eigen::last), Eigen::seq(0, Eigen::last - 1))
                = (dPsi_vx * v_psix(Eigen::seq(1, Eigen::last), Eigen::all).array()).matrix()
                + (dPsi_vy * v_psiy(Eigen::seq(1, Eigen::last), Eigen::all).array()).matrix();

        return {u_next, v_next};
    }

    std::tuple<MatrixXd, MatrixXd> CovectorAdvectionBFECC(StGrid::StGrid grid, const  MatrixXd& u, const MatrixXd& v, const MatrixXd& u_f, const MatrixXd& v_f, double dt) {
        auto [u_1, v_1] = CovectorAdvectionSL(grid, u, v, u_f, v_f, dt);
        auto [u_0, v_0] = CovectorAdvectionSL(grid,u_1, v_1, u_f, v_f, -dt);

        auto e_u = u_0 - u_1;
        auto e_v = v_0 - v_1;
        return CovectorAdvectionSL(grid, e_u, e_v, u_f, v_f, dt);
    }

    void CovectorFluid1(StGrid::StGrid grid) {
        auto u = grid.getU();
        auto v = grid.getV();
        double dt = grid.getDt();
        auto [u_star, v_star] = CovectorAdvectionSLVec(grid, u, v, u, v, dt);
        grid.setU(u_star);
        grid.setV(v_star);
        // grid.ApplyObstacles();
        grid.SolvePressurePoisson();
        grid.SolveMomentumEquation();
        //grid.ApplyVelocityBoundaries();
    }
}