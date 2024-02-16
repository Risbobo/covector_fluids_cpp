//
// Created by boris on 04.04.23.
//

#include "StGrid.h"
#include <iostream>
#include <utility>
#include "cmath"
#include <tuple>
//#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;

// /!\ Obsolete
namespace StGrid {
    // Constructor
    StGrid::StGrid(int rows, int cols, double length, double breadth, double CFL, double left_velocity, double rho) :
                   rows(rows), cols(cols), length(length), breadth(breadth),
                   CFL(CFL), left_velocity(left_velocity), rho(rho) {
        dx = length / cols;
        dy = breadth / rows;
        dt = CFL * (dx + dy);

        u = MatrixXd::Zero(rows, cols + 1);
        v = MatrixXd::Zero(rows + 1, cols);

        psi_x = MatrixXd::Zero(rows, cols);
        psi_y = MatrixXd::Zero(rows, cols);

        p = MatrixXd::Zero(rows, cols);
    }

    void StGrid::ApplyVelocityBoundaries() {
        // Nothing for u(0, all) ?
        // Left
        u(Eigen::all, 0) = MatrixXd::Constant(u.rows(), 1, left_velocity);

        // Top
        v(0, Eigen::all) = MatrixXd::Constant(1, v.cols(), 0);

        // Bottom
        v(Eigen::last, Eigen::all) = MatrixXd::Constant(1, v.cols(), 0);
    }

    void StGrid::ApplyPressureBoundaries() {
        // Top
        p(0, Eigen::all) = p(1, Eigen::all);

        // Bottom
        p(Eigen::last, Eigen::all) = p(Eigen::last - 1, Eigen::all);
    }

    void StGrid::SetTimeStep() {
        double factor = u.maxCoeff() / dx + v.maxCoeff() / dy;
        if (factor == 0) {
            dt = CFL * (dx + dy);
        } else {
            dt = CFL / factor;
        }
    }

    int StGrid::getRows() const {
        return rows;
    }

    int StGrid::getCols() const {
        return cols;
    }

    double StGrid::getDx() const {
        return dx;
    }

    double StGrid::getDy() const {
        return dy;
    }

    double StGrid::getDt() const {
        return dt;
    }

    double StGrid::getPsiX_ij(int i, int j) const {
        return this->psi_x(i,j);
    }

    double StGrid::getPsiY_ij(int i, int j) const {
        return this->psi_y(i,j);
    }

    const MatrixXd &StGrid::getU() const {
        return u;
    }

    const MatrixXd &StGrid::getV() const {
        return v;
    }

    const MatrixXd &StGrid::getP() const {
        return p;
    }

    void StGrid::setPsiX_ij(int i, int j, double val) {
        this->psi_x(i, j) = val;
    }

    void StGrid::setPsiY_ij(int i, int j, double val) {
        this->psi_y(i, j) = val;
    }

    void StGrid::setU(const MatrixXd &u) {
        StGrid::u = u;
    }

    void StGrid::setV(const MatrixXd &v) {
        StGrid::v = v;
    }

    void StGrid::SolvePressurePoisson() {
        double factor = 1.0 / (2.0 / (dx * dx) + 2.0 / (dy * dy));

        double error = 1.0;
        double tol = 1.0e-3;

        MatrixXd ustar1_x = (u(Eigen::seq(1, Eigen::last - 1), Eigen::seq(2, Eigen::last - 1)) - u(Eigen::seq(1, Eigen::last - 1), Eigen::seq(0, Eigen::last - 3))) / (2 * dx);
        MatrixXd vstar1_y = (v(Eigen::seq(2, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) - u(Eigen::seq(0, Eigen::last - 3), Eigen::seq(1, Eigen::last - 1))) / (2 * dy);

        int i = 0;
        while (error > tol){
            ++i;
            MatrixXd p_old = p;

            MatrixXd p2_xy = (p_old(Eigen::seq(1, Eigen::last - 1), Eigen::seq(2, Eigen::last)) + p_old(Eigen::seq(1, Eigen::last - 1), Eigen::seq(0, Eigen::last - 2))) / (dx * dx) + (p_old(Eigen::seq(2, Eigen::last), Eigen::seq(1, Eigen::last - 1)) + p_old(Eigen::seq(0, Eigen::last - 2), Eigen::seq(1, Eigen::last - 1))) / (dy *dy);
            p(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) = p2_xy * factor - (rho * factor / dt) * (ustar1_x + vstar1_y);

            error = (p - p_old).array().abs().maxCoeff();

            ApplyPressureBoundaries();

            if (i > 500){
                tol *= 10;
            }
        }
    }

    void StGrid::SolveMomentumEquation() {
        MatrixXd p1_x = (p(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last)) - p(Eigen::seq(1, Eigen::last - 1), Eigen::seq(0, Eigen::last - 1))) / dx;
        u(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) = u(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) - (dt / rho) * p1_x;

        MatrixXd p1_y = (p(Eigen::seq(1, Eigen::last), Eigen::seq(1, Eigen::last - 1)) - p(Eigen::seq(0, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1))) / dy;
        v(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) = v(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)) - (dt / rho) * p1_y;

        ApplyVelocityBoundaries();
    }
} // StGrid