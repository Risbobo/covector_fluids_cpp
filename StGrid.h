//
// Created by boris on 04.04.23.
//

#ifndef COVECTOR_FLUIDS_CPP_STGRID_H
#define COVECTOR_FLUIDS_CPP_STGRID_H
#include <eigen3/Eigen/Dense>

using Eigen::MatrixXd;

namespace StGrid {
    class StGrid {
        // Grid Dimension
        // along y
        int rows;
        // along x
        int cols;
        double length;
        double breadth;
        double dx;
        double dy;
        double dt;

        // Velocity Matrix
        MatrixXd u;
        MatrixXd v;

        // Inverse Flow Map
        MatrixXd psi_x;
        MatrixXd psi_y;

        // Pressure Matrix
        MatrixXd p;

        // Fluid Properties
        double CFL;
        double left_velocity;
        double rho;


    public:
        StGrid(int rows, int cols, double length, double breadth, double CFL, double left_velocity, double rho);

        void SetTimeStep();

        void ApplyVelocityBoundaries();
        void ApplyPressureBoundaries();

        // Pressure Projection
        void SolvePressurePoisson();
        void SolveMomentumEquation();

        int getRows() const;

        int getCols() const;

        double getDx() const;

        double getDy() const;

        double getDt() const;

        double getPsiX_ij(int i, int j) const;
        double getPsiY_ij(int i, int j) const;

        const MatrixXd &getU() const;

        const MatrixXd &getV() const;

        const MatrixXd &getP() const;

        void setPsiX_ij(int i, int j, double val);
        void setPsiY_ij(int i, int j, double val);

        void setU(const MatrixXd &u);
        void setV(const MatrixXd &v);

    };

} // StGrid

#endif //COVECTOR_FLUIDS_CPP_STGRID_H
