//
// Created by boris on 11.01.24.
//

#ifndef COVECTOR_FLUIDS_CPP_STGRIDV2_H
#define COVECTOR_FLUIDS_CPP_STGRIDV2_H
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCholesky>
#include <eigen3/Eigen/IterativeLinearSolvers>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <unordered_set>
#include <set>
#include <cmath>

namespace StGridv2{
    struct Boundary {
        std::string type;
        double value;
    };

    struct Boundaries {
        Boundary left;
        Boundary right;
        Boundary up;
        Boundary down;
    };

    class StGridv2 {
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

        // Coordinates in the Grid
        Eigen::ArrayXd xs;
        Eigen::ArrayXd ys;

        // Velocity Matrix
        Eigen::ArrayXd u;
        Eigen::ArrayXd v;

        // Pressure Matrix
        Eigen::ArrayXd p;

        // Differential operator Matrices
        Eigen::SparseMatrix<double> laplacian;
        Eigen::SparseMatrix<double> div_u;
        Eigen::SparseMatrix<double> div_v;
        Eigen::SparseMatrix<double> div_cx;
        Eigen::SparseMatrix<double> div_cy;
        Eigen::SparseVector<double> boundary_pressure;

        // Fluid Properties
        double CFL;

        // Boundaries
        Boundaries u_bound;
        Boundaries v_bound;
        Boundaries p_bound;

        // Obstacles
        std::unordered_set<int> obstacles_c;
        std::unordered_set<int> obstacles_u;
        std::unordered_set<int> obstacles_v;

    public:
        StGridv2(int rows, int cols, double length, double breadth, Boundaries u_bound, Boundaries v_bound, Boundaries p_bound, double CFL);

        void SetTimeStep();

        void ApplyVelocityBoundaries();
        void ApplyVelocityObstacles();
        void ApplyPressureBoundaries();

        void GenerateOperators();

        void AddObstacle(std::unordered_set<std::tuple<int, int>> indices);
        void AddCircleObstacle(double cx, double cy, double r);
        void AddSquareObstacle(double cx, double cy, double l);

        // Pressure Projection
        void SolvePressureLin();
        void SolvePressureIter();
        void SolveMomentumEquation();

        Eigen::VectorXd ComputeDivergence();

        int getRows() const;

        int getCols() const;

        double getBreadth() const;

        double getLength() const;

        double getDx() const;

        double getDy() const;

        double getDt() const;

        Eigen::ArrayXd getXs();

        Eigen::ArrayXd getYs();

        Eigen::SparseMatrix<double> getDiv_x();

        Eigen::SparseMatrix<double> getDiv_y();

        const Eigen::ArrayXd &getU() const;

        const Eigen::ArrayXd &getV() const;

        const Eigen::ArrayXd &getP() const;

        void setDt(const double dt);
        void setU(const Eigen::ArrayXd &u);
        void setV(const Eigen::ArrayXd &v);

        void WriteToFile(std::string name);

        bool checkPvalid(int i, int j) const;
        int indexP(int i, int j);
        int indexU(int i, int j);
        int indexV(int i, int j);
        int indexCU(int i, int j) const;
        int indexCV(int i, int j);
        Eigen::ArrayXi indexFUVec(Eigen::ArrayXi I, Eigen::ArrayXi J);
        Eigen::ArrayXi indexFVVec(Eigen::ArrayXi I, Eigen::ArrayXi J);
    };
}


#endif //COVECTOR_FLUIDS_CPP_STGRIDV2_H
