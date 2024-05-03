//
// Created by boris on 22.02.24.
//

#include "run.h"
#include "FlowCppv2.h"
#include <chrono>

namespace run {
    // Simulation applied to the grid with the Mode (Basic or Covector Fluids)
    void run::Simulation(StGridv2::StGridv2 grid, double sim_time, int interval, std::string Name, std::string Mode) {
        // ================ General Simulation info =============
        std::cout << "============== Begining Simulation ================" << std::endl;
        std::cout << "===================================================" << std::endl;
        std::cout << "=== Simulation time : " << sim_time << " ===" << std::endl;
        double max_V = std::max(grid.getU().maxCoeff(), grid.getV().maxCoeff());
        //std::cout << "=== Re : " << max_V * std::max(grid.getBreadth(), grid.getLength()) << " ===" << std::endl;
        double t = 0.;
        int i = 0;
        int barWidth = 60;
        auto start = std::chrono::high_resolution_clock::now();

        while (t < sim_time) {
            //grid.SetTimeStep();
            if (Mode == "Basic") {
                FlowCppv2::BasicFlow(grid);
            } else if (Mode == "CovectorFluids1") {
                FlowCppv2::CovectorFluids1(grid);
            } else {
                std::cerr << "Wrong mode" << std::endl;
                return;
            }
            double time_step = grid.getDt();
            if (i % interval == 0) {
                grid.WriteToFile(Name + "/data_" + std::to_string(i) + ".vtk");
            }

            t += time_step;
            ++i;
            auto now = std::chrono::high_resolution_clock::now();
            double elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
            double time_left = elapsed_time * (sim_time / t - 1.0);
            double progress = t / sim_time;
            int pos = barWidth * progress;
            std::cout << "[";
            for (int n = 0; n < barWidth; ++n) {
                if (n < pos) std::cout << "=";
                else if (n == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout.precision(3);
            std::cout << "] Estimated Time Left : " << int(time_left / 60) << " min (dt : " << grid.getDt() << ") | max div : " << grid.ComputeDivergence().cwiseAbs().maxCoeff() << std::endl;
            std::cout.flush();
        }

        auto now = std::chrono::high_resolution_clock::now();
        double elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(now - start).count();
        std::cout << "===================================================" << std::endl;
        std::cout << "================ Simulation Ended =================" << std::endl;
        std::cout << "Simulation ended in " << int(elapsed_time / 60) << " min" << std::endl;
    }

    // test for the ray tracing subroutine (FlowMap function)
    void run::test1(std::string Name) {
        std::cout << "============== Begining Test 1 ================" << std::endl;
        std::cout << "================ Ray Tracing ==================" << std::endl;
        StGridv2::Boundary null = {"null", 0};
        StGridv2::Boundaries nulls = {null, null, null, null};
        int COLS = 257;
        int ROWS = 257;
        double LENGTH = 4.0;
        double BREADTH = 4.0;
        int N = 40;

        StGridv2::StGridv2 grid = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, nulls, nulls, nulls, 6.0);
        double cx = 2.;
        double cy = 2.;
        Eigen::ArrayXd ui = Eigen::ArrayXd::Zero((ROWS + 2) * (COLS + 1));
        Eigen::ArrayXd vi = Eigen::ArrayXd::Zero((ROWS + 1) * (COLS + 2));
        // V(x, y) = (-y, x)
        for (int i = 0; i < COLS + 1; ++i) {
            for (int j = 0; j < ROWS + 2; ++j) {
                double x = i * grid.getDx() - cx;
                double y = (j - 0.5) * grid.getDy() - cy;
                double dist = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
                if (dist <= 1.5) {
                    ui.reshaped(ROWS + 2, COLS + 1)(j, i) = -y / dist;
                }
            }
        }
        for (int i = 0; i < COLS + 2; ++i) {
            for (int j = 0; j < ROWS + 1; ++j) {
                double x = (i - 0.5) * grid.getDx() - cx;
                double y = j * grid.getDy() - cy;
                double dist = std::sqrt(std::pow(x, 2) + std::pow(y, 2));
                if (dist <= 1.5) {
                    vi.reshaped(ROWS + 1, COLS + 2)(j, i) = x / dist;
                }
            }
        }
        grid.setU(ui);
        grid.setV(vi);
        Eigen::ArrayXd particles_x = Eigen::ArrayXd::LinSpaced(N, -1.4, 1.4) + 2.0;
        Eigen::ArrayXd particles_y = Eigen::ArrayXd::Zero(N) + 2.0;

        std::ofstream myfile;
        myfile.open(Name + "/particle_" + std::to_string(0) + ".vtk");
        myfile << "# vtk DataFile Version 2.0" << std::endl;
        myfile << "First test" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET STRUCTURED_GRID" << std::endl;
        myfile << "DIMENSIONS " << N << " " << 1 << " 1" << std::endl;

        myfile << "POINTS " << N << " double" << std::endl;
        for(int j = 0; j < N; ++j){
            myfile << particles_x(j) << " " << particles_y(j) << " " << 0 << std::endl;
        }

        myfile << "POINT_DATA " << N << std::endl;
        myfile << "VECTORS velocity double" << std::endl;
        auto [u, v] = FlowCppv2::GetVelocity(grid, grid.getU(), grid.getV(), particles_x, particles_y);
        for (int j = 0; j < N; ++j) {
            myfile << u(j) << " " << v(j) << " " << 0.0 << std::endl;
        }
        myfile.close();
        grid.WriteToFile(Name + "/grid_" + std::to_string(0) + ".vtk");

        int barWidth = 60;
        double EPSILON = 0.001;
        std::cout << "Test starting..." << std::endl;
        for (int i = 1; i < 120; ++i) {
            grid.SetTimeStep();
            double dt = grid.getDt();

            //auto [nx, ny] = FlowCppv2::FlowMapPsi4(grid, grid.getU(), grid.getV(), particles_x, particles_y, dt);
            Eigen::ArrayXd psix, psiy;
            double ndt = dt;
            double te;
            do {
                auto [npsix, npsiy, nte] = FlowCppv2::FlowMapPsi45(grid, grid.getU(), grid.getV(), particles_x, particles_y, ndt);
                psix = npsix;
                psiy = npsiy;
                te = nte;
                ndt = 0.9 * ndt * std::pow(EPSILON / te, 1./5.);
            } while (te > EPSILON);
            grid.setDt(ndt);
            particles_x = psix;
            particles_y = psiy;

            std::ofstream myfile;
            myfile.open(Name + "/particle_" + std::to_string(i) + ".vtk");
            myfile << "# vtk DataFile Version 2.0" << std::endl;
            myfile << "First test" << std::endl;
            myfile << "ASCII" << std::endl;
            myfile << "DATASET STRUCTURED_GRID" << std::endl;
            myfile << "DIMENSIONS " << N << " " << 1 << " 1" << std::endl;

            myfile << "POINTS " << N << " double" << std::endl;
            for(int j = 0; j < N; ++j){
                myfile << particles_x(j) << " " << particles_y(j) << " " << 0 << std::endl;
            }

            myfile << "POINT_DATA " << N << std::endl;
            myfile << "VECTORS velocity double" << std::endl;
            auto [u, v] = FlowCppv2::GetVelocity(grid, grid.getU(), grid.getV(), particles_x, particles_y);
            for (int j = 0; j < N; ++j) {
                myfile << u(j) << " " << v(j) << " " << 0.0 << std::endl;
            }
            myfile.close();
            grid.WriteToFile(Name + "/grid_" + std::to_string(i) + ".vtk");

            int pos = barWidth * i / 120;
            std::cout << "[";
            for (int n = 0; n < barWidth; ++n) {
                if (n < pos) std::cout << "=";
                else if (n == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "] - dt : " << grid.getDt() << std::endl;
            std::cout.flush();
        }
    }

    // test for the pressure solvers
    void run::test2(std::string Name) {
        std::cout << "============== Begining Test 2 ================" << std::endl;
        std::cout << "===============================================" << std::endl;
        StGridv2::Boundary null = {"null", 0};
        StGridv2::Boundaries nulls = {null, null, null, null};
        int COLS = 257;
        int ROWS = 257;
        double LENGTH = 4.0;
        double BREADTH = 4.0;

        StGridv2::StGridv2 grid = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, nulls, nulls, nulls, 6.0);
        grid.GenerateOperators();
        Eigen::ArrayXd vi = Eigen::ArrayXd::Zero((ROWS + 1) * (COLS + 2));
        vi.reshaped(ROWS + 1, COLS + 2)(int(ROWS / 2), int(COLS / 2) + 1) = 1.;
        grid.setV(vi);

        int barWidth = 60;
        std::cout << "Test starting..." << std::endl;
        for (int i = 0; i < 60; ++i) {
            //grid.SolvePressureLin();
            grid.SolvePressureIter();
            grid.SolveMomentumEquation();
            grid.WriteToFile(Name + "/test_" + std::to_string(i) + ".vtk");

            int pos = barWidth * i / 60;
            std::cout << "[";
            for (int n = 0; n < barWidth; ++n) {
                if (n < pos) std::cout << "=";
                else if (n == pos) std::cout << ">";
                else std::cout << " ";
            }
            std::cout << "]";
            std::cout << " Max div : " << grid.ComputeDivergence().maxCoeff() << std::endl;
            std::cout.flush();
        }
    }
} // run