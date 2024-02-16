#include <iostream>
#include <fstream>
#include <filesystem>
#include <eigen3/Eigen/Dense>
#include "StGrid.h"
#include "StGridv2.h"
#include "FlowCpp.h"
#include "FlowCppv2.h"
#include <chrono>

using namespace std;

int main()
{
    // ================= Grid Parameters ====================
    const int COLS = 257; // Number of grid pts in the x direction
    const int ROWS = 257; // Number of grid pts in the y direction
    const double LENGTH = 4.; // Length of computational domain in the x direction
    const double BREADTH = 4.; // Breadth of computational domain in the y direction

    // ================= Fluid Parameters ====================
    const double CFL = 6.; // For Covector Fluid

    // ================= Simulation Parameters ===============
    const int sim_time = 150;
    const int interval = 10;

    // ================= File setup ==========================
    string NAME = "test1";
    string PATH = "../data/" + NAME;
    //filesystem::path test = filesystem::current_path();
    //cout << test << endl;
    filesystem::create_directory(PATH);

    // ================= Boundary Conditions =================
    // Basic Boundary conditions
    StGridv2::Boundary noslip = {"D", 0};
    StGridv2::Boundary zeroflux = {"N", 0};
    StGridv2::Boundary pressure = {"D", 0};
    StGridv2::Boundary freeflow = {"N", 0};

    // Lid-cavity problem setup
    double lidSpeed = 1;
    StGridv2::Boundary flow_lid = {"D", lidSpeed};
    StGridv2::Boundaries u_lid = {noslip, noslip, flow_lid, noslip};
    StGridv2::Boundaries v_lid = {noslip, noslip, noslip, noslip};
    StGridv2::Boundaries p_lid = {zeroflux, zeroflux, pressure, zeroflux};

    // Tunnel
    StGridv2::Boundary flow_tunnel = {"C", 0};
    StGridv2::Boundaries u_tunnel = {flow_tunnel, freeflow, noslip, noslip};
    StGridv2::Boundaries v_tunnel = {flow_tunnel, freeflow, noslip, noslip};
    StGridv2::Boundaries p_tunnel = {flow_tunnel, freeflow, pressure, pressure};

    // No conditions
    StGridv2::Boundary null = {"null", 0};
    StGridv2::Boundaries nulls = {null, null, null, null};

    // ================ General Simulation info =============
    cout << "============== Begining Simulation ================" << endl;
    cout << "===================================================" << endl;
    cout << "=== Simulation time : " << sim_time << " ===" << endl;

    // Lid-driven
    //StGridv2::StGridv2 grid = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, u_lid, v_lid, p_lid, CFL);

    // Tunnel
    StGridv2::StGridv2 grid = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, u_tunnel, v_tunnel, p_tunnel, CFL);
    Eigen::ArrayXd u1 = Eigen::ArrayXd::Ones(ROWS * (COLS + 1));
    grid.setU(u1);

    double t = 0.;
    int i = 0;
    int barWidth = 60;
    auto start = std::chrono::high_resolution_clock::now();

    while (t < sim_time){
        grid.SetTimeStep();
        double time_step = grid.getDt();
        FlowCppv2::BasicFlow(grid, time_step);
        //FlowCppv2::CovectorFluids1(grid, time_step);

        if (i % interval == 0){
            grid.WriteToFile(PATH + "/data_" + to_string(i) + ".vtk");
        }

        t += time_step;
        ++i;
        auto now = std::chrono::high_resolution_clock::now();
        double elapsed_time = chrono::duration_cast<chrono::seconds>(now - start).count();
        double time_left = elapsed_time * (sim_time / t - 1.0);
        double progress = t / sim_time;
        int pos = barWidth * progress;
        cout << "[";
        for (int n = 0; n < barWidth; ++n){
            if (n < pos) cout << "=";
            else if (n == pos) cout << ">";
            else cout << " ";
        }
        cout << "] Estimated Time Left : " << int(time_left / 60) << " min" << endl;
        cout.flush();
    }
    // Here goes the whole thing

    auto now = std::chrono::high_resolution_clock::now();
    double elapsed_time = chrono::duration_cast<chrono::seconds>(now - start).count();
    cout << "===================================================" << endl;
    cout << "================ Simulation Ended =================" << endl;
    cout << "Simulation ended in " << int(elapsed_time / 60) << " min" << endl;

//    StGridv2::StGridv2 grid = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, u_lid, v_lid, p_lid, CFL);
//    StGridv2::StGridv2 grid = StGridv2::StGridv2(11, 11, 1, 1, u_lid, v_lid, p_lid, CFL);
//    grid.SetTimeStep();
//    FlowCppv2::BasicFlow(grid, grid.getDt());
//    FlowCppv2::CovectorFluids1(grid, grid.getDt());
//    cout << "divergence : \n" << grid.ComputeDivergence().reshaped(grid.getRows(), grid.getCols()) << endl;
//    cout << "u : \n" << grid.getU().reshaped(grid.getRows(), grid.getCols() + 1) << endl << endl;
//    Eigen::ArrayXd u1 = Eigen::ArrayXd::Ones(ROWS * (COLS + 1));
//    grid.setU(u1);
//    cout << "test of GetVelocity on u : \n" << FlowCppv2::GetVelocityU(grid, grid.getU(), grid.getXs(), grid.getYs()).reshaped(11, 11) << endl;

//    grid.WriteToFile(PATH + "/test.vtk");
}

