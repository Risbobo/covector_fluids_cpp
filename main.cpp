#include <iostream>
#include <fstream>
#include <filesystem>
#include <eigen3/Eigen/Dense>
#include "StGridv2.h"
#include "FlowCppv2.h"
#include "run.h"

using namespace std;

int main()
{
    // ================= Grid Parameters ====================
    const int COLS = 257; // Number of grid pts in the x direction (for ex: 257)
    const int ROWS = 257; // Number of grid pts in the y direction
    const double LENGTH = 4.; // Length of computational domain in the x direction
    const double BREADTH = 4.; // Breadth of computational domain in the y direction

    // ================= Simulation Modes ====================

    std::string BASIC = "Basic";
    std::string COVFLUIDS1 = "CovectorFluids1";

    // ================= Fluid Parameters ====================
    const double CFL = 6.; // For Covector Fluid

    // ================= Simulation Parameters ===============
    const int sim_time = 20;
    const int interval = 10;

    // ================= File setup ==========================
    string NAME = "test";
    string PATH = "../data/" + NAME;
    //filesystem::path test = filesystem::current_path();
    //cout << test << endl;
    filesystem::create_directory(PATH);

    // ================= Boundary Conditions =================
    // Basic Boundary conditions
    StGridv2::Boundary noslip = {"D", 0};
    StGridv2::Boundary nopressure = {"D", 0};
    StGridv2::Boundary zeroflux = {"N", 0};
    StGridv2::Boundary pressure = {"D", 1};
    StGridv2::Boundary freeflow = {"N", 0};

    // Lid-cavity problem setup
    double lidSpeed = 1;
    StGridv2::Boundary flow_lid = {"D", lidSpeed};
    StGridv2::Boundaries u_lid = {noslip, noslip, flow_lid, noslip};
    StGridv2::Boundaries v_lid = {noslip, noslip, noslip, noslip};
    StGridv2::Boundaries p_lid = {zeroflux, zeroflux, zeroflux, zeroflux};

    // Tunnel
    double entry_speed = 20;
    StGridv2::Boundary flow_tunnel = {"D", entry_speed};
    StGridv2::Boundaries u_tunnel = {flow_tunnel, freeflow, freeflow, freeflow};
    StGridv2::Boundaries v_tunnel = {freeflow, freeflow, freeflow, freeflow};
    StGridv2::Boundaries p_tunnel = {zeroflux, nopressure, zeroflux, zeroflux};

    // Poiseuille
    double entry_pressure = 1;
    StGridv2::Boundary flow_poiseuille = {"D", entry_pressure};
    StGridv2::Boundaries u_poiseuille = {freeflow, freeflow, noslip, noslip};
    StGridv2::Boundaries v_poiseuille = {freeflow, freeflow, noslip, noslip};
    StGridv2::Boundaries p_poiseuille = {flow_poiseuille, nopressure, zeroflux, zeroflux};

    // No conditions
    StGridv2::Boundary null = {"null", 0};
    StGridv2::Boundaries nulls = {null, null, null, null};

    // Lid-driven
    StGridv2::StGridv2 grid_lid = StGridv2::StGridv2(11, 11, LENGTH, BREADTH, u_lid, v_lid, p_lid, CFL);
    //grid_lid.GenerateOperators();
    // Small perturbation so it can start
    Eigen::ArrayXd ui = grid_lid.getU();
    ui.reshaped(grid_lid.getRows() + 2, grid_lid.getCols() + 1)(Eigen::last - 1, Eigen::all) = Eigen::ArrayXd::Zero(grid_lid.getRows() + 2) + 0.1;
    grid_lid.setU(ui);

    // Tunnel
    StGridv2::StGridv2 grid_tunnel = StGridv2::StGridv2(ROWS, COLS, LENGTH, BREADTH, u_tunnel, v_tunnel, p_tunnel, CFL);
    grid_tunnel.AddCircleObstacle(LENGTH / 3.0, BREADTH / 2.0, BREADTH / 10.0);
    //grid_tunnel.AddSquareObstacle(LENGTH / 3.0, BREADTH / 2.0, LENGTH / 16.0);
    grid_tunnel.GenerateOperators();

    // Poiseuille
    StGridv2::StGridv2 grid_poiseuille =StGridv2::StGridv2(11, 11, LENGTH, BREADTH, u_poiseuille, v_poiseuille, p_poiseuille, CFL);
    //grid_poiseuille.GenerateOperators();

    //run::run::Simulation(grid_lid, sim_time, interval, PATH, BASIC);
    //FlowCppv2::BasicFlow(grid_lid);

    run::run::Simulation(grid_tunnel, sim_time, interval, PATH, BASIC);
    //FlowCppv2::BasicFlow(grid_tunnel);

    //run::run::Simulation(grid_poiseuille, sim_time, interval, PATH, BASIC);
    //FlowCppv2::BasicFlow(grid_poiseuille);


    //run::run::test1(PATH);
    //run::run::test2(PATH);
}

