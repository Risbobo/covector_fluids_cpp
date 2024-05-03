//
// Created by boris on 19.01.24.
//

#include "FlowCppv2.h"

namespace FlowCppv2{
    // error tolerance for RG45
    double EPSILON = 0.001;

    // Bound the values of an array Ind between a l(ower)  and u(pper) bound
    Eigen::ArrayXi bounded(int l, int u, Eigen::ArrayXi Ind){
        return Ind.max(l).min(u);
    }

    // Linear interpolation
    Eigen::ArrayXd linint(const Eigen::ArrayXd& u0, const Eigen::ArrayXd& u1, const Eigen::ArrayXd& a){
        return (1 - a) * u0 + a * u1;
    }

    // Bilinear interpolation - Apply the horizontal linear interpolation after two vertical linear interpolation
    Eigen::ArrayXd bilinint(const Eigen::ArrayXd& u00, const Eigen::ArrayXd& u01, const Eigen::ArrayXd& u10, const Eigen::ArrayXd u11, const Eigen::ArrayXd a, const Eigen::ArrayXd b){
        return linint(linint(u00, u01, b), linint(u10, u11, b), a);
    }

    Eigen::ArrayXd GetVelocityU(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys){
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getCols();

        // for u (horizontal component)
        // offset for the staggered grid
        const Eigen::ArrayXd& uxs = xs;
        // With the offset on y, we can have the four u "around" the coord (x, y) -> u00, u01, u10, u11
        Eigen::ArrayXd uys = ys - 0.5 * dy;
        // Compute the indexes for each point
        // Offset of 1 because of the virtual row at j = 0
        Eigen::ArrayXi I_u = (uxs / dx).floor().cast<int>();
        Eigen::ArrayXi J_u = (uys / dy).floor().cast<int>() + 1;

        // Done check resize u
        // lower left
        Eigen::ArrayXd u00 = u_f(grid.indexFUVec(bounded(0, cols, I_u), bounded(0, rows + 1, J_u)));
        // upper left
        Eigen::ArrayXd u01 = u_f(grid.indexFUVec(bounded(0, cols, I_u), bounded(0, rows + 1, J_u + 1)));
        // lower right
        Eigen::ArrayXd u10 = u_f(grid.indexFUVec(bounded(0, cols, I_u + 1), bounded(0, rows + 1, J_u)));
        // upper right
        Eigen::ArrayXd u11 = u_f(grid.indexFUVec(bounded(0, cols, I_u + 1), bounded(0, rows + 1, J_u + 1)));

        // Bilinear interpolation factors
        // remove offset
        Eigen::ArrayXd A_u = (uxs / dx) - I_u.cast<double>();
        Eigen::ArrayXd B_u = (uys / dy) - (J_u.cast<double>() - 1);

        return bilinint(u00, u01, u10, u11, A_u, B_u);
    }

    Eigen::ArrayXd GetVelocityV(StGridv2::StGridv2 grid, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys){
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getCols();

        // for v (vertical component)
        Eigen::ArrayXd vxs = xs - 0.5 * dx;
        Eigen::ArrayXd vys = ys;
        Eigen::ArrayXi I_v = (vxs / dx).floor().cast<int>() + 1;
        Eigen::ArrayXi J_v = (vys / dy).floor().cast<int>();

        // DONE check v resize
        Eigen::ArrayXd v00 = v_f(grid.indexFVVec(bounded(0, cols + 1, I_v), bounded(0, rows, J_v)));
        Eigen::ArrayXd v01 = v_f(grid.indexFVVec(bounded(0, cols + 1, I_v), bounded(0, rows, J_v + 1)));
        Eigen::ArrayXd v10 = v_f(grid.indexFVVec(bounded(0, cols + 1, I_v + 1), bounded(0, rows, J_v)));
        Eigen::ArrayXd v11 = v_f(grid.indexFVVec(bounded(0, cols + 1, I_v + 1), bounded(0, rows, J_v + 1)));

        Eigen::ArrayXd A_v = (vxs / dx) - (I_v.cast<double>() - 1);
        Eigen::ArrayXd B_v = (vys / dy) - J_v.cast<double>();

        return bilinint(v00,  v01, v10, v11, A_v, B_v);
    }

    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    GetVelocity(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys) {
        return {GetVelocityU(grid, u_f, xs, ys), GetVelocityV(grid, v_f, xs, ys)};
    }

    // Using RK4
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    FlowMapPsi4(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys,
                double dt) {
        double c1 = 1.0 / 6.0 * dt;
        double c2 = 1.0 / 3.0 * dt;
        double c3 = c2;
        double c4 = c1;

        auto [k1u, k1v] = GetVelocity(grid, u_f, v_f, xs, ys);
        Eigen::ArrayXd mid1x = xs + 0.5 * dt * k1u;
        Eigen::ArrayXd mid1y = ys + 0.5 * dt * k1v;
        auto [k2u, k2v] = GetVelocity(grid, u_f, v_f, mid1x, mid1y);
        Eigen::ArrayXd mid2x = xs + 0.5 * dt * k2u;
        Eigen::ArrayXd mid2y = ys + 0.5 * dt * k2v;
        auto [k3u, k3v] = GetVelocity(grid, u_f, v_f, mid2x, mid2y);
        Eigen::ArrayXd mid3x = xs + dt * k3u;
        Eigen::ArrayXd mid3y = ys + dt * k3v;
        auto [k4u, k4v] = GetVelocity(grid, u_f, v_f, mid3x, mid3y);

        return {xs + c1 * k1u + c2 * k2u + c3 * k3u + c4 * k4u, ys + c1 * k1v + c2 * k2v + c3 * k3v + c4 * k4v};
    }

    // Using adaptative RK (45)
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd, double>
    FlowMapPsi45(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys,
                 double dt) {
        // COEFFICIENTS FOR Sarafyan's RK4(5), Table IV in Fehlberg
        double b21 = 1./2.;
        double b31 = 1./4.,     b32 = 1./4.;
        double b41 = 0.,        b42 = -1.,          b43 = 2.;
        double b51 = 7./27.,    b52 = 10./27.,      b53 = 0,        b54 = 1./27.;
        double b61 = 28./625.,  b62 = -1./5.,       b63 = 546./625.,b64 = 54./625., b65 = -378./625;
        double ch1 = 1./24.,    ch2 = 0.,           ch3 = 0.,       ch4 = 5./48.,   ch5 = 27./56.,  ch6 = 125./336;
        double ct1 = 1./8.,     ct2 = 0.,           ct3 = 2./3.,    ct4 = 1./16.,   ct5 = -27./56,  ct6 = -125./336;

        auto [k1u, k1v] = GetVelocity(grid, u_f, v_f, xs, ys);
        Eigen::ArrayXd mid1x = xs + b21 * dt * k1u;
        Eigen::ArrayXd mid1y = ys + b21 * dt * k1v;
        auto [k2u, k2v] = GetVelocity(grid, u_f, v_f, mid1x, mid1y);
        Eigen::ArrayXd mid2x = xs + b31 * dt * k1u + b32 * dt * k2u;
        Eigen::ArrayXd mid2y = ys + b31 * dt * k1v + b32 * dt * k2v;
        auto [k3u, k3v] = GetVelocity(grid, u_f, v_f, mid2x, mid2y);
        Eigen::ArrayXd mid3x = xs + b41 * dt * k1u + b42 * dt * k2u + b43 * dt * k3u;
        Eigen::ArrayXd mid3y = ys + b41 * dt * k1v + b42 * dt * k2v + b43 * dt * k3v;
        auto [k4u, k4v] = GetVelocity(grid, u_f, v_f, mid3x, mid3y);
        Eigen::ArrayXd mid4x = xs + b51 * dt * k1u + b52 * dt * k2u + b53 * dt * k3u + b54 * dt * k4u;
        Eigen::ArrayXd mid4y = ys + b51 * dt * k1v + b52 * dt * k2v + b53 * dt * k3v + b54 * dt * k4v;
        auto [k5u, k5v] = GetVelocity(grid, u_f, v_f, mid4x, mid4y);
        Eigen::ArrayXd mid5x = xs + b61 * dt * k1u + b62 * dt * k2u + b63 * dt * k3u + b64 * dt * k4u + b65 * dt * k5u;
        Eigen::ArrayXd mid5y = ys + b61 * dt * k1v + b62 * dt * k2v + b63 * dt * k3v + b64 * dt * k4v + b65 * dt * k5v;
        auto [k6u, k6v] = GetVelocity(grid, u_f, v_f, mid5x, mid5y);

        Eigen::ArrayXd nextx = xs + ch1 * dt * k1u + ch2 * dt * k2u + ch3 * dt * k3u + ch4 * dt * k4u + ch5 * dt * k5u + ch6 * dt * k6u;
        Eigen::ArrayXd nexty = ys + ch1 * dt * k1v + ch2 * dt * k2v + ch3 * dt * k3v + ch4 * dt * k4v + ch5 * dt * k5v + ch6 * dt * k6v;

        // Troncature error
        double te_x = (ct1 * dt * k1u + ct2 * dt * k2u + ct3 * dt * k3u + ct4 * dt * k4u + ct5 * dt * k5u + ct6 * dt * k6u).abs().maxCoeff();
        double te_y = (ct1 * dt * k1v + ct2 * dt * k2v + ct3 * dt * k3v + ct4 * dt * k4v + ct5 * dt * k5v + ct6 * dt * k6v).abs().maxCoeff();
        double te = std::max(te_x, te_y);

        return {nextx, nexty, te};
    }

    void BasicFlow(StGridv2::StGridv2 &grid) {
        grid.ApplyVelocityBoundaries();
        grid.ApplyVelocityObstacles();
        double dt = grid.getDt();
        // u_f <- u | Freeze Flow Velocity
        Eigen::ArrayXd u_f = grid.getU();
        Eigen::ArrayXd v_f = grid.getV();
        // u(x) <- u(psi(x))
        auto [nu, nv] = BasicAdvectionSL(grid, grid.getU(), grid.getV(), u_f, v_f, dt);
        grid.setU(nu);
        grid.setV(nv);

        // u <- P(u) | Pressure Projection
        //grid.SolvePressureLin();
        grid.SolvePressureIter();
        grid.SolveMomentumEquation();
    }

    // First order Covector Fluids (Algorithm 1)
    void CovectorFluids1(StGridv2::StGridv2 &grid) {
        grid.ApplyVelocityBoundaries();
        grid.ApplyVelocityObstacles();
        double dt = grid.getDt();
        // u_f <- u | Freeze Flow Velocity
        Eigen::ArrayXd u_f = grid.getU();
        Eigen::ArrayXd v_f = grid.getV();
        // u <- A_covec(u_f, u, dt) | Covector Lie Advection
        auto [nu, nv] = CovectorAdvectionSL(grid, grid.getU(), grid.getV(), u_f, v_f, dt);
        grid.setU(nu);
        grid.setV(nv);
        // u <- P(u) | Pressure Projection
        //grid.SolvePressureLin();
        grid.SolvePressureIter();
        grid.SolveMomentumEquation();
    }

    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    BasicAdvectionSL(StGridv2::StGridv2 &grid, Eigen::ArrayXd u, Eigen::ArrayXd v, Eigen::ArrayXd u_f,
                     Eigen::ArrayXd v_f, double dt) {
        Eigen::ArrayXd xs = grid.getXs();
        Eigen::ArrayXd ys = grid.getYs();
        int rows = grid.getRows();
        int cols = grid.getCols();
        Eigen::ArrayXd nu = u;
        Eigen::ArrayXd nv = v;
        // Compute the Flow Map at the center of every cell
        // Applying adaptive RG
        Eigen::ArrayXd psix, psiy;
        double ndt = dt;
        double te;
        do {
            auto [npsix, npsiy, nte] = FlowMapPsi45(grid, u_f, v_f, xs, ys, ndt);
            psix = npsix;
            psiy = npsiy;
            te = nte;

            if(te) {
                ndt = 0.9 * ndt * std::pow(EPSILON / te, 1. / 5.);
            }
        } while (te > EPSILON);
        grid.setDt(ndt);
        // Linear interpolation of the value of psi on the interior vertical faces (not on the boundary as it is applied later)
        // This compute psi(x)
        Eigen::ArrayXd psix_Fu = linint(psix.reshaped(rows, cols)(Eigen::all, Eigen::seq(0, Eigen::last - 1)).reshaped(),
                                        psix.reshaped(rows, cols)(Eigen::all, Eigen::seq(1, Eigen::last)).reshaped(),
                                        Eigen::ArrayXd::Ones(rows * (cols - 1)) * 0.5);
        Eigen::ArrayXd psiy_Fu = linint(psiy.reshaped(rows, cols)(Eigen::all, Eigen::seq(0, Eigen::last - 1)).reshaped(),
                                        psiy.reshaped(rows, cols)(Eigen::all, Eigen::seq(1, Eigen::last)).reshaped(),
                                        Eigen::ArrayXd::Ones(rows * (cols - 1)) * 0.5);
        Eigen::ArrayXd psix_Fv = linint(psix.reshaped(rows, cols)(Eigen::seq(0, Eigen::last - 1), Eigen::all).reshaped(),
                                        psix.reshaped(rows, cols)(Eigen::seq(1, Eigen::last), Eigen::all).reshaped(),
                                        Eigen::ArrayXd::Ones((rows - 1) * cols) * 0.5);
        Eigen::ArrayXd psiy_Fv = linint(psiy.reshaped(rows, cols)(Eigen::seq(0, Eigen::last - 1), Eigen::all).reshaped(),
                                        psiy.reshaped(rows, cols)(Eigen::seq(1, Eigen::last), Eigen::all).reshaped(),
                                        Eigen::ArrayXd::Ones((rows - 1) * cols) * 0.5);
        // Compute the new velocity on the interior vertical and horizontal faces for the new u and new v respectively.
        // this compute u(psi(x))
        // Update the central part of u (without the boundaries)
        nu.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() = GetVelocityU(grid, u_f, psix_Fu, psiy_Fu);
        nv.reshaped(rows + 1, cols + 2)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() = GetVelocityV(grid, v_f, psix_Fv, psiy_Fv);

        return {nu, nv};
    }

    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    CovectorAdvectionSL(StGridv2::StGridv2 &grid, Eigen::ArrayXd u, Eigen::ArrayXd v,
                        Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, double dt) {
        Eigen::ArrayXd xs = grid.getXs();
        Eigen::ArrayXd ys = grid.getYs();
        int rows = grid.getRows();
        int cols = grid.getCols();
        Eigen::ArrayXd nu = u;
        Eigen::ArrayXd nv = v;
        // Compute the Flow Map at the center of every cell
        // Applying adaptive RG
        Eigen::ArrayXd psix, psiy;
        double ndt = dt;
        double te;
        do {
            auto [npsix, npsiy, nte] = FlowMapPsi45(grid, u_f, v_f, xs, ys, ndt);
            psix = npsix;
            psiy = npsiy;
            te = nte;
            ndt = 0.9 * ndt * std::pow(EPSILON / te, 1./5.);
        } while (te > EPSILON);
        grid.setDt(ndt);
        // Linear interpolation of the value of psi on the interior vertical faces (not on the boundary as it is applied later)
        // This compute psi(x)
        Eigen::ArrayXd psix_Fu = linint(psix.reshaped(rows, cols)(Eigen::all, Eigen::seq(0, Eigen::last - 1)).reshaped(),
                                        psix.reshaped(rows, cols)(Eigen::all, Eigen::seq(1, Eigen::last)).reshaped(),
                                        Eigen::ArrayXd::Ones(rows * (cols - 1)) * 0.5);
        Eigen::ArrayXd psiy_Fu = linint(psiy.reshaped(rows, cols)(Eigen::all, Eigen::seq(0, Eigen::last - 1)).reshaped(),
                                        psiy.reshaped(rows, cols)(Eigen::all, Eigen::seq(1, Eigen::last)).reshaped(),
                                        Eigen::ArrayXd::Ones(rows * (cols - 1)) * 0.5);
        Eigen::ArrayXd psix_Fv = linint(psix.reshaped(rows, cols)(Eigen::seq(0, Eigen::last - 1), Eigen::all).reshaped(),
                                        psix.reshaped(rows, cols)(Eigen::seq(1, Eigen::last), Eigen::all).reshaped(),
                                        Eigen::ArrayXd::Ones((rows - 1) * cols) * 0.5);
        Eigen::ArrayXd psiy_Fv = linint(psiy.reshaped(rows, cols)(Eigen::seq(0, Eigen::last - 1), Eigen::all).reshaped(),
                                        psiy.reshaped(rows, cols)(Eigen::seq(1, Eigen::last), Eigen::all).reshaped(),
                                        Eigen::ArrayXd::Ones((rows - 1) * cols) * 0.5);
        // Compute the new velocity on the interior vertical and horizontal faces for the new u and new v respectively.
        // this compute (div psi(x)) * u(psi(x))
        // Update the central part of u (without the boundaries)
        nu.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() = (grid.getDiv_x() * psix.matrix()).array() * GetVelocityU(grid, u_f, psix_Fu, psiy_Fu);
        nv.reshaped(rows + 1, cols + 2)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() = (grid.getDiv_y() * psiy.matrix()).array() * GetVelocityV(grid, v_f, psix_Fv, psiy_Fv);

        return {nu, nv};
    }
}