//
// Created by boris on 19.01.24.
//

#include "FlowCppv2.h"

namespace FlowCppv2{
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
        //std::cout << "u : \n" << u_f.reshaped(grid.getRows(), grid.getCols()) << std::endl;

        // for u (horizontal component)
        // offset for the staggered grid
        Eigen::ArrayXd uxs = xs;
        // With the offset on y, we can have the four u "around" the coord (x, y) -> u00, u01, u10, u11
        Eigen::ArrayXd uys = ys - 0.5 * dy;
        // Compute the indexes for each point
        Eigen::ArrayXi I_u = (uxs / dx).floor().cast<int>();
        //std::cout << "I_u : \n" << bounded(0, cols - 1, I_u).reshaped(grid.getRows(), grid.getCols()) << std::endl;
        Eigen::ArrayXi J_u = (uys / dy).floor().cast<int>();
        // lower left
        Eigen::ArrayXd u00 = u_f(grid.indexUVec(bounded(0, cols, I_u), bounded(0, rows - 1, J_u)));
        //std::cout << "u00 : \n" << u00.reshaped(grid.getRows(), grid.getCols()) << std::endl;

        // upper left
        Eigen::ArrayXd u01 = u_f(grid.indexUVec(bounded(0, cols, I_u), bounded(0, rows - 1, J_u + 1)));
        //std::cout << "u01 : \n" << u01.reshaped(grid.getRows(), grid.getCols()) << std::endl;

        // lower right
        Eigen::ArrayXd u10 = u_f(grid.indexUVec(bounded(0, cols, I_u + 1), bounded(0, rows - 1, J_u)));
        //std::cout << "u10 : \n" << u10.reshaped(grid.getRows(), grid.getCols()) << std::endl;

        // upper right
        Eigen::ArrayXd u11 = u_f(grid.indexUVec(bounded(0, cols, I_u + 1), bounded(0, rows - 1, J_u + 1)));
        //std::cout << "u11 : \n" << u11.reshaped(grid.getRows(), grid.getCols()) << std::endl << "return value : " << std::endl;

        //Bilinear interpolation factors
        Eigen::ArrayXd A_u = (uxs / dx) - I_u.cast<double>();
        Eigen::ArrayXd B_u = (uys / dy) - J_u.cast<double>();
        // I think it is correct since Eigen Arrays operation are coefficient-wise
        // Maybe check the math, shouldn't be too hard to modify (check bilerp in OG covector)
        return bilinint(u00, u01, u10, u11, A_u, B_u); //(1 - A_u) * (1 - B_u) * u00 + (1 - A_u) * B_u * u01 + A_u * (1 - B_u) * u10 + A_u * B_u * u11;
    }

    Eigen::ArrayXd GetVelocityV(StGridv2::StGridv2 grid, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys){
        double dx = grid.getDx();
        double dy = grid.getDy();
        int rows = grid.getRows();
        int cols = grid.getCols();

        // for v (vertical component)
        Eigen::ArrayXd vxs = xs - 0.5 * dx;
        Eigen::ArrayXd vys = ys;
        Eigen::ArrayXi I_v = (vxs / dx).floor().max(0).min(cols - 2).cast<int>();
        Eigen::ArrayXi J_v = (vys / dy).floor().max(0).min(rows - 1).cast<int>();

        Eigen::ArrayXd v00 = v_f(grid.indexVVec(bounded(0, cols - 1, I_v), bounded(0, rows, J_v)));
        Eigen::ArrayXd v01 = v_f(grid.indexVVec(bounded(0, cols - 1, I_v), bounded(0, rows, J_v + 1)));
        Eigen::ArrayXd v10 = v_f(grid.indexVVec(bounded(0, cols - 1, I_v + 1), bounded(0, rows, J_v)));
        Eigen::ArrayXd v11 = v_f(grid.indexVVec(bounded(0, cols - 1, I_v + 1), bounded(0, rows, J_v + 1)));

        Eigen::ArrayXd A_v = (vxs / dx) - I_v.cast<double>();
        Eigen::ArrayXd B_v = (vys / dy) - J_v.cast<double>();

        return bilinint(v00,  v01, v10, v11, A_v, B_v);//(1 - A_v) * (1 - B_v) * v00 + (1 - A_v) * B_v * v01 + A_v * (1 - B_v) * v10 + A_v * B_v * v11;
    }

    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    GetVelocity(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, const Eigen::ArrayXd& xs, const Eigen::ArrayXd& ys) {
        return {GetVelocityU(grid, u_f, xs, ys), GetVelocityV(grid, v_f, xs, ys)};
    }

    // Using RK4
    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    FlowMapPsi(StGridv2::StGridv2 grid, Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, Eigen::ArrayXd xs, Eigen::ArrayXd ys,
               double dt) {
        double c1 = 1.0 / 6.0 * dt;
        double c2 = 1.0 / 3.0 * dt;
        double c3 = c2;
        double c4 = c1;

        auto [k1u, k1v] = GetVelocity(grid, u_f, v_f, xs, ys);
        //std::cout << "k1u : \n" << k1u.reshaped(xs.rows(), xs.cols()) << std::endl;
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

    void BasicFlow(StGridv2::StGridv2& grid, double dt) {
        grid.ApplyVelocityBoundaries();
        grid.ApplyPressureBoundaries();
        Eigen::ArrayXd u_f = grid.getU();
        Eigen::ArrayXd v_f = grid.getV();
        Eigen::ArrayXd xs = grid.getXs();
        Eigen::ArrayXd ys = grid.getYs();
        int rows = grid.getRows();
        int cols = grid.getCols();
        Eigen::ArrayXd nu = grid.getU();
        Eigen::ArrayXd nv = grid.getV();
        // Compute the Flow Map at the center of every cell
        auto [psix, psiy] = FlowMapPsi(grid, u_f, v_f, xs, ys, dt);
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
        nu.reshaped(rows, cols + 1)(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped() = GetVelocityU(grid, u_f, psix_Fu, psiy_Fu);
        nv.reshaped(rows + 1, cols)(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped() = GetVelocityV(grid, v_f, psix_Fv, psiy_Fv);

        // u(x) <- u(psi(x))
        grid.setU(nu);
        grid.setV(nv);
        grid.SolvePressure();
        grid.SolveMomentumEquation();
    }

    // First order Covector Fluids (Algorithm 1)
    void CovectorFluids1(StGridv2::StGridv2 &grid, double dt) {
        grid.ApplyVelocityBoundaries();
        grid.ApplyPressureBoundaries();
        // u_f <- u | Freeze Flow Velocity
        Eigen::ArrayXd u_f = grid.getU();
        Eigen::ArrayXd v_f = grid.getV();
        // u <- A_covec(u_f, u, dt) | Covector Lie Advection
        auto [nu, nv] = CovectorAdvectionSL(grid, grid.getU(), grid.getV(), u_f, v_f, dt);
        grid.setU(nu);
        grid.setV(nv);
        // u <- P(u) | Pressure Projection
        grid.SolvePressure();
        grid.SolveMomentumEquation();
    }

    std::tuple<Eigen::ArrayXd, Eigen::ArrayXd>
    CovectorAdvectionSL(StGridv2::StGridv2 grid, Eigen::ArrayXd u, Eigen::ArrayXd v,
                        Eigen::ArrayXd u_f, Eigen::ArrayXd v_f, double dt) {
        Eigen::ArrayXd xs = grid.getXs();
        Eigen::ArrayXd ys = grid.getYs();
        int rows = grid.getRows();
        int cols = grid.getCols();
        Eigen::ArrayXd nu = u;
        Eigen::ArrayXd nv = v;
        // Compute the Flow Map at the center of every cell
        auto [psix, psiy] = FlowMapPsi(grid, u_f, v_f, xs, ys, dt);
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
        nu.reshaped(rows, cols + 1)(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped() = (grid.getDiv_x() * psix.matrix()).array() * GetVelocityU(grid, u_f, psix_Fu, psiy_Fu);
        nv.reshaped(rows + 1, cols)(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped() = (grid.getDiv_y() * psiy.matrix()).array() * GetVelocityV(grid, v_f, psix_Fv, psiy_Fv);

        return {nu, nv};
    }
}