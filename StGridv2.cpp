//
// Created by boris on 11.01.24.
//

#include "StGridv2.h"

namespace StGridv2 {
    StGridv2::StGridv2(int rows, int cols, double length, double breadth, Boundaries u_bound, Boundaries v_bound, Boundaries p_bound, double CFL) {
        this->rows = rows;
        this->cols = cols;
        this->length = length;
        this->breadth = breadth;
        this->u_bound = u_bound;
        this->v_bound = v_bound;
        this->p_bound = p_bound;
        this->CFL = CFL;

        dx = length / cols;
        dy = breadth / rows;
        // TODO check if fine ?
        dt = CFL * (dx + dy);

        // Arrays of the coordinates of the cells' center (made in column major)
        xs = Eigen::ArrayXd::LinSpaced(rows, 0, rows - 1).replicate(1, cols).transpose().reshaped() * dx + 0.5 * dx;
        ys = Eigen::ArrayXd::LinSpaced(cols, 0, cols - 1).replicate(1, rows).reshaped() * dy + 0.5 * dy;

        u = Eigen::ArrayXd::Zero(rows * (cols + 1));
        v = Eigen::ArrayXd::Zero((rows + 1) * cols);

        p = Eigen::ArrayXd::Zero(rows * cols);

        ApplyVelocityBoundaries();
        ApplyPressureBoundaries();

        // Compute laplacian and divergence operator matrix for L p = D V
        // Size analysis :
        // L -> rows * cols x rows * cols
        // p -> rows * cols x 1
        // div_u -> rows * cols x rows * (cols + 1)
        // div_v -> rows * cols x (rows + 1) * cols
        // D -> rows * cols x rows * (cols + 1) + (rows + 1) * cols (append colwise div_u and div_v)
        // V -> rows * (cols + 1) + (rows + 1) * cols x 1 (V = u and v rowwise appended)
        // L p -> rows * cols x 1
        // D v -> rows * cols x 1
        laplacian.resize(cols * rows, cols * rows);
        div_u.resize(cols * rows, (cols + 1) * rows);
        div_v.resize(cols * rows, cols * (rows + 1));
        div_cx.resize((cols - 1) * rows, cols * rows);
        div_cy.resize(cols * (rows - 1), cols * rows);
        // Warning : Have to be changed if obstacles(create a check validity for u and v) and updated dynamically if dynamic obstacles
        typedef Eigen::Triplet<double> T;
        std::vector<T> trip_lap;
        std::vector<T> trip_div_u;
        std::vector<T> trip_div_v;
        std::vector<T> trip_div_cx;
        std::vector<T> trip_div_cy;
        for(int i = 0; i < cols; ++i){
            for(int j = 0; j < rows; ++j){
                //std::cout << "i : " << i << " j : " << j << std::endl;

                int index = j + i * rows; //i + j * cols; <- Row Major

                // for div op along u -> (div u)_ij = (u_i+1j - u_ij) / dx
                int index_u = indexU(i ,j);
                //div_u.insert(index, index_u) = -1.0 / dx;
                trip_div_u.push_back(T(index, index_u, -1.0 / dx));
                int indNextU = indexU(i + 1, j);
                //div_u.insert(index, indNextU) = 1.0 / dx;
                trip_div_u.push_back(T(index, indNextU, 1.0 / dx));
                // for div op along v -> (div v)_ij = (v_ij+1 - v_ij) / dy
                int index_v = indexV(i, j);
                //div_v.insert(index, index_v) = -1.0 / dy;
                trip_div_v.push_back(T(index, index_v, -1.0 / dy));
                int indNextV = indexV(i, j + 1);
                //div_v.insert(index, indNextV) = 1.0 / dy;
                trip_div_v.push_back(T(index, indNextV, 1.0 / dy));

                if (i < cols - 1) {
                    // for div op for values at the center of cells (p and psi) -> (div_x p)_ij = (p_i+1j - p_ij) / dx
                    int ind = j + i * rows;
                    int index_p = j + i * rows; //indexP(i, j);
                    //div_cx.insert(index, index_p) = -1.0 / dx;
                    trip_div_cx.push_back(T(ind, index_p, -1.0 / dx));
                    int indexNextP = j + (i + 1) * rows; //indexP(i + 1, j);
                    //div_cx.insert(index, indexNextP) = 1.0 / dx;
                    trip_div_cx.push_back(T(ind, indexNextP, 1.0 / dx));
                    //std::cout << "In div_cx -> index : " << index << " | index_p : "<< index_p << " | indexNextP : " << indexNextP << std::endl;
                }
                if (j < rows - 1){
                    // (div_y p)_ij = (p_ij+1 - p_ij) / dy
                    int ind = j + i * (rows - 1);
                    int index_p = j + i * rows; //indexP(i, j);
                    //div_cy.insert(index, index_p) = -1.0 / dy;
                    trip_div_cy.push_back(T(ind, index_p, -1.0 / dy));
                    int indexNextP = j + 1 + i * rows; //indexP(i, j + 1);
                    //div_cy.insert(index, indexNextP) = 1.0 / dy;
                    trip_div_cy.push_back(T(ind, indexNextP, 1.0 / dy));
                    //std::cout << "In div_cy -> index : " << index << " | index_p : "<< index_p << " | indexNextP : " << indexNextP << std::endl;
                }

                // for Laplacian op -> (Lap p)_ij = (p_i-1j + p_i+1j + p_ij-1 + p_ij+1 - 4 p_ij) / (dx * dy)
                double factor = dx * dy;
                int index_p = indexP(i, j);
                // Check validity for i,j ? Check if in obstacle ?
                //laplacian.insert(index_p, index_p) = -4.0;
                trip_lap.push_back(T(index_p, index_p, -4.0 / factor));
                if (checkPvalid(i - 1, j)){
                    int indexLeft = indexP(i-1, j);
                    //laplacian.insert(index_p, indexLeft) = 1.0;
                    trip_lap.push_back(T(index_p, indexLeft, 1.0 / factor));
                } else {
                    //laplacian.coeffRef(index_p, index_p) += 1.0;
                    trip_lap.push_back(T(index_p, index_p, 1.0 / factor));
                }
                if (checkPvalid(i + 1, j)){
                    int indexRight = indexP(i+1, j);
                    //laplacian.insert(index_p, indexRight) = 1.0;
                    trip_lap.push_back(T(index_p, indexRight, 1.0 / factor));
                } else {
                    //laplacian.coeffRef(index_p, index_p) += 1.0;
                    trip_lap.push_back(T(index_p, index_p, 1.0 / factor));
                }
                if (checkPvalid(i, j - 1)){
                    int indexBelow = indexP(i, j-1);
                    //laplacian.insert(index_p, indexBelow) = 1.0;
                    trip_lap.push_back(T(index_p, indexBelow, 1.0 / factor));
                } else {
                    //laplacian.coeffRef(index_p, index_p) += 1.0;
                    trip_lap.push_back(T(index_p, index_p, 1.0 / factor));
                }
                if (checkPvalid(i, j + 1)){
                    int indexUp = indexP(i, j+1);
                    //laplacian.insert(index_p, indexUp) = 1.0;
                    trip_lap.push_back(T(index_p, indexUp, 1.0 / factor));
                } else {
                    //laplacian.coeffRef(index_p, index_p) += 1.0;
                    trip_lap.push_back(T(index_p, index_p, 1.0 / factor));
                }
            }
        }
        laplacian.setFromTriplets(trip_lap.begin(), trip_lap.end());
        div_u.setFromTriplets(trip_div_u.begin(), trip_div_u.end());
        div_v.setFromTriplets(trip_div_v.begin(), trip_div_v.end());
        div_cx.setFromTriplets(trip_div_cx.begin(), trip_div_cx.end());
        div_cy.setFromTriplets(trip_div_cy.begin(), trip_div_cy.end());
        //std::cout << "Div x : " << div_u << std::endl;
    }

    void StGridv2::SetTimeStep() {
        double factor = u.maxCoeff() / dx + v.maxCoeff() / dy;
        if (factor == 0) {
            dt = CFL * (dx + dy);
        } else {
            dt = CFL / factor;
        }
    }

    // Apply the Boundary conditions for u and v
    void StGridv2::ApplyVelocityBoundaries() {
        // For u
        Boundary u_left = u_bound.left;
        Boundary u_right = u_bound.right;
        Boundary u_up = u_bound.up;
        Boundary u_down = u_bound.down;

        // D for Dirichlet
        // N for Von Neumann | (u_imax-1 - u_imax) / dx = val -> u_imax = -val*dx + u_imax-1
        // C for Cyclic | U_imax = U_0
        if (u_up.type == "D"){
            // Since the bottom and top line of u are not on the boundary, we do a linear interpolation between the boundary value and the line above or below respectively
            // u_imax = (2 * top value + u_imax-1)/ 3
            u.reshaped(rows, cols + 1)(Eigen::last, Eigen::all) = (2.0 * u_up.value + u.reshaped(rows, cols + 1)(Eigen::last - 1, Eigen::all)) / 3.0;
            //std::cout << "u : \n" << u.reshaped(rows, cols + 1) << std::endl;
        } else if (u_up.type == "N"){
            u.reshaped(rows, cols + 1)(Eigen::last, Eigen::all) = - u_up.value * dy + u.reshaped(rows, cols + 1)(Eigen::last - 1, Eigen::all);
        } else if (u_up.type == "C"){
            u.reshaped(rows, cols + 1)(Eigen::last, Eigen::all) = u.reshaped(rows, cols + 1)(0, Eigen::all);
        }
        if (u_down.type == "D"){
            u.reshaped(rows, cols + 1)(0, Eigen::all) = (2.0 * u_down.value + u.reshaped(rows, cols + 1)(1, Eigen::all)) / 3.0;
        } else if (u_down.type == "N"){
            u.reshaped(rows, cols + 1)(0, Eigen::all) = - u_down.value * dy + u.reshaped(rows, cols + 1)(1, Eigen::all);
        } else if (u_down.type == "C"){
            u.reshaped(rows, cols + 1)(0, Eigen::all) = u.reshaped(rows, cols + 1)(Eigen::last, Eigen::all);
        }
        if (u_left.type == "D"){
            u.reshaped(rows, cols + 1)(Eigen::all, 0) = u_left.value;
        } else if (u_left.type == "N"){
            u.reshaped(rows, cols + 1)(Eigen::all, 0) = - u_left.value * dx + u.reshaped(rows, cols + 1)(Eigen::all, 1);
        } else if (u_left.type == "C"){
            u.reshaped(rows, cols + 1)(Eigen::all, 0) = u.reshaped(rows, cols + 1)(Eigen::all, Eigen::last);
        }
        if (u_right.type == "D") {
            u.reshaped(rows, cols + 1)(Eigen::all, Eigen::last) = u_right.value;
        } else if (u_right.type == "N"){
            u.reshaped(rows, cols + 1)(Eigen::all, Eigen::last) = - u_right.value * dx + u.reshaped(rows, cols + 1)(Eigen::all, Eigen::last - 1);
        } else if (u_right.type =="C"){
            u.reshaped(rows, cols + 1)(Eigen::all, Eigen::last) = u.reshaped(rows, cols + 1)(Eigen::all, 0);
        }

        // For v
        Boundary v_left = v_bound.left;
        Boundary v_right = v_bound.right;
        Boundary v_up = v_bound.up;
        Boundary v_down = v_bound.down;

        if (v_left.type == "D"){
            v.reshaped(rows + 1, cols)(Eigen::all, 0) = (2.0 * v_left.value + v.reshaped(rows + 1, cols)(Eigen::all, 1)) / 3.0;
        } else if (v_left.type == "N"){
            v.reshaped(rows + 1, cols)(Eigen::all, 0) = - v_left.value * dx + v.reshaped(rows + 1, cols)(Eigen::all, 1);
        } else if (v_left.type == "C"){
            v.reshaped(rows + 1, cols)(Eigen::all, 0) = v.reshaped(rows + 1, cols)(Eigen::all, Eigen::last);
        }
        if (v_right.type == "D"){
            v.reshaped(rows + 1, cols)(Eigen::all, Eigen::last) = (2.0 * v_right.value + v.reshaped(rows + 1, cols)(Eigen::all, Eigen::last - 1)) / 3.0;
        } else if (v_right.type == "N"){
            v.reshaped(rows + 1, cols)(Eigen::all, Eigen::last) = - v_right.value * dx + v.reshaped(rows + 1, cols)(Eigen::all,Eigen::last - 1);
        } else if (v_right.type =="C"){
            v.reshaped(rows + 1, cols)(Eigen::all, Eigen::last) = v.reshaped(rows + 1, cols)(Eigen::all, 0);
        }
        if (v_up.type == "D"){
            v.reshaped(rows + 1, cols)(Eigen::last, Eigen::all) = v_up.value;
        } else if (v_up.type == "N"){
            v.reshaped(rows + 1, cols)(Eigen::last, Eigen::all) = - v_up.value * dy + v.reshaped(rows + 1, cols)(Eigen::last - 1, Eigen::all);
        } else if (v_up.type == "C"){
            v.reshaped(rows + 1, cols)(Eigen::last, Eigen::all) = v.reshaped(rows + 1, cols)(0, Eigen::all);
        }
        if (v_down.type == "D"){
            v.reshaped(rows + 1, cols)(0, Eigen::all) = v_down.value;
        } else if (v_down.type == "N"){
            v.reshaped(rows + 1, cols)(0, Eigen::all) = - v_down.value * dy + v.reshaped(rows + 1, cols)(1, Eigen::all);
        } else if (v_down.type == "C"){
            v.reshaped(rows + 1, cols)(0, Eigen::all) = v.reshaped(rows + 1, cols)(Eigen::last, Eigen::all);
        }
    }

    // apply the boundary condition for the pressure
    void StGridv2::ApplyPressureBoundaries() {
        Boundary p_left = p_bound.left;
        Boundary p_right = p_bound.right;
        Boundary p_up = p_bound.up;
        Boundary p_down = p_bound.down;

        if (p_left.type == "D"){
            p.reshaped(rows, cols)(Eigen::all, 0) = p_left.value;
        } else if (p_left.type == "N"){
            p.reshaped(rows, cols)(Eigen::all, 0) = - p_left.value * dx + p.reshaped(rows, cols)(Eigen::all ,1);
        } else if (p_left.type == "C"){
            p.reshaped(rows, cols)(Eigen::all, 0) = p.reshaped(rows, cols)(Eigen::all, Eigen::last);
        }
        if (p_right.type == "D"){
            p.reshaped(rows, cols)(Eigen::all, Eigen::last) = p_right.value;
        } else if (p_right.type == "N"){
            p.reshaped(rows, cols)(Eigen::all, Eigen::last) = -p_right.value * dx + p.reshaped(rows, cols)(Eigen::all, Eigen::last - 1);
        } else if (p_right.type == "C"){
            p.reshaped(rows, cols)(Eigen::all, Eigen::last) = p.reshaped(rows, cols)(Eigen::all, 0);
        }
        if (p_up.type == "D"){
            p.reshaped(rows, cols)(Eigen::last, Eigen::all) = p_up.value;
        } else if(p_up.type == "N"){
            p.reshaped(rows, cols)(Eigen::last, Eigen::all) = - p_up.value * dy + p.reshaped(rows, cols)(Eigen::last - 1, Eigen::all);
        }  else if (p_up.type == "C"){
            p.reshaped(rows, cols)(Eigen::last, Eigen::all) = p.reshaped(rows, cols)(0, Eigen::all);
        }
        if (p_down.type == "D"){
            p.reshaped(rows, cols)(0, Eigen::all) = p_down.value;
        } else if (p_down.type == "N"){
            p.reshaped(rows, cols)(0, Eigen::all) = - p_down.value * dy + p.reshaped(rows, cols)(1, Eigen::all);
        } else if (p_down.type == "C"){
            p.reshaped(rows, cols)(0, Eigen::all) = p.reshaped(rows, cols)(Eigen::last, Eigen::all);
        }
    }

    // compute the divergence of the velocity D V
    // TODO - Check if it works (check size of return value)
    Eigen::VectorXd StGridv2::ComputeDivergence() {
//        // append div_u and div_v
//        Eigen::MatrixXd div(rows * cols, rows * (cols + 1) + (rows + 1) * cols);
//        div << div_u, div_v;
//        // append u and v
//        Eigen::VectorXd V(rows * (cols + 1) + (rows + 1) * cols);
//        V << u, v;
//        return div * V;
        return div_u * u.matrix() + div_v * v.matrix();
    }

    // Solve L p = D V for p with LU factorization
    // L is the Laplacian operator, p is the pressure, D is the divergence operator and V is the velocity
    void StGridv2::SolvePressure() {
        // Solve A x = b where L -> A, p -> x, D V -> b
        Eigen::VectorXd b = ComputeDivergence();
        Eigen::VectorXd x(rows * cols);
        Eigen::SparseMatrix<double> A = laplacian;
        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        x = solver.solve(b);
        p = x.array();
        //std::cout << "p size : " << p.size() << std::endl;
        //ApplyPressureBoundaries();
    }

    // Pressure projection
    // Note : we pad div p because it only gives the divergence of the pressure on the interior faces
    void StGridv2::SolveMomentumEquation() {
        // u <- u - dt * div_x p
        //std::cout << "Max u value before "
        //Eigen::MatrixXd paded_divx_p = Eigen::MatrixXd::Zero(rows, cols + 1);
        //std::cout << "p : \n" << p.reshaped(rows, cols) << std::endl;
        //paded_divx_p(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped() = div_cx * p.matrix();
        //std::cout << "div_x p : \n" << paded_divx_p.reshaped(rows, cols + 1) << std::endl;
        //std::cout << "u before projection : " << u.reshaped(rows, cols + 1) << std::endl;
        //u = u - dt * paded_divx_p.reshaped().array();
        u.reshaped(rows, cols + 1)(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped() += -dt * (div_cx * p.matrix()).array().reshaped();
        //std::cout << "u after projection : " << u.reshaped(rows, cols + 1) << std::endl;

        // v <- v - dt * div_y p
        //Eigen::MatrixXd paded_divy_p = Eigen::MatrixXd::Zero(rows + 1, cols);
        //paded_divy_p(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped() = div_cy * p.matrix();
        //v = v - dt * paded_divy_p.reshaped().array();
        v.reshaped(rows + 1, cols)(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped() += -dt * (div_cy * p.matrix()).array().reshaped();
        // ApplyVelocityBoundaries(); ??
    }


    // ======= Getter and Setter =======
    int StGridv2::getRows() const {
        return rows;
    }

    int StGridv2::getCols() const {
        return cols;
    }

    double StGridv2::getDx() const {
        return dx;
    }

    double StGridv2::getDy() const {
        return dy;
    }

    double StGridv2::getDt() const {
        return dt;
    }

    Eigen::ArrayXd StGridv2::getXs() {
        return xs;
    }

    Eigen::SparseMatrix<double> StGridv2::getDiv_x() {
        return div_cx;
    }

    Eigen::SparseMatrix<double> StGridv2::getDiv_y() {
        return div_cy;
    }

    Eigen::ArrayXd StGridv2::getYs() {
        return ys;
    }

    const Eigen::ArrayXd &StGridv2::getU() const {
        return u;
    }

    const Eigen::ArrayXd &StGridv2::getV() const {
        return v;
    }

    const Eigen::ArrayXd &StGridv2::getP() const {
        return p;
    }

    void StGridv2::setDt(const double dt) {
        this->dt = dt;
    }

    void StGridv2::setU(const Eigen::ArrayXd &u) {
        this->u = u;
    }

    void StGridv2::setV(const Eigen::ArrayXd &v) {
        this->v = v;

    }

    // TODO : - Obstacles
    bool StGridv2::checkPvalid(int i, int j) const {
        if (i < 0 || i >= cols || j < 0 || j >= rows /* || in obstacle */) {
            return false;
        } else {
            return true;
        }
    }

    // The indices are computed in column major (as it is the default used by Eigen)
    int StGridv2::indexP(int i, int j) {
        return j + i * rows; //i + j * cols; <- in row major
    }

    int StGridv2::indexU(int i, int j) {
        // u is rows x cols+1
        return j + i * rows; //i + j * (cols + 1); <- in row major
    }

    int StGridv2::indexV(int i, int j) {
        return j + i * (rows + 1); //i + j * cols; <- in row major
    }

    Eigen::ArrayXi StGridv2::indexUVec(Eigen::ArrayXi I, Eigen::ArrayXi J) {
        return J + I * rows; //I + J * (cols + 1); <- in row major
    }

    Eigen::ArrayXi StGridv2::indexVVec(Eigen::ArrayXi I, Eigen::ArrayXi J) {
        return J + I * (rows + 1); //I + J * cols; <- in row major
    }

    void StGridv2::WriteToFile(std::string name){
        int N = cols * rows;
        std::ofstream myfile;
        myfile.open(name);
        myfile << "# vtk DataFile Version 2.0" << std::endl;
        myfile << "First test" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET STRUCTURED_GRID" << std::endl;
        myfile << "DIMENSIONS " << cols << " " << rows << " 1" << std::endl;

        myfile << "POINTS " << N << " double" << std::endl;
        for(int i = 0; i < N; ++i){
            myfile << xs(i) << " " << ys(i) << " " << 0 << std::endl;
            //std::cout << i << " | xs : " << xs(i) << " ys : " << ys(i) << std::endl;
        }
        myfile << "POINT_DATA " << N << std::endl;
        myfile << "SCALARS pressure double 1" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < cols; ++i){
            for (int j = 0; j < rows; ++j){
                myfile << p(indexP(i,j)) << std::endl;
            }
        }
        myfile << "VECTORS velocity double" << std::endl;
        //myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < cols; ++i){
            for (int j = 0; j < rows; ++j){
                // Is it ok to average at the center between two faces ?
                double uij = (u(indexU(i, j)) + u(indexU(i + 1, j))) * 0.5;
                double vij = (v(indexV(i, j)) + v(indexV(i, j + 1))) * 0.5;
                myfile << uij << " " << vij << " " << 0.0 << std::endl;
                //std::cout << indexP(i, j) << " | uij : " << uij << std::endl;
            }
        }
        myfile.close();
    }
}