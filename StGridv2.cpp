//
// Created by boris on 11.01.24.
//

#include "StGridv2.h"

namespace StGridv2 {
    StGridv2::StGridv2(int rows, int cols, double length, double breadth, Boundaries u_bound, Boundaries v_bound,
                       Boundaries p_bound, double CFL) {
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

        std::unordered_set<int> obstacles_c = {};
        std::unordered_set<int> obstacles_u = {};
        std::unordered_set<int> obstacles_v = {};

        // Arrays of the coordinates of the cells' center (made in column major)
        xs = Eigen::ArrayXd::LinSpaced(rows, 0, rows - 1).replicate(1, cols).transpose().reshaped() * dx + 0.5 * dx;
        ys = Eigen::ArrayXd::LinSpaced(cols, 0, cols - 1).replicate(1, rows).reshaped() * dy + 0.5 * dy;

        u = Eigen::ArrayXd::Zero((rows + 2) * (cols + 1));
        v = Eigen::ArrayXd::Zero((rows + 1) * (cols + 2));

        p = Eigen::ArrayXd::Zero(rows * cols);

        //ApplyVelocityBoundaries();
        //ApplyPressureBoundaries();
    }

    void StGridv2::GenerateOperators() {
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
        boundary_pressure.resize(p.size());
        // Warning : Have to be changed if obstacles_c(create a check validity for u and v) and updated dynamically if dynamic obstacles_c
        typedef Eigen::Triplet<double> T;
        std::vector<T> trip_lap;
        std::vector<T> trip_div_u;
        std::vector<T> trip_div_v;
        std::vector<T> trip_div_cx;
        std::vector<T> trip_div_cy;
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                int index = j + i * rows; //i + j * cols; <- Row Major

                // for div op along u -> (div u)_ij = (u_i+1j - u_ij) / dx
                int index_u = indexCU(i, j);
                if (!obstacles_u.count(indexU(i, j))){
                    trip_div_u.push_back(T(index, index_u, -1.0 / dx));
                }
                int indNextU = indexCU(i + 1, j);
                if (!obstacles_u.count(indexU(i + 1, j))){
                    trip_div_u.push_back(T(index, indNextU, 1.0 / dx));
                }

                // for div op along v -> (div v)_ij = (v_ij+1 - v_ij) / dy
                int index_v = indexCV(i, j);
                if (!obstacles_v.count(indexV(i, j))){
                    trip_div_v.push_back(T(index, index_v, -1.0 / dy));
                }
                int indNextV = indexCV(i, j + 1);
                if (!obstacles_v.count(indexV(i, j + 1))){
                    trip_div_v.push_back(T(index, indNextV, 1.0 / dy));
                }

                if (i < cols - 1) {
                    // for div op for values at the center of cells (p and psi) -> (div_x p)_ij = (p_i+1j - p_ij) / dx
                    int ind = j + i * rows;
                    int index_p = j + i * rows;
                    if (!obstacles_c.count(index_p)) {
                        trip_div_cx.push_back(T(ind, index_p, -1.0 / dx));
                    }
                    int indexNextP = j + (i + 1) * rows;
                    if(!obstacles_c.count(indexNextP)){
                        trip_div_cx.push_back(T(ind, indexNextP, 1.0 / dx));
                    }
                }
                if (j < rows - 1) {
                    // (div_y p)_ij = (p_ij+1 - p_ij) / dy
                    int ind = j + i * (rows - 1);
                    int index_p = j + i * rows;
                    if (!obstacles_c.count(index_p)){
                        trip_div_cy.push_back(T(ind, index_p, -1.0 / dy));
                    }
                    int indexNextP = j + 1 + i * rows;
                    if (!obstacles_c.count(indexNextP)){
                        trip_div_cy.push_back(T(ind, indexNextP, 1.0 / dy));
                    }
                }

                // for Laplacian op -> (Lap p)_ij = (p_i-1j - 2p_ij + p_i+1j) / dx² + (p_ij-1 - 2p_ij + p_ij+1) / dy²
                int index_p = indexP(i, j);
                // if ij is not in an obstacle
                if (!obstacles_c.count(index_p)) {
                    trip_lap.push_back(T(index_p, index_p, -2.0 / (dx * dx) - 2.0 / (dy * dy)));

                    // Check if the celle on the left is out of boundary or in an obstacle
                    if (i - 1 < 0) {
                        if (p_bound.left.type == "N") {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dx * dx)));
                        } else if (p_bound.left.type == "D") {
                            boundary_pressure.coeffRef(index_p) = p_bound.left.value / (dx * dx);
                        }
                    } else {
                        int indexLeft = indexP(i - 1, j);
                        if (obstacles_c.count(indexLeft)) {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dx * dx)));
                        } else {
                            trip_lap.push_back(T(index_p, indexLeft, 1.0 / (dx * dx)));
                        }
                    }
                    // on the right
                    if (i + 1 >= cols) {
                        if (p_bound.right.type == "N") {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dx * dx)));
                        } else if (p_bound.right.type == "D") {
                            boundary_pressure.coeffRef(index_p) = p_bound.right.value / (dx * dx);
                        }
                    } else {
                        int indexRight = indexP(i + 1, j);
                        if (obstacles_c.count(indexRight)) {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dx * dx)));
                        } else {
                            trip_lap.push_back(T(index_p, indexRight, 1.0 / (dx * dx)));
                        }
                    }

                    // on the bottom
                    if (j - 1 < 0) {
                        if (p_bound.down.type == "N") {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dy * dy)));
                        } else if (p_bound.down.type == "D") {
                            boundary_pressure.coeffRef(index_p) = p_bound.down.value / (dy * dy);
                        }
                    } else {
                        int indexBelow = indexP(i, j - 1);
                        if (obstacles_c.count(indexBelow)) {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dy * dy)));
                        } else {
                            trip_lap.push_back(T(index_p, indexBelow, 1.0 / (dy * dy)));
                        }
                    }

                    // on the top
                    if (j + 1 >= rows) {
                        if (p_bound.up.type == "N") {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dy * dy)));
                        } else if (p_bound.up.type == "D") {
                            boundary_pressure.coeffRef(index_p) = p_bound.up.value / (dy * dy);
                        }
                    } else {
                        int indexUp = indexP(i, j + 1);
                        if (obstacles_c.count(indexUp)) {
                            trip_lap.push_back(T(index_p, index_p, 1.0 / (dy * dy)));
                        } else {
                            trip_lap.push_back(T(index_p, indexUp, 1.0 / (dy * dy)));
                        }
                    }
                } else {
                    trip_lap.push_back(T(index_p, index_p, 1.0 ));
                }
            }
        }
        laplacian.setFromTriplets(trip_lap.begin(), trip_lap.end());
        div_u.setFromTriplets(trip_div_u.begin(), trip_div_u.end());
        div_v.setFromTriplets(trip_div_v.begin(), trip_div_v.end());
        div_cx.setFromTriplets(trip_div_cx.begin(), trip_div_cx.end());
        div_cy.setFromTriplets(trip_div_cy.begin(), trip_div_cy.end());
    }

    void StGridv2::AddObstacle(std::unordered_set<std::tuple<int, int>> indices) {
        for (std::tuple<int, int> index : indices){
            int i = std::get<0>(index), j = std::get<1>(index);

            this->obstacles_c.insert(indexP(i, j));
            // we add the indices of every faces of the cells in an obstacle
            this->obstacles_u.insert(indexU(i, j));
            this->obstacles_u.insert(indexU(i + 1, j));
            this->obstacles_v.insert(indexV(i, j));
            this->obstacles_v.insert(indexV(i, j + 1));
        }
    }

    void StGridv2::AddCircleObstacle(double cx, double cy, double r) {
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                if (std::pow(i * dx - cx, 2) + std::pow(j * dy - cy, 2) <= std::pow(r, 2)) {
                    this->obstacles_c.insert(indexP(i, j));
                    // we add the indices of every faces of the cells in an obstacle
                    this->obstacles_u.insert(indexU(i, j));
                    this->obstacles_u.insert(indexU(i + 1, j));
                    this->obstacles_v.insert(indexV(i, j));
                    this->obstacles_v.insert(indexV(i, j + 1));
                }
            }
        }
    }

    void StGridv2::AddSquareObstacle(double cx, double cy, double l) {
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                if ((i * dx - cx > -l/2) and (i * dx - cx < l/2) and (j * dy - cy > -l/2) and (j * dy - cy < l/2)) {
                    this->obstacles_c.insert(indexP(i, j));
                    // we add the indices of every faces of the cells in an obstacle
                    this->obstacles_u.insert(indexU(i, j));
                    this->obstacles_u.insert(indexU(i + 1, j));
                    this->obstacles_v.insert(indexV(i, j));
                    this->obstacles_v.insert(indexV(i, j + 1));
                }
            }
        }
        std::cout << obstacles_c.size() << std::endl;
    }

    // Obsolete with RG45
    void StGridv2::SetTimeStep() {
        double factor = u.abs().maxCoeff() / dx + v.abs().maxCoeff() / dy;
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
        if (u_up.type == "D") {
            // Done u_imax+1 = 2 * value - u_imax
            // Since the bottom and top line of u are not on the boundary, we do a linear interpolation between the boundary value and the line above or below respectively
            // u_imax = (2 * top value + u_imax-1)/ 3

            // Since the bottom and top line of u are not on the boundary, we had a virtual line at the top such that the value at the boundary would be as desired
            // 0.5 * (u_imax+1 + u_imax) = value -> u_imax+1 = 2 * value - u_imax
            u.reshaped(rows + 2, cols + 1)(Eigen::last, Eigen::all) =
                    2.0 * u_up.value - u.reshaped(rows + 2, cols + 1)(Eigen::last - 1, Eigen::all);
        } else if (u_up.type == "N") {
            u.reshaped(rows + 2, cols + 1)(Eigen::last, Eigen::all) =
                    -u_up.value * dy + u.reshaped(rows + 2, cols + 1)(Eigen::last - 1, Eigen::all);
        } else if (u_up.type == "C") {
            u.reshaped(rows + 2, cols + 1)(Eigen::last, Eigen::all) = u.reshaped(rows + 2, cols + 1)(0, Eigen::all);
        }
        if (u_down.type == "D") {
            u.reshaped(rows + 2, cols + 1)(0, Eigen::all) =
                    2.0 * u_down.value - u.reshaped(rows + 2, cols + 1)(1, Eigen::all);
        } else if (u_down.type == "N") {
            u.reshaped(rows + 2, cols + 1)(0, Eigen::all) =
                    -u_down.value * dy + u.reshaped(rows + 2, cols + 1)(1, Eigen::all);
        } else if (u_down.type == "C") {
            u.reshaped(rows + 2, cols + 1)(0, Eigen::all) = u.reshaped(rows + 2, cols + 1)(Eigen::last, Eigen::all);
        }
        if (u_left.type == "D") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), 0) = u_left.value;
        } else if (u_left.type == "N") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), 0) =
                    -u_left.value * dx + u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), 1);
        } else if (u_left.type == "C") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), 0) = u.reshaped(rows + 2, cols + 1)(
                    Eigen::seq(1, Eigen::last - 1), Eigen::last);
        }
        if (u_right.type == "D") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::last) = u_right.value;
        } else if (u_right.type == "N") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::last) = -u_right.value * dx + u.reshaped(rows + 2,cols + 1)(Eigen::seq(1,Eigen::last - 1),Eigen::last - 1);
        } else if (u_right.type == "C") {
            u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::last) = u.reshaped(rows + 2,cols + 1)(Eigen::seq(1, Eigen::last - 1), 0);
        }

        // For v
        Boundary v_left = v_bound.left;
        Boundary v_right = v_bound.right;
        Boundary v_up = v_bound.up;
        Boundary v_down = v_bound.down;

        if (v_left.type == "D") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, 0) =
                    2.0 * v_left.value - v.reshaped(rows + 1, cols + 2)(Eigen::all, 1);
        } else if (v_left.type == "N") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, 0) =
                    -v_left.value * dx + v.reshaped(rows + 1, cols + 2)(Eigen::all, 1);
        } else if (v_left.type == "C") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, 0) = v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last);
        }
        if (v_right.type == "D") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last) =
                    2.0 * v_right.value - v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last - 1);
        } else if (v_right.type == "N") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last) =
                    -v_right.value * dx + v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last - 1);
        } else if (v_right.type == "C") {
            v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::last) = v.reshaped(rows + 1, cols + 2)(Eigen::all, 0);
        }
        if (v_up.type == "D") {
            v.reshaped(rows + 1, cols + 2)(Eigen::last, Eigen::seq(1, Eigen::last - 1)) = v_up.value;
        } else if (v_up.type == "N") {
            v.reshaped(rows + 1, cols + 2)(Eigen::last, Eigen::seq(1, Eigen::last - 1)) =
                    -v_up.value * dy + v.reshaped(rows + 1, cols + 2)(Eigen::last - 1, Eigen::seq(1, Eigen::last - 1));
        } else if (v_up.type == "C") {
            v.reshaped(rows + 1, cols + 2)(Eigen::last, Eigen::seq(1, Eigen::last - 1)) = v.reshaped(rows + 1,cols + 2)(0,Eigen::seq(1,Eigen::last - 1));
        }
        if (v_down.type == "D") {
            v.reshaped(rows + 1, cols + 2)(0, Eigen::seq(1, Eigen::last - 1)) = v_down.value;
        } else if (v_down.type == "N") {
            v.reshaped(rows + 1, cols + 2)(0, Eigen::seq(1, Eigen::last - 1)) =
                    -v_down.value * dy + v.reshaped(rows + 1, cols + 2)(1, Eigen::seq(1, Eigen::last - 1));
        } else if (v_down.type == "C") {
            v.reshaped(rows + 1, cols + 2)(0, Eigen::seq(1, Eigen::last - 1)) = v.reshaped(rows + 1, cols + 2)(Eigen::last, Eigen::seq(1, Eigen::last - 1));
        }
    }

    void StGridv2::ApplyVelocityObstacles() {
        // set to 0 the velocity inside obstacles
        for(const auto index : obstacles_u){
            u(index) = 0;
        }
        for(const auto index : obstacles_v){
            v(index) = 0;
        }
    }

    // compute the divergence of the velocity D V
    Eigen::VectorXd StGridv2::ComputeDivergence() {
        return div_u * u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped().matrix()
               + div_v * v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped().matrix();
    }

    // Solve L p = D V for p with Cholesky solver
    // L is the Laplacian operator, p is the pressure, D is the divergence operator and V is the velocity
    void StGridv2::SolvePressureLin() {
        // Solve A x = b where L -> A, p -> x, D V -> b
        // dividing by dt reduce div (?)
        Eigen::VectorXd b = 1.0 / dt * ComputeDivergence() - boundary_pressure;
        Eigen::VectorXd x(rows * cols);
        // Correction to make the system solvable
        // L <- L + a I where a = 10E-6 * trace(L) / n (size of L = n x n)
        Eigen::SparseMatrix<double> A = laplacian;
        A.diagonal().array() += 10E-6 * A.diagonal().sum() / A.rows();
        // Cholesky solver
        Eigen::SimplicialCholesky<Eigen::SparseMatrix<double>> solver;
        solver.analyzePattern(A);
        solver.factorize(A);
        x = solver.solve(b);
        p = x.array();
    }

    void StGridv2::SolvePressureIter() {
        // Solve A x = b where L -> A, p -> x, D V -> b
        // dividing by dt reduce div (?)
        Eigen::VectorXd b = 1.0 / dt * ComputeDivergence() - boundary_pressure;
        Eigen::VectorXd x(rows * cols);
        // Correction to make the system solvable
        // L <- L + a I where a = 10E-6 * trace(L) / n (size of L = n x n)
        Eigen::SparseMatrix<double> A = laplacian;
        A.diagonal().array() += 10E-6 * A.diagonal().sum() / A.rows();
        // Iterative solver using Conjugate Gradient from Eigen lib
        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper> iter;
        //iter.setMaxIterations(100); // TODO check good value ?
        iter.setTolerance(0.1);
        iter.compute(A);
        x = iter.solve(b);
        //std::cout << "#iterations:     " << iter.iterations() << std::endl;
        //std::cout << "estimated error: " << iter.error()      << std::endl;
        p = x.array();
    }

    // Pressure projection
    void StGridv2::SolveMomentumEquation() {
        //DONE resize u
        // u <- u - dt * div_x p
        u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() +=
                -dt * (div_cx * p.matrix()).array().reshaped();

        // v <- v - dt * div_y p
        v.reshaped(rows + 1, cols + 2)(Eigen::seq(1, Eigen::last - 1), Eigen::seq(1, Eigen::last - 1)).reshaped() +=
                -dt * (div_cy * p.matrix()).array().reshaped();
    }


    // ======= Getter and Setter =======
    int StGridv2::getRows() const {
        return rows;
    }

    int StGridv2::getCols() const {
        return cols;
    }

    double StGridv2::getBreadth() const {
        return breadth;
    }

    double StGridv2::getLength() const {
        return length;
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

    // Return the index of the velocity u on the left face of a cell at coord i,j
    int StGridv2::indexU(int i, int j) {
        return j + 1 + i * (rows + 2);
    }

    // Return the index of the velocity v on the bottom face of a cell at coord i,j
    int StGridv2::indexV(int i, int j) {
        return j + (i + 1) * (rows + 1);
    }

    // Done resize u
    // Gives the index k for an array (rows * (cols + 1)) corresponding to the index i, j of a matrix (rows X cols + 1)
    // This is used to compute the index i j on the centered u array (i.e. without the virtual rows)
    int StGridv2::indexCU(int i, int j) const {
        // u is rows x cols+1
        return j + i * rows; //i + j * (cols + 1); <- in row major
    }

    // Gives the index k for an array ((rows + 1) * cols) corresponding to the index i, j of a matrix (rows + 1 X cols)
    // This is used to compute the index i j on the centered v array (i.e. without the virtual cols)
    // The difference is less important for v than for u since it is column major
    int StGridv2::indexCV(int i, int j) {
        return j + i * (rows + 1); //i + j * cols; <- in row major
    }

    // Gives the index k for an array ((rows + 2) * (cols + 1)) corresponding to the index i, j of a matrix (rows + 2 X cols + 1)
    // This is used to compute the index i j on the full u array (i.e. with the virtual rows)
    Eigen::ArrayXi StGridv2::indexFUVec(Eigen::ArrayXi I, Eigen::ArrayXi J) {
        return J + I * (rows + 2); //I + J * (cols + 1); <- in row major
    }

    // Gives the index k for an array ((rows + 1) * (cols + 2)) corresponding to the index i, j of a matrix (rows + 1 X cols + 2)
    // This is used to compute the index i j on the full v array (i.e. with the virtual cols)
    // The difference is less important for v than for u since it is column major
    Eigen::ArrayXi StGridv2::indexFVVec(Eigen::ArrayXi I, Eigen::ArrayXi J) {
        return J + I * (rows + 1); //I + J * cols; <- in row major
    }

    // VTK file writer for Paraview vizualisation
    void StGridv2::WriteToFile(std::string name) {
        int N = cols * rows;
        std::ofstream myfile;
        myfile.open(name);
        myfile << "# vtk DataFile Version 2.0" << std::endl;
        myfile << "First test" << std::endl;
        myfile << "ASCII" << std::endl;
        myfile << "DATASET STRUCTURED_GRID" << std::endl;
        myfile << "DIMENSIONS " << cols << " " << rows << " 1" << std::endl;

        myfile << "POINTS " << N << " double" << std::endl;
        for (int i = 0; i < N; ++i) {
            myfile << xs(i) << " " << ys(i) << " " << 0 << std::endl;
        }
        myfile << "POINT_DATA " << N << std::endl;
        myfile << "SCALARS pressure double 1" << std::endl;
        myfile << "LOOKUP_TABLE default" << std::endl;
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                myfile << p(indexP(i, j)) << std::endl;
            }
        }
        myfile << "VECTORS velocity double" << std::endl;
        // u_c is the centered u array (i.e. without the virtual rows)
        Eigen::ArrayXd u_c = u.reshaped(rows + 2, cols + 1)(Eigen::seq(1, Eigen::last - 1), Eigen::all).reshaped();
        Eigen::ArrayXd v_c = v.reshaped(rows + 1, cols + 2)(Eigen::all, Eigen::seq(1, Eigen::last - 1)).reshaped();
        for (int i = 0; i < cols; ++i) {
            for (int j = 0; j < rows; ++j) {
                //DONE resize u
                //std::cout << "u size : " << u_c.rows() << " " << u_c.cols() << std::endl;
                double uij = (u_c(indexCU(i, j)) + u_c(indexCU(i + 1, j))) * 0.5;
                double vij = (v_c(indexCV(i, j)) + v_c(indexCV(i, j + 1))) * 0.5;
                myfile << uij << " " << vij << " " << 0.0 << std::endl;
            }
        }
        myfile.close();
    }
}