#include "HSIM.h"
using namespace std;
using namespace Eigen;
using namespace Spectra;

pair<VectorXd, MatrixXd> SIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    MatrixXd &Phi,
    const int &p,
    const double &epsilon,
    const double &mu,
    const string &metric)
{
    SparseMatrix<double> shifted_S = S - mu * M;
    auto start = chrono::high_resolution_clock::now();
    CholmodSupernodalLLT<SparseMatrix<double>> solver;
    solver.compute(shifted_S);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time SIM LDLT shifted_S taken: " << duration.count() << " milliseconds" << endl;

    bool converged = false;
    VectorXd eigenvalues;

    while (!converged)
    {
        MatrixXd temp_Mphi = M * Phi;
        MatrixXd temp_phi = solver.solve(temp_Mphi);
        MatrixXd Reduced_S = temp_phi.transpose() * shifted_S * temp_phi;
        MatrixXd Reduced_M = temp_phi.transpose() * M * temp_phi;

        start = chrono::high_resolution_clock::now();
        GeneralizedSelfAdjointEigenSolver<MatrixXd> eigen_solver(Reduced_S, Reduced_M);
        eigenvalues = eigen_solver.eigenvalues();
        MatrixXd eigenvectors = eigen_solver.eigenvectors();
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time SIM GeneralizedSelfAdjointEigenSolver taken: " << duration.count() << " milliseconds" << endl;

        Phi = temp_phi * eigenvectors;
        // Check converged
        bool c = true;
        start = chrono::high_resolution_clock::now();
        if (metric == "I")
        {
            int i = 0;
            for (i = 0; i < p; i++)
            {
                VectorXd Sv = shifted_S * Phi.col(i);
                VectorXd Mv = M * Phi.col(i);
                double norm_a = (Sv - eigenvalues(i) * Mv).norm();
                double norm_b = Mv.norm();
                double residual = norm_a / norm_b;
                if (residual >= epsilon)
                {
                    c = false;
                    break;
                }
            }
        }
        /*M^-1 Norm*/
        /*LDLT for M*/
        else
        {

            auto start_time1 = std::chrono::high_resolution_clock::now();
            Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double>> M_solver;
            M_solver.compute(M);
            auto end_time1 = std::chrono::high_resolution_clock::now();
            auto duration_time1 = std::chrono::duration_cast<std::chrono::milliseconds>(end_time1 - start_time1);
            cout << "Time  SimplicialLDLT M_solver taken: " << duration_time1.count() << " milliseconds" << endl;
            for (int i = 0; i < p; i++)
            {
                auto start_time = std::chrono::high_resolution_clock::now();
                VectorXd Sv = shifted_S * Phi.col(i);
                VectorXd Mv = M * Phi.col(i);
                VectorXd a = Sv - eigenvalues(i) * Mv;
                VectorXd M_inv_a = M_solver.solve(a);

                double norm_a = sqrt(a.dot(M_inv_a));
                VectorXd M_inv_Sv = M_solver.solve(Sv);
                double norm_b = sqrt(Sv.dot(M_inv_Sv));
                double residual = norm_a / norm_b;
                auto end_time = std::chrono::high_resolution_clock::now();
                auto duration_time = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
                cout << "Time SIM residual taken: " << duration_time.count() << " milliseconds" << endl;

                if (residual >= epsilon)
                {
                    c = false;
                    break;
                }
            }
        }
        end = chrono::high_resolution_clock::now();
        duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "Time SIM converged taken: " << duration.count() << " milliseconds" << endl;
        if (c)
        {
            converged = true;
        }
    }
    eigenvalues = eigenvalues.array() + mu;
    return make_pair(eigenvalues, Phi);
}

void computeDistances(
    const vector<set<pair<int, double>>> &adjacency,
    const int &source,
    vector<double> &distances)
{
    distances[source] = 0;
    std::priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> pq;
    pq.push({0.0, source});
    while (!pq.empty())
    {
        auto [curr_dist, curr_vertex] = pq.top();
        pq.pop();
        if (curr_dist > distances[curr_vertex])
        {
            continue;
        }
        for (const auto &[neighbor, edge_length] : adjacency[curr_vertex])
        {
            double new_dist = curr_dist + edge_length;

            if (new_dist < distances[neighbor])
            {
                distances[neighbor] = new_dist;
                pq.push({new_dist, neighbor});
            }
        }
    }
}

vector<vector<int>> Construct_hierarchy(
    const std::vector<Vector3d> &vertices,
    const std::vector<Vector3i> &faces,
    const int &T, const int &p)
{
    // Adjacency Matrix
    vector<set<pair<int, double>>> adjacency;
    adjacency.resize(vertices.size());
    for (int i = 0; i < faces.size(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            int v1 = faces[i](j);
            int v2 = faces[i]((j + 1) % 3);
            double dist = (vertices[v1] - vertices[v2]).norm();
            adjacency[v1].insert({v2, dist});
            adjacency[v2].insert({v1, dist});
        }
    }
    // Init
    vector<vector<int>> hierarchy;
    hierarchy.resize(T);
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(0, vertices.size() - 1);
    int start_v = dis(gen);
    vector<int> VT = {start_v};
    int n = max((int)(1.5 * p), 1000);
    // For test
    if (n >= vertices.size())
    {
        n = min((int)(1.5 * p), 1000);
        if (n >= vertices.size())
            n = vertices.size();
    }
    double mu = pow((double)(vertices.size()) / (double)n, 1.0 / (double)T);
    //   distance vector
    vector<double>
        distances(adjacency.size(), numeric_limits<double>::infinity());
    computeDistances(adjacency, start_v, distances);
    for (int i = T - 1; i > 0; i--)
    {
        vector<int> V_current = VT;
        for (int j = V_current.size(); j < n; j++)
        {
            // find max
            auto max_dist = max_element(distances.begin(), distances.end());
            int max_index = distance(distances.begin(), max_dist);
            V_current.push_back(max_index);
            // update distances
            computeDistances(adjacency, max_index, distances);
        }
        hierarchy[i] = V_current;
        VT = V_current;
        n = n * mu;
    }
    for (int i = 0; i < vertices.size(); i++)
    {
        hierarchy[0].push_back(i);
    }
    return hierarchy;
}

vector<SparseMatrix<double>> Build_Prolongation(
    const vector<vector<int>> &Hierarchy,
    const vector<Vector3d> &vertices,
    const vector<Vector3i> &faces,
    const double &sigma)
{
    double A = 0;

    vector<SparseMatrix<double>> result;
    result.resize(Hierarchy.size() - 1);
    // area of the surface
    for (int i = 0; i < faces.size(); i++)
    {
        Vector3d v1 = vertices[faces[i](0)];
        Vector3d v2 = vertices[faces[i](1)];
        Vector3d v3 = vertices[faces[i](2)];
        Vector3d e1 = v2 - v1;
        Vector3d e2 = v3 - v2;
        Vector3d e3 = v1 - v3;
        double area = 0.5 * (e1.cross(-e3)).norm();
        A += area;
    }
    // For all level
    for (int i = Hierarchy.size() - 1; i > 0; i--)
    {
        double rho = sqrt((sigma * A) / Hierarchy[i].size() * M_PI);
        SparseMatrix<double> U;
        U.resize(Hierarchy[i - 1].size(), Hierarchy[i].size());
        vector<Triplet<double>> tripletList;
        // Normalization
        vector<double> rowsum;
        rowsum.resize(Hierarchy[i - 1].size());
        // Each row
        for (int j = 0; j < Hierarchy[i - 1].size(); j++)
        {
            // Each col
            rowsum[j] = 0;
            for (int k = 0; k < Hierarchy[i].size(); k++)
            {
                Eigen::Vector3d vi = vertices[Hierarchy[i - 1][j]];
                Eigen::Vector3d vj = vertices[Hierarchy[i][k]];
                double dist = (vi - vj).norm();
                if (dist <= rho)
                {
                    double temp_u = 1 - dist / rho;
                    rowsum[j] += temp_u;
                    tripletList.push_back(Eigen::Triplet<double>(j, k, temp_u));
                }
            }
        }
        for (auto &iter : tripletList)
        {
            int row = iter.row();
            int col = iter.col();
            double value = iter.value();
            iter = Triplet<double>(row, col, value / rowsum[row]);
        }
        U.setFromTriplets(tripletList.begin(), tripletList.end());
        result[i - 1] = U;
    }
    return result;
}

vector<SparseMatrix<double>> Build_Prolongation_rigid(
    const vector<vector<int>> &Hierarchy,
    const vector<Vector3d> &vertices,
    const vector<Vector3i> &faces,
    const double &sigma)
{
    // Build prolongation for Laplace
    vector<SparseMatrix<double>> U;
    U = Build_Prolongation(Hierarchy, vertices, faces, sigma);
    for (int i = 0; i < U.size(); i++)
    {
        SparseMatrix<double, RowMajor> u = U[i];
        vector<Triplet<double>> tripletList;
        for (int j = 0; j < u.rows(); j++)
        {
            // Get each row nozero element
            VectorXd u_row(u.row(j).nonZeros());
            MatrixXd C(4, u_row.rows());
            // Build d
            int indexd = Hierarchy[i][j];
            Vector4d d(1, vertices[indexd](0), vertices[indexd](1), vertices[indexd](2));
            int idx = 0;
            for (SparseMatrix<double, RowMajor>::InnerIterator it(u, j); it; ++it)
            {
                u_row(idx) = it.value();
                // Build C Matrix C*UT = d
                int index = Hierarchy[i + 1][it.col()];
                C(0, idx) = 1;
                C(1, idx) = vertices[index](0);
                C(2, idx) = vertices[index](1);
                C(3, idx) = vertices[index](2);
                idx++;
            }
            /*
                Solve C * CT lambda = C * UT - d
            */
            MatrixXd A = C * C.transpose();
            MatrixXd b = C * u_row - d;
            LDLT<MatrixXd> ldlt(A);
            Vector4d lambda = ldlt.solve(b);
            /*
                Build new U  = U - lambdaT*C
            */
            VectorXd n_u = u_row - C.transpose() * lambda;
            idx = 0;
            for (SparseMatrix<double, RowMajor>::InnerIterator it(u, j); it; ++it)
            {
                tripletList.push_back(Triplet<double>(j * 3, it.col() * 3, n_u(idx)));
                tripletList.push_back(Triplet<double>(j * 3 + 1, it.col() * 3 + 1, n_u(idx)));
                tripletList.push_back(Triplet<double>(j * 3 + 2, it.col() * 3 + 2, n_u(idx)));
                idx++;
            }
        }
        SparseMatrix<double> new_u(U[i].rows() * 3, U[i].cols() * 3);
        new_u.setFromTriplets(tripletList.begin(), tripletList.end());
    }
    return U;
}

pair<VectorXd, MatrixXd> HSIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    const int &p,
    const int &T,
    const double &epsilon,
    const vector<SparseMatrix<double>> &U,
    const string &metric,
    const int &zero)
{
    vector<SparseMatrix<double>> ST;
    vector<SparseMatrix<double>> MT;
    ST.resize(T);
    MT.resize(T);
    ST[0] = S;
    MT[0] = M;
    // Build stiffness matrix and mass matrix for all level
    auto start = chrono::high_resolution_clock::now();
    for (int i = 1; i <= T - 1; i++)
    {
        ST[i] = U[i - 1].transpose() * ST[i - 1] * U[i - 1];
        MT[i] = U[i - 1].transpose() * MT[i - 1] * U[i - 1];
    }
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time1 taken: " << duration.count() << " milliseconds" << endl;

    int q = max((int)(1.5 * p), p + 8);
    MatrixXd DS = MatrixXd(ST[T - 1]);
    MatrixXd DM = MatrixXd(MT[T - 1]);

    // Compute first q eigenpairs
    start = chrono::high_resolution_clock::now();
    DenseSymMatProd<double> opS(DS);
    DenseCholesky<double> opM(DM);
    q = min(q, (int)DS.rows() - 1);
    int nev = min(q * 2, (int)DS.rows());
    SymGEigsSolver<DenseSymMatProd<double>, DenseCholesky<double>, GEigsMode::Cholesky> eigs(opS, opM, q, nev);
    eigs.init();
    int nconv = eigs.compute(SortRule::SmallestAlge);
    VectorXd eigenvalues = eigs.eigenvalues();
    MatrixXd eigenvectors = eigs.eigenvectors();
    eigenvalues.reverseInPlace();
    eigenvectors.rowwise().reverseInPlace();
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time2 dense spectra taken: " << duration.count() << " milliseconds" << endl;

    start = chrono::high_resolution_clock::now();
    for (int i = T - 2; i >= 0; i--)
    {
        cout << i << endl;
        cout << eigenvalues << endl;
        MatrixXd temp_eigenvectors = U[i] * eigenvectors;
        double mu = -eigenvalues(zero);
        pair<VectorXd, MatrixXd> result = SIM(ST[i], MT[i], temp_eigenvectors, p, epsilon, mu, metric);
        eigenvalues = result.first;
        eigenvectors = result.second;
    }
    end = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Time3 iteration taken: " << duration.count() << " milliseconds" << endl;
    return make_pair(eigenvalues, eigenvectors);
}
