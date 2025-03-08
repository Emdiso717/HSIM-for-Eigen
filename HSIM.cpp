#include "HSIM.h"

std::pair<Eigen::VectorXd, Eigen::MatrixXd> SIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    MatrixXd &Phi,
    int p,
    double epsilon,
    double mu)
{
    SparseMatrix<double> shifted_S = S - mu * M;

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(shifted_S);

    bool converged = false;
    Eigen::VectorXd eigenvalues;

    while (!converged)
    {

        MatrixXd temp_Mphi = M * Phi;
        MatrixXd temp_phi = solver.solve(temp_Mphi);

        MatrixXd Reduced_S = temp_phi.transpose() * MatrixXd(shifted_S) * temp_phi;

        MatrixXd Reduced_M = temp_phi.transpose() * MatrixXd(M) * temp_phi;

        GeneralizedSelfAdjointEigenSolver<MatrixXd> eigen_solver(Reduced_S, Reduced_M);
        eigenvalues = eigen_solver.eigenvalues();
        MatrixXd eigenvectors = eigen_solver.eigenvectors();

        Phi = temp_phi * eigenvectors;
        bool c = true;
        for (int i = 0; i < p; i++)
        {
            VectorXd Sv = shifted_S * Phi.col(i);
            VectorXd Mv = M * Phi.col(i);
            VectorXd a = Sv - eigenvalues(i) * Mv;
            SimplicialLDLT<SparseMatrix<double>> M_solver;
            M_solver.compute(M);
            VectorXd M_inv_a = M_solver.solve(a);

            double norm_a = sqrt(a.dot(M_inv_a));
            VectorXd M_inv_Sv = M_solver.solve(Sv);
            double norm_b = sqrt(Sv.dot(M_inv_Sv));
            double residual = norm_a / norm_b;

            if (residual >= epsilon)
            {
                c = false;
                break;
            }
        }
        if (c)
        {
            converged = true;
        }
    }
    eigenvalues = eigenvalues.array() + mu;
    return make_pair(eigenvalues, Phi);
}

void computeDistances(
    const std::vector<std::set<std::pair<int, double>>> &adjacency,
    int source,
    std::vector<double> &distances)
{
    distances[source] = 0;
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> pq;
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

std::vector<vector<int>> Construct_hierarchy(
    std::vector<Eigen::Vector3d> vertices,
    std::vector<Eigen::Vector3i> faces,
    int T, int p)
{
    // Adjacency Matrix
    std::vector<std::set<std::pair<int, double>>> adjacency;
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
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, vertices.size() - 1);
    int start_v = dis(gen);
    std::vector<int> VT = {start_v};
    int n = max((int)(1.5 * p), 1000);
    // int n = 2; //For Test
    double mu = pow((double)(vertices.size()) / (double)n, 1.0 / (double)T);
    //   distance vector
    std::vector<double>
        distances(adjacency.size(), std::numeric_limits<double>::infinity());
    computeDistances(adjacency, start_v, distances);
    for (int i = T - 1; i > 0; i--)
    {
        std::vector<int> V_current = VT;
        for (int j = V_current.size(); j < n; j++)
        {
            // find max
            auto max_dist = max_element(distances.begin(), distances.end());
            int max_index = std::distance(distances.begin(), max_dist);
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

std::vector<SparseMatrix<double>> Build_Prolongation(
    vector<vector<int>> Hierarchy,
    std::vector<Eigen::Vector3d> vertices,
    std::vector<Eigen::Vector3i> faces,
    double sigma)
{
    double A = 0;

    vector<SparseMatrix<double>> result;
    result.resize(Hierarchy.size());
    // area of the surface
    for (int i = 0; i < faces.size(); i++)
    {
        Eigen::Vector3d v1 = vertices[faces[i](0)];
        Eigen::Vector3d v2 = vertices[faces[i](1)];
        Eigen::Vector3d v3 = vertices[faces[i](2)];
        Eigen::Vector3d e1 = v2 - v1;
        Eigen::Vector3d e2 = v3 - v2;
        Eigen::Vector3d e3 = v1 - v3;
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

std::pair<Eigen::VectorXd, Eigen::MatrixXd> HSIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    int p,
    int T,
    double epsilon,
    std::vector<SparseMatrix<double>> U)
{
    vector<SparseMatrix<double>> ST;
    vector<SparseMatrix<double>> MT;
    ST.resize(T);
    MT.resize(T);
    ST[0] = S;
    MT[0] = M;
    // Build stiffness matrix and mass matrix for all level
    for (int i = 1; i <= T - 1; i++)
    {
        ST[i] = U[i - 1].transpose() * ST[i - 1] * U[i - 1];
        MT[i] = U[i - 1].transpose() * MT[i - 1] * U[i - 1];
    }
    double q = max((int)(1.5 * p), p + 8);
    // Compute first q eigenpairs
    MatrixXd DS = MatrixXd(ST[T - 1]);
    MatrixXd DM = MatrixXd(MT[T - 1]);
    GeneralizedSelfAdjointEigenSolver<MatrixXd> eigen_solver(DS, DM);
    VectorXd temp_value = eigen_solver.eigenvalues();
    MatrixXd temp_vector = eigen_solver.eigenvectors();
    VectorXd eigenvalues = temp_value.head(q);
    MatrixXd eigenvectors = temp_vector.leftCols(q);

    for (int i = T - 2; i >= 0; i--)
    {
        MatrixXd temp_eigenvectors = U[i] * eigenvectors;
        int j = (p + 9) / 10;
        double mu = eigenvalues(j);
        pair<VectorXd, MatrixXd> result = SIM(ST[i], MT[i], temp_eigenvectors, p, epsilon, mu);
        eigenvalues = result.first;
        eigenvectors = result.second;
    }
    return make_pair(eigenvalues, eigenvectors);
}