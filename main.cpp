
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "HSIM.h"
#include "laplace_3d.cpp"

using namespace Eigen;
using namespace std;

void readOFF(vector<Vector3d> &vertices, vector<Vector3i> &faces)
{
    std::ifstream file("../cube_mesh10000.off");
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << std::endl;
    }
    std::string line;
    // OFF
    std::getline(file, line);
    int nv = 0, nf = 0, ne = 0;
    std::getline(file, line);
    std::istringstream iss(line);
    iss >> nv >> nf >> ne;
    vertices.resize(nv);
    faces.resize(nf);
    for (int i = 0; i < nv; i++)
    {
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            iss >> vertices[i](0) >> vertices[i](1) >> vertices[i](2);
            break;
        }
    }
    for (int i = 0; i < nf; i++)
    {
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue;
            std::istringstream iss(line);
            int numVerticesInFace;
            iss >> numVerticesInFace;
            iss >> faces[i](0) >> faces[i](1) >> faces[i](2);
            break;
        }
    }
    file.close();
}

int main()
{
    vector<Vector3d> vertices;
    vector<Vector3i> faces;

    readOFF(vertices, faces);
    int p = 10;
    LaplaceBeltrami3D k(vertices, faces);
    Eigen::SparseMatrix<double> stiffness;
    Eigen::SparseMatrix<double> mass;
    k.computeMatrices(stiffness, mass);

    // build hierarchy
    vector<vector<int>> Hierarchy;
    Hierarchy = Construct_hierarchy(vertices, faces, 3, p);
    for (auto H : Hierarchy)
    {
        cout << H.size() << endl;
    }
    vector<SparseMatrix<double>> U;
    U = Build_Prolongation(Hierarchy, vertices, faces, 7);
    cout << "Finished Prolongation" << endl;

    auto start = std::chrono::high_resolution_clock::now();
    pair<VectorXd, MatrixXd> result = HSIM(stiffness, mass, p, 3, 0.01, U);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Time HSIM taken: " << duration.count() << " milliseconds" << endl;
    cout << "\n特征值:" << endl;
    cout << result.first.transpose() << endl;

    // start = std::chrono::high_resolution_clock::now();
    // MatrixXd Phi = MatrixXd::Random(stiffness.cols(), (int)(p + 8));
    // result = SIM(stiffness, mass, Phi, p, 0.01, 2);
    // end = std::chrono::high_resolution_clock::now();
    // duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    // std::cout << "Time SIM taken: " << duration.count() << " milliseconds" << std::endl;
    // cout << "\n特征值：" << endl;
    // cout << result.first.transpose() << endl;

    // {
    //     stiffness.makeCompressed();
    //     mass.makeCompressed();
    //     auto start = std::chrono::high_resolution_clock::now();
    //     SparseSymMatProd<double> opS(stiffness);
    //     SparseCholesky<double> opM(mass);
    //     SymGEigsSolver<SparseSymMatProd<double>, SparseCholesky<double>, GEigsMode::Cholesky> eigs(opS, opM, p, p * 2);
    //     eigs.init();
    //     int nconv = eigs.compute(SortRule::SmallestAlge);
    //     VectorXd eigenvalues = eigs.eigenvalues();
    //     MatrixXd eigenvectors = eigs.eigenvectors();
    //     cout << eigenvalues.transpose() << endl;
    //     auto end = std::chrono::high_resolution_clock::now();
    //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    //     std::cout << "Time sparse taken: " << duration.count() << " milliseconds" << std::endl;
    // }
}