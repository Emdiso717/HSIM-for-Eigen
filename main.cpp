
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
    std::ifstream file("../cube_mesh6000.off");
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

    // std::ofstream file("../test/1000.txt");
    // GeneralizedSelfAdjointEigenSolver<MatrixXd> eigen_solver(stiffness, mass);
    // VectorXd eigenvalues = eigen_solver.eigenvalues();
    // MatrixXd eigenvectors = eigen_solver.eigenvectors();
    // file << eigenvalues << endl;
    // file << eigenvectors << endl;

    // MatrixXd Phi = MatrixXd::Random(stiffness.cols(), p + 100);
    pair<VectorXd, MatrixXd> result = HSIM(stiffness, mass, p, 3, 0.01, U);
    // pair<VectorXd, MatrixXd> result = SIM(stiffness, mass, Phi, p, 0.01, 2);
    cout << "\n特征值：" << endl;
    cout << result.first.transpose() << endl;
}