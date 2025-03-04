
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
    std::ifstream file("../bunny.off");
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << std::endl;
    }
    std::string line;
    // # Created by boundingmesh
    std::getline(file, line);
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

    // // 创建一个单位正方体的顶点
    // vector<Vector3d> vertices = {
    //     Vector3d(-0.5, -0.5, -0.5), // 0: 左下后
    //     Vector3d(0.5, -0.5, -0.5),  // 1: 右下后
    //     Vector3d(0.5, 0.5, -0.5),   // 2: 右上后
    //     Vector3d(-0.5, 0.5, -0.5),  // 3: 左上后
    //     Vector3d(-0.5, -0.5, 0.5),  // 4: 左下前
    //     Vector3d(0.5, -0.5, 0.5),   // 5: 右下前
    //     Vector3d(0.5, 0.5, 0.5),    // 6: 右上前
    //     Vector3d(-0.5, 0.5, 0.5)    // 7: 左上前
    // };

    // // 创建正方体的三角形面片（每个面由两个三角形组成）
    // vector<Vector3i> faces = {
    //     // 前面
    //     Vector3i(4, 5, 6),
    //     Vector3i(4, 6, 7),
    //     // 后面
    //     Vector3i(1, 0, 3),
    //     Vector3i(1, 3, 2),
    //     // 左面
    //     Vector3i(0, 4, 7),
    //     Vector3i(0, 7, 3),
    //     // 右面
    //     Vector3i(5, 1, 2),
    //     Vector3i(5, 2, 6),
    //     // 上面
    //     Vector3i(7, 6, 2),
    //     Vector3i(7, 2, 3),
    //     // 下面
    //     Vector3i(0, 1, 5),
    //     Vector3i(0, 5, 4)};

    // build hierarchy
    vector<vector<int>> Hierarchy;
    Hierarchy = Construct_hierarchy(vertices, faces, 3, 5);
    vector<SparseMatrix<double>> U;
    U = Build_Prolongation(Hierarchy, vertices, faces, 7);
    // for (int i = 0; i < 3; i++)
    // {
    //     cout << i << ":" << endl;
    //     for (int num : Hierarchy[i])
    //     {
    //         std::cout << num << " ";
    //     }
    //     cout << endl;
    //     cout << Hierarchy[i].size() << endl;
    // }

    // LaplaceBeltrami3D k(vertices, faces);
    // Eigen::SparseMatrix<double> stiffness;
    // Eigen::SparseMatrix<double> mass;
    // k.computeMatrices(stiffness, mass);

    // int p = 5;
    // double epsilon = 1e-2;
    // double mu = -1;

    // GeneralizedSelfAdjointEigenSolver<MatrixXd> eigen_solver(stiffness, mass);
    // VectorXd eigenvalues = eigen_solver.eigenvalues();
    // MatrixXd eigenvectors = eigen_solver.eigenvectors();

    // // MatrixXd Phi = eigenvectors.leftCols(5);
    // MatrixXd Phi = MatrixXd::Random(stiffness.cols(), p);

    // try
    // {
    //     auto result = SIM(stiffness, mass, Phi, p, epsilon, mu);

    //     cout << "\n特征值：" << endl;
    //     cout << result.first.transpose() << endl;

    //     cout << "\n特征向量：" << endl;
    //     cout << result.second << endl;
    // }
    // catch (const std::runtime_error &e)
    // {
    //     cout << "计算出错: " << e.what() << endl;
    // }

    // return 0;
}