
#include <Eigen/Dense>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include "HSIM.h"
#include "laplace_3d.cpp"

using namespace Eigen;
using namespace std;

bool readOFF(const std::string &filepath, vector<Vector3d> vertices, vector<Vector3i> faces)
{
    std::ifstream file(filepath);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filepath << std::endl;
        return false;
    }

    std::string line;
    std::getline(file, line); // 读取第一行，应该是 "OFF"

    if (line != "OFF")
    {
        std::cerr << "Error: File is not in OFF format." << std::endl;
        return false;
    }

    // 读取顶点数、面数和边数
    int numVertices = 0, numFaces = 0, numEdges = 0;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
            continue; // 跳过空行和注释
        std::istringstream iss(line);
        iss >> numVertices >> numFaces >> numEdges;
        break;
    }

    if (numVertices <= 0 || numFaces <= 0)
    {
        std::cerr << "Error: Invalid number of vertices or faces." << std::endl;
        return false;
    }

    // 读取顶点
    vertices.resize(numVertices, 3);
    for (int i = 0; i < numVertices; ++i)
    {
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue; // 跳过空行和注释
            std::istringstream iss(line);
            iss >> vertices(i, 0) >> vertices(i, 1) >> vertices(i, 2);
            break;
        }
    }

    // 读取面
    faces.resize(numFaces, 3);
    for (int i = 0; i < numFaces; ++i)
    {
        while (std::getline(file, line))
        {
            if (line.empty() || line[0] == '#')
                continue; // 跳过空行和注释
            std::istringstream iss(line);
            int numVerticesInFace;
            iss >> numVerticesInFace;
            if (numVerticesInFace != 3)
            {
                std::cerr << "Error: Only triangular faces are supported." << std::endl;
                return false;
            }
            iss >> faces(i, 0) >> faces(i, 1) >> faces(i, 2);
            break;
        }
    }

    file.close();
    return true;
}

int main()
{
    // 创建一个单位正方体的顶点
    vector<Vector3d> vertices = {
        Vector3d(-0.5, -0.5, -0.5), // 0: 左下后
        Vector3d(0.5, -0.5, -0.5),  // 1: 右下后
        Vector3d(0.5, 0.5, -0.5),   // 2: 右上后
        Vector3d(-0.5, 0.5, -0.5),  // 3: 左上后
        Vector3d(-0.5, -0.5, 0.5),  // 4: 左下前
        Vector3d(0.5, -0.5, 0.5),   // 5: 右下前
        Vector3d(0.5, 0.5, 0.5),    // 6: 右上前
        Vector3d(-0.5, 0.5, 0.5)    // 7: 左上前
    };

    // 创建正方体的三角形面片（每个面由两个三角形组成）
    vector<Vector3i> faces = {
        // 前面
        Vector3i(4, 5, 6),
        Vector3i(4, 6, 7),
        // 后面
        Vector3i(1, 0, 3),
        Vector3i(1, 3, 2),
        // 左面
        Vector3i(0, 4, 7),
        Vector3i(0, 7, 3),
        // 右面
        Vector3i(5, 1, 2),
        Vector3i(5, 2, 6),
        // 上面
        Vector3i(7, 6, 2),
        Vector3i(7, 2, 3),
        // 下面
        Vector3i(0, 1, 5),
        Vector3i(0, 5, 4)};

    Construct_hierarchy(vertices, faces, 3, 5);

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