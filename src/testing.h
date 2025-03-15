#pragma once
#include <fstream>
#include "laplace_3d.h"
#include "HSIM.h"
#include "save_spm.hpp"
namespace testing
{
    void readOFF(std::vector<Eigen::Vector3d> &vertices, std::vector<Eigen::Vector3i> &faces)
    {
        std::ifstream file("../test/cube_mesh10000.off");
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

    void Eigen_solver()
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> faces;

        readOFF(vertices, faces);
        int p = 10;
        int layer = 3;
        double epsilon = 1e-6;
        LaplaceBeltrami3D k(vertices, faces);
        Eigen::SparseMatrix<double> stiffness;
        Eigen::SparseMatrix<double> mass;
        k.computeMatrices(stiffness, mass);

        // // build hierarchy
        std::vector<std::vector<int>> Hierarchy;
        Hierarchy = Construct_hierarchy(vertices, faces, layer, p);
        for (auto H : Hierarchy)
        {
            std::cout << H.size() << std::endl;
        }
        std::vector<Eigen::SparseMatrix<double>> U;
        U = Build_Prolongation(Hierarchy, vertices, faces, 7);
        std::cout << "Finished Prolongation" << std::endl;

        auto start = std::chrono::high_resolution_clock::now();
        std::pair<Eigen::VectorXd, Eigen::MatrixXd> result = HSIM(stiffness, mass, p, layer, epsilon, U, "I");
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        std::cout << "Time HSIM taken: " << duration.count() << " milliseconds" << std::endl;
        std::cout << "\n特征值:" << std::endl;
        std::cout << result.first.transpose() << std::endl;

        // auto start = std::chrono::high_resolution_clock::now();
        // MatrixXd Phi = MatrixXd::Random(stiffness.cols(), (int)(p + 8));
        // auto result = SIM(stiffness, mass, Phi, p, 0.01, 10);
        // auto end = std::chrono::high_resolution_clock::now();
        // auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        // std::cout << "Time SIM taken: " << duration.count() << " milliseconds" << std::endl;
        // cout << "\n特征值：" << endl;
        // cout << result.first.transpose() << endl;

        // {
        //     stiffness.makeCompressed();
        //     mass.makeCompressed();
        //     MatrixXd M = MatrixXd(mass);
        //     MatrixXd S = MatrixXd(stiffness);
        //     auto start = std::chrono::high_resolution_clock::now();
        //     DenseSymMatProd<double> opS(S);
        //     DenseCholesky<double> opM(M);
        //     SymGEigsSolver<DenseSymMatProd<double>, DenseCholesky<double>, GEigsMode::Cholesky> eigs(opS, opM, p, p * 2);
        //     eigs.init();
        //     int nconv = eigs.compute(SortRule::SmallestAlge);
        //     VectorXd eigenvalues = eigs.eigenvalues();
        //     MatrixXd eigenvectors = eigs.eigenvectors();
        //     eigenvalues.reverseInPlace();
        //     eigenvectors.rowwise().reverseInPlace();
        //     auto end = std::chrono::high_resolution_clock::now();
        //     auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
        //     std::cout << "Time sparse taken: " << duration.count() << " milliseconds" << std::endl;
        //     cout << eigenvalues.transpose() << endl;
        // }
    }

    void Get_Matrix_MSU()
    {
        // Read vertices and faces
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3i> faces;
        readOFF(vertices, faces);
        // Basic parameters
        int p = 10;
        int layer = 3;
        // Make stiffness and mass
        LaplaceBeltrami3D k(vertices, faces);
        Eigen::SparseMatrix<double> stiffness;
        Eigen::SparseMatrix<double> mass;
        k.computeMatrices(stiffness, mass);
        // write to *.txt
        std::string S_path = "../data/S.bin";
        std::string M_path = "../data/M.bin";
        zcy::io::write_spm(S_path.c_str(), stiffness);
        zcy::io::write_spm(M_path.c_str(), mass);
        // build hierarchy
        std::vector<std::vector<int>> Hierarchy;
        Hierarchy = Construct_hierarchy(vertices, faces, layer, p);
        for (auto H : Hierarchy)
        {
            std::cout << H.size() << std::endl;
        }
        // Make prolongation
        std::vector<Eigen::SparseMatrix<double>> U;
        U = Build_Prolongation(Hierarchy, vertices, faces, 7);
        std::vector<std::string> U_paths;
        for (int i = 0; i < layer - 1; i++)
        {
            std::string u_path = "../data/U" + std::to_string(i) + ".bin";
            zcy::io::write_spm(u_path.c_str(), U[i]);
        }
        std::cout << "Finished Prolongation" << std::endl;
    }

}
