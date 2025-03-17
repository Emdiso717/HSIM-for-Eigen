#pragma once
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>
#include <queue>
#include <random>
#include <cmath>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/DenseSymMatProd.h>
#include <Spectra/MatOp/DenseCholesky.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <chrono>
#include <Eigen/CholmodSupport>
#include <fstream>

std::pair<Eigen::VectorXd, Eigen::MatrixXd> SIM(
    const Eigen::SparseMatrix<double> &S,
    const Eigen::SparseMatrix<double> &M,
    Eigen::MatrixXd &Phi,
    const int &p,
    const double &epsilon,
    const double &mu,
    const std::string &metric);

void computeDistances(
    const std::vector<std::set<std::pair<int, double>>> &adjacency,
    const int &source,
    std::vector<double> &distances);

std::vector<std::vector<int>> Construct_hierarchy(
    const std::vector<Eigen::Vector3d> &vertices,
    const std::vector<Eigen::Vector3i> &faces,
    const int &T, const int &p);

std::vector<Eigen::SparseMatrix<double>> Build_Prolongation(
    const std::vector<std::vector<int>> &Hierarchy,
    const std::vector<Eigen::Vector3d> &vertices,
    const std::vector<Eigen::Vector3i> &faces,
    const double &sigma);

Eigen::SparseMatrix<double> kroneckerProduct(const Eigen::SparseMatrix<double> &A, const Eigen::SparseMatrix<double> &B);

std::vector<Eigen::SparseMatrix<double>> Build_Prolongation_rigid(
    const std::vector<std::vector<int>> &Hierarchy,
    const std::vector<Eigen::Vector3d> &vertices,
    const std::vector<Eigen::Vector3i> &faces,
    const double &sigma);

std::pair<Eigen::VectorXd, Eigen::MatrixXd> HSIM(
    const Eigen::SparseMatrix<double> &S,
    const Eigen::SparseMatrix<double> &M,
    const int &p,
    const int &T,
    const double &epsilon,
    const std::vector<Eigen::SparseMatrix<double>> &U,
    const std::string &metric);