#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>
#include <queue>
#include <random>
#include <cmath>
using namespace Eigen;
using namespace std;

std::pair<Eigen::VectorXd, Eigen::MatrixXd> SIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    MatrixXd &Phi,
    int p,
    double epsilon,
    double mu);
void computeDistances(
    const std::vector<std::set<std::pair<int, double>>> &adjacency,
    int source,
    std::vector<double> &distances);
std::vector<vector<int>> Construct_hierarchy(std::vector<Eigen::Vector3d> vertices, std::vector<Eigen::Vector3i> faces, int T, int p);

std::vector<SparseMatrix<double>> Build_Prolongation(
    vector<vector<int>> Hierarchy,
    std::vector<Eigen::Vector3d> vertices,
    std::vector<Eigen::Vector3i> faces,
    double sigma);

std::pair<Eigen::VectorXd, Eigen::MatrixXd> HSIM(
    const SparseMatrix<double> &S,
    const SparseMatrix<double> &M,
    int p,
    int T,
    double epsilon,
    std::vector<SparseMatrix<double>> U);