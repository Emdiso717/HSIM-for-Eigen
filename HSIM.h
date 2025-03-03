#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <set>
#include <queue>
#include <random>
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
std::vector<set<int>> Construct_hierarchy(std::vector<Eigen::Vector3d> vertices, std::vector<Eigen::Vector3i> faces, int T, int p);