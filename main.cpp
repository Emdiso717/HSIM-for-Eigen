#include "HSIM.h"
#include <CLI/CLI.hpp>
#include "save_spm.hpp"
#include <cstring>
#include "testing.h"
#define DISABLE_COUT std::cout.setstate(std::ios_base::failbit)
#define ENABLE_COUT std::cout.clear()
using namespace std;
using namespace Eigen;

int main(int argc, char **argv)
{
    CLI::App app;
    std::map<int, std::pair<int, std::string>> map;
    int p = 10;
    SparseMatrix<double, Eigen::RowMajor> S_row;
    SparseMatrix<double> S;
    string S_path, M_path;
    SparseMatrix<double, Eigen::RowMajor> M_row;
    SparseMatrix<double> M;
    vector<SparseMatrix<double>> U;
    vector<string> U_paths;
    string metric = "I";
    double epsilon = 1e-6;
    int zero = 6;
    app.add_option("--num", p, "Number of Eigenpairs");
    app.add_option("-M", M_path, "Path of the mass matrix")->check(CLI::ExistingFile);
    app.add_option("-K", S_path, "Path of the stiffness matrix")->check(CLI::ExistingFile);
    app.add_option("-U", U_paths, "Path of prolongation matrix")->check(CLI::ExistingFile);
    app.add_option("--metric", metric, "Form of metric (I or Minv)")->check(CLI::IsMember{std::vector<std::string>{"I", "Minv"}});
    app.add_option("--epsilon", epsilon, "Convergence factor");
    app.add_option("--zero", zero, "The number of zero eigenvalues");
    CLI11_PARSE(app, argc, argv);
    zcy::io::read_spm(M_path.c_str(), M_row);
    zcy::io::read_spm(S_path.c_str(), S_row);
    M = M_row;
    S = S_row;
    int layer = 1;
    for (auto &U_path : U_paths)
    {
        SparseMatrix<double> temp_u;
        zcy::io::read_spm(U_path.c_str(), temp_u);
        U.push_back(temp_u);
        layer++;
    }
    DISABLE_COUT;
    //  Eigen solver
    auto start = std::chrono::high_resolution_clock::now();
    pair<VectorXd, MatrixXd> result = HSIM(S, M, p, layer, epsilon, U, metric, zero);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Time HSIM taken: " << duration.count() << " milliseconds" << std::endl;
    ENABLE_COUT;
    cout << "Eigen values:" << endl;
    cout << result.first.transpose() << endl;
}
