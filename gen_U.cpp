#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <CLI/CLI.hpp>
#include "save_spm.hpp"
#include "HSIM.h"
using namespace std;
using namespace Eigen;

bool ReadVTKFile(const string &filename, vector<Vector3d> &vertices, vector<Vector3i> &faces)
{
    vtkSmartPointer<vtkUnstructuredGridReader> reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid = reader->GetOutput();

    vtkSmartPointer<vtkPoints> points = unstructuredGrid->GetPoints();

    for (vtkIdType i = 0; i < points->GetNumberOfPoints(); i++)
    {
        double p[3];
        points->GetPoint(i, p);
        vertices.emplace_back(p[0], p[1], p[2]);
    }

    vtkSmartPointer<vtkCellArray> cells = unstructuredGrid->GetCells();

    vtkIdType npts;
    const vtkIdType *pts;
    cells->InitTraversal();
    while (cells->GetNextCell(npts, pts))
    {
        if (npts == 3)
        {
            faces.emplace_back(pts[0], pts[1], pts[2]);
        }
    }

    return true;
}

int main(int argc, char **argv)
{
    CLI::App app;
    string filename;
    int layer = 3;
    int n_eigen = 10;
    double sigma = 7;
    string target = "../data/";
    app.add_option("--mesh", filename, "Path of the mesh file(.vtk,unstructred grid,binary).")->check(CLI::ExistingFile);
    app.add_option("-t", target, "Location of the storage U. ")->check(CLI::ExistingDirectory);
    app.add_option("--div", layer, "The number of layers for hierarchy.");
    app.add_option("--num", n_eigen, "The number of eigenpairs.");
    app.add_option("--sigma", sigma, "Control factor of hierarchy.");
    CLI11_PARSE(app, argc, argv);
    vector<Vector3d> vertices;
    vector<Vector3i> faces;
    ReadVTKFile(filename, vertices, faces);
    vector<vector<int>> hierarchy = Construct_hierarchy(vertices, faces, layer, n_eigen);
    vector<SparseMatrix<double>> U = Build_Prolongation_rigid(hierarchy, vertices, faces, sigma);
    for (int i = 0; i < layer - 1; i++)
    {
        std::string u_path = target + "U" + std::to_string(i) + ".bin";
        zcy::io::write_spm(u_path.c_str(), U[i]);
    }
    cout << "Already stored U in" << target << "U*.bin" << endl;
    return 0;
}