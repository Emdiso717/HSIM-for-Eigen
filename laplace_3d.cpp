#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <iostream>
#include <cmath>

class LaplaceBeltrami3D
{
private:
    std::vector<Eigen::Vector3d> vertices;
    std::vector<Eigen::Vector3i> faces;

public:
    LaplaceBeltrami3D(const std::vector<Eigen::Vector3d> &vertices,
                      const std::vector<Eigen::Vector3i> &faces)
        : vertices(vertices), faces(faces) {}

    void computeTriangleGeometry(const Eigen::Vector3d &v1,
                                 const Eigen::Vector3d &v2,
                                 const Eigen::Vector3d &v3,
                                 double &area,
                                 std::vector<double> &cotangents)
    {
        Eigen::Vector3d e1 = v2 - v1;
        Eigen::Vector3d e2 = v3 - v2;
        Eigen::Vector3d e3 = v1 - v3;

        double l1_squared = e1.squaredNorm();
        double l2_squared = e2.squaredNorm();
        double l3_squared = e3.squaredNorm();

        area = 0.5 * (e1.cross(-e3)).norm();

        cotangents.resize(3);

        double l1 = std::sqrt(l1_squared);
        double l2 = std::sqrt(l2_squared);
        double l3 = std::sqrt(l3_squared);

        double cos_1 = (-l3_squared + l1_squared + l2_squared) / (2 * l1 * l2);
        double sin_1 = 2 * area / (l1 * l2);
        cotangents[0] = cos_1 / sin_1;

        double cos_2 = (-l1_squared + l2_squared + l3_squared) / (2 * l2 * l3);
        double sin_2 = 2 * area / (l2 * l3);
        cotangents[1] = cos_2 / sin_2;

        double cos_3 = (-l2_squared + l3_squared + l1_squared) / (2 * l3 * l1);
        double sin_3 = 2 * area / (l3 * l1);
        cotangents[2] = cos_3 / sin_3;
    }
    void computeMatrices(Eigen::SparseMatrix<double> &stiffness,
                         Eigen::SparseMatrix<double> &mass)
    {
        int numVertices = vertices.size();
        stiffness.resize(numVertices, numVertices);
        mass.resize(numVertices, numVertices);

        std::vector<Eigen::Triplet<double>> stiffnessTriplets;
        std::vector<Eigen::Triplet<double>> massTriplets;

        stiffnessTriplets.reserve(faces.size() * 9);
        massTriplets.reserve(faces.size() * 9);

        for (const auto &face : faces)
        {
            const Eigen::Vector3d &v1 = vertices[face[0]];
            const Eigen::Vector3d &v2 = vertices[face[1]];
            const Eigen::Vector3d &v3 = vertices[face[2]];

            double area;
            std::vector<double> cotangents;
            computeTriangleGeometry(v1, v2, v3, area, cotangents);

            int indices[3] = {face[0], face[1], face[2]};

            // Stiffness Matrix
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    if (i == j)
                    {
                        double sum_cot = 0.0;
                        for (int k = 0; k < 3; k++)
                        {
                            if (k != i)
                            {
                                sum_cot += cotangents[3 - i - k];
                            }
                        }
                        stiffnessTriplets.push_back(
                            Eigen::Triplet<double>(indices[i], indices[i], sum_cot / 2.0));
                    }
                    else
                    {
                        int k = 3 - i - j;
                        stiffnessTriplets.push_back(
                            Eigen::Triplet<double>(indices[i], indices[j], -cotangents[k] / 2.0));
                    }
                }
            }
            // Mass Matrix
            double coeff = area / 12.0;
            for (int i = 0; i < 3; i++)
            {
                massTriplets.push_back(
                    Eigen::Triplet<double>(indices[i], indices[i], 2.0 * coeff));
                for (int j = i + 1; j < 3; j++)
                {
                    massTriplets.push_back(
                        Eigen::Triplet<double>(indices[i], indices[j], coeff));
                    massTriplets.push_back(
                        Eigen::Triplet<double>(indices[j], indices[i], coeff));
                }
            }
        }
        stiffness.setFromTriplets(stiffnessTriplets.begin(), stiffnessTriplets.end());
        mass.setFromTriplets(massTriplets.begin(), massTriplets.end());
    }
};