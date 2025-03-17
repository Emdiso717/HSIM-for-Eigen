#pragma once
#include <Eigen/SparseCore>
#include <cstdint>
#include <fstream>
#include <iostream>
namespace zcy
{
    namespace io
    {
        template <int Major>
        int write_spm(const char *path, const Eigen::SparseMatrix<double, Major> &A)
        {
            std::ofstream ofs(path, std::ios::binary);
            if (ofs.fail())
            {
                return __LINE__;
            }
            const int64_t mat_size[4] = {A.rows(), A.cols(), A.nonZeros(), Major};
            ofs.write((const char *)&mat_size[0], 4 * sizeof(int64_t));
            // write value, Innerindex, OuterIndex
            ofs.write((const char *)A.valuePtr(), A.nonZeros() * sizeof(double));
            ofs.write((const char *)A.innerIndexPtr(), A.nonZeros() * sizeof(int));
            if (Major == Eigen::ColMajor)
            {
                ofs.write((const char *)A.outerIndexPtr(), (A.cols() + 1) * sizeof(int));
            }
            else
            {
                ofs.write((const char *)A.outerIndexPtr(), (A.rows() + 1) * sizeof(int));
            }
            ofs.close();
            return 0;
        }

        template <int Major>
        int read_spm(const char *path, Eigen::SparseMatrix<double, Major> &A)
        {
            std::ifstream ifs(path, std::ios::binary);
            if (ifs.fail())
            {
                return __LINE__;
            }
            int64_t mat_size[4];
            ifs.read((char *)&mat_size[0], 4 * sizeof(int64_t));
            A.resize(mat_size[0], mat_size[1]);
            Eigen::VectorXd val;
            Eigen::VectorXi idx, offset;
            val.resize(mat_size[2]);
            idx.resize(mat_size[2]);
            if (Major != mat_size[3])
            {
                std::cerr << "major not match. A Major: " << Major << ", file Major: " << mat_size[3] << std::endl;
                return __LINE__;
            }
            if (Major == Eigen::ColMajor)
            {
                offset.resize(mat_size[1] + 1);
            }
            else
            {
                offset.resize(mat_size[0] + 1);
            }
            ifs.read((char *)val.data(), val.size() * sizeof(double));
            ifs.read((char *)idx.data(), idx.size() * sizeof(int));
            ifs.read((char *)offset.data(), offset.size() * sizeof(int));
            A = Eigen::Map<Eigen::SparseMatrix<double, Major>>(mat_size[0], mat_size[1], mat_size[2], offset.data(), idx.data(), val.data());

            ifs.close();
            return 0;
        }
    }
}
