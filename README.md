# HSIM  Command Line

## Usage of HSIM Eigen_solver

```bash
$ ./HSIM --num N -M path_of_M -K path_of_K -U path_of_U0 path_of U_1 ... --epsilon eps --metric I or Minv --zero n_zero
```

- `-num N` : Number of eigen pairs need to get.
- `-M path_of_M` :  The path of the mass matrix.
- `-K path_of_K` : The path of the Stiffness matrix.
- `-U path_of_U0 path_of U_1 ...`:The path of prolongation matrixes.
- `--epsilon eps`:Convergence factor.
- `--metric I or Minv`:Form of metric (I or Minv).
- `--zero n_zero`: The number of zero eigenvalues.

```bash
$ ./HSIM -h
Usage: ./HSIM [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --num INT                   Number of Eigenpairs
  -M TEXT:FILE                Path of the mass matrix
  -K TEXT:FILE                Path of the stiffness matrix
  -U TEXT:FILE ...            Path of prolongation matrix
  --metric TEXT:{I,Minv}      Form of metric (I or Minv)
  --epsilon FLOAT             Convergence factor
  --zero INT                  The number of zero eigenvalues
```

## Usage of `gen_U`

```bash
 $./gen_U --mesh path_mesh -t path_U --div layer --num num_e --sigma sigma
```

- `--mesh path_mesh`:Default .vtk file, format requirements unstructred grid,binary.
- `-t path_U` :Location of the storage U.
- `--div layer` : The number of layers for hierarchy.
- `--num num_e` : The number of eigenpairs.
- `--sigma sigma` :Control factor of hierarchy.Optional, default value is 7.

```bash
$ ./gen_U -h
Usage: ./gen_U [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --mesh TEXT:FILE            Path of the mesh file(.vtk,unstructred grid,binary).
  -t TEXT:DIR                 Location of the storage U.
  --div INT                   The number of layers for hierarchy.
  --num INT                   The number of eigenpairs.
  --sigma FLOAT               Control factor of hierarchy.
```

