# HSIM  Command Line

## Usage

```bash
$ ./HSIM --num N -M path_of_M -K path_of_K -U path_of_U0 path_of U_1 ... --epsilon eps --metric I or Minv
```

- `-num N` : Number of eigen pairs need to get.
- `-M path_of_M` :  The path of the mass matrix.
- `-K path_of_K` : The path of the Stiffness matrix.
- `-U path_of_U0 path_of U_1 ...`:The path of prolongation matrixes.
- `--epsilon eps`:Convergence factor.
- `--metric I or Minv`:Form of metric (I or Minv).

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
```

