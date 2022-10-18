# DisCoTec-combischeme-utilities

Python utilities to create (combination technique schemes) input files for the [DisCoTec](https://github.com/SGpp/DisCoTec) code.

![pytest workflow](https://github.com/SGpp/DisCoTec-combischeme-utilities/actions/workflows/python-package.yml/badge.svg)

## Combination Schemes

The [DisCoTec](https://github.com/SGpp/DisCoTec) code is a black-box framework for PDE solvers that uses the sparse grid combination technique.
For an overview on the topic of Sparse Grids, read [Sparse Grids in a Nutshell](ftp://nozdr.ru/biblio/kolxoz/M/MN/MNd/Garcke%20J.,%20Griebel%20M.%20(eds.)%20Sparse%20grids%20and%20applications%20(Springer,%202013)(ISBN%209783642317026)(O)(290s)_MNd_.pdf#page=68) by Jochen Garcke.
A lot of proofs and assumptions about combination technique schemes are also [summarized here](https://link.springer.com/chapter/10.1007/978-3-319-28262-6_4) by Brendan Harding.

The resulting files will have `.json` structure: a long list of dictionaries that each contains a coefficient "coeff" and a level vector "level".

## Get the code

Only available on GitHub currently. In your shell, clone the repository:

```sh
git clone git@github.com:SGpp/DisCoTec-combischeme-utilities.git
```

or

```sh
git clone https://github.com/SGpp/DisCoTec-combischeme-utilities.git
```

The Python modules required to use all functionality can be found in the [requirements.txt file](requirements.txt)


## Run the tests

Enter the directory you just cloned and execute `pytest` on the test file

```sh
pytest test_combischeme_utils.py
```

You should see that all tests are passing. If not, please [create an issue](https://github.com/SGpp/DisCoTec-combischeme-utilities/issues/new)

## Create combination schemes

To create combination schemes in Python, have a close look at the [test file](test_combischeme_utils.py) -- there are quite a few examples there.

For use with the command line, there are some thin scripts for creating combination scheme files.
For instance, you can initially create a combination scheme from the minimum and maximum levels with

```sh
python3 create_combischeme.py --lmin 1 2 3 4 --lmax 5 6 7 8
```

which will place a `.json` file in your folder that has the combination technique memory requirement in its file name. For the minimum and maximum level parameters above, it will be `scheme_10.34_MiB.json`.
Along the same lines, you can then split this scheme into two parts with

```sh
python3 divide_combischeme_by_level_sum.py --file_name scheme_10.34_MiB.json
```

and the resulting files would have `_split1` and `_split2` in their file names.

If you want to use the static task assignment functionality in DisCoTec, you can also assign schemes to process groups with a fair share of degrees of freedom for each process group

```sh
python3 assign_combischeme_to_groups.py --file_name=scheme_10.34_MiB_split1.json --num_groups=32
```

which will add an entry "group_no" into each level in the `.json` file that tells DisCoTec which process group should run the particular component grid.

For really large regular combination schemes, it can make sense to do it all in one step

```sh
python3 create_large_assigned_schemes.py --lmin 1 1 1 1 1 1 --lmax 18 18 18 18 18 18 --num_groups 15 20
```

because this allows the code to use some optimizations which will make it run a lot faster.

## Current limitations

- always assumes both boundary points on level 0
- only 50/50 splits of divided combination scheme

## Acknowledgements

Some of the code originated in the [EXAHD project](https://ipvs.informatik.uni-stuttgart.de/SGS/EXAHD/) (at the time, in [this repo](https://gitlab.lrz.de/sparse_grids/gene_python_interface_clean))
and was written by Christoph Kowitz, Mario Heene, and Alfredo Parra Hinojosa.

If you use the code academically, please cite it.
