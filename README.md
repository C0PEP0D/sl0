# sl0 : Simple C++ Lagrangian Objects

[![standard-readme compliant](https://img.shields.io/badge/readme%20style-standard-brightgreen.svg?style=flat-square)](https://github.com/RichardLitt/standard-readme)

A simple C++ library that helps advecting various lagrangian objects into a flow.

This simple C++ library provides passive lagrangian objects that are able to be advected in a flow.
Simple examples are provided and should be enough to explain how to use this library.

Note that a similar library providing active lagrangian objects advection is available: [**sa0**](https://github.com/C0PEP0D/sa0)

This repository contains:

1. The software itself provided as a header only library in the directory [include/sl0](./include/sl0)
2. A few [examples](./examples).

This library currently supports the following objects:

* passive point particles
* spheroids
* groups of such particles
* chains

## Table of Contents

- [Background](#background)
- [Install](#install)
- [License](#license)

## Background

This library has been produced during my PhD thesis and as part as the European Research Council project: [C0PEP0D](https://c0pep0d.github.io/)
This library is used as part of [SHELD0N](https://github.com/C0PEP0D/sheld0n), a lagrangian particle advection software.

## Install

### Dependencies

The dependencies are standard softwares that may already be installed on your system.
If not, you should be able to install these dependencies with your package manager.

* [**CMake** `v?`](https://cmake.org/download/)
* a c++17 compliant compiler, such as [**gcc** `v9`](https://gcc.gnu.org/) or higher must be installed. For Ubuntu users: [ppa](https://launchpad.net/%7Ejonathonf/+archive/ubuntu/gcc?field.series_filter=bionic).
* the [**Threading Building Block Library** `v2018`](https://github.com/ibaned/tbb) or higher must be installed ([this version](https://github.com/wjakob/tbb) that enables is installing it using CMake is advised)

Examples:
* [**Eigen**](https://eigen.tuxfamily.org) must be installed
* [**s0s**](https://github.com/C0PEP0D/s0s) must be installed

Chain example:
* [**m0sh**](https://github.com/C0PEP0D/m0sh) must be installed
* [**p0l**](https://github.com/C0PEP0D/p0l) must be installed

The examples assume the following directory tree structure:
```bash
..
 ├── .
 │   │── sl0
 │   │── s0s
 │   │── (m0sh)
 │   └── (p0l)
 └── thirdparty
     └── eigen
```
One should either install these dependencies accordingly, or adapt their path in the **CMakeList.txt** file of the examples.

### Installing

Start by cloning this repository.

```sh
$ git clone https://gitlab.com/rmonthil/c0pep0d.git
```

### Examples

Running an example:

```bash
$ cd examples/point
$ mkdir build
$ cd build
$ cmake ..
$ make
$ ./point
```

### Updating

A simple pull should be enough.

```sh
$ git pull
```

## Maintainers

Rémi Monthiller - [@rmonthil](https://gitlab.com/rmonthil) - remi.monthiller@gmail.com

## Contributing

Feel free to dive in! [Open an issue](https://github.com/rmonthil/c0pep0d/issues/new) or submit PRs.

## License

[MIT © Centrale Marseille, Rémi MONTHILLER.](./LICENSE)
