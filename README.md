# KIM: k-mer based genomic variation identification in raw genomic sequencing data

[TOC]

[^ignore TOC]: @tableofcontents

KIM stand for k-mer Identification Metric.
It is an alignment-free software which allows to predict genomic variants
from raw sequencing data.

Documentation
-------------

Usage: `kim [options] <file> [<file> ...]`


Installation
------------

### Requirements

KIM requires:

* A C++ compiler such as `g++` or `clang`.
* A lexical anlayser compiler such as `flex`.
* A grammatical anlayser compiler such as `bison`.
* The `zlib` development files (on Debian/Ubuntu it corresponds to the package `zlib1g-dev`).
* `Doxygen`, `cloc` and `ohcount` (recommanded but not mandatory)


### Single user installation

To download and install `KIM` into some user local directory (e.g.,
`${HOME}/local_install`) , use the following commands:

First, it is needed to clone the `kim` repository
```sh
git clone --recurse-submodules https://github.com/CLOCK-PhD/kim.git
```

Once cloned, go to the newly created directory and artificially
restore the relative order of creation/modification dates for some
files. Indeed, creation dates and last modification dates are not
preserved by the `git clone` operation, and quite often it leads to an
infinite loop or an error during the built.

```sh
cd kim
./config/fix-timestamp.sh
```

Now, run the `configure` script, build the program and install it.
```sh
./configure --prefix=${HOME}/local_install
make
make install
```

As an alternative, in order to get built files to be in a separated
dedicated directory you also can run the following (instead of the previous
`./configure` command)
```sh
mkdir build
cd build
../configure --prefix=${HOME}/local_install
```

If you obtain an error as below:
```sh
prompt$ ./configure
Checking for git-version-gen... config/git-version-gen
...
=== configuring in external/htslib (/home/doccy/Work/kim/external/htslib)
configure: running /bin/bash ./configure --disable-option-checking '--prefix=/usr/local'  --cache-file=/dev/null --srcdir=.
checking for gcc... gcc
checking whether the C compiler works... yes
checking for C compiler default output file name... a.out
checking for suffix of executables...
checking whether we are cross compiling... no
checking for suffix of object files... o
checking whether the compiler supports GNU C... yes
checking whether gcc accepts -g... yes
checking for gcc option to enable C11 features... none needed
checking for ranlib... ranlib
checking for grep that handles long lines and -e... /usr/bin/grep
checking for C compiler warning flags... -Wall
checking for pkg-config... /usr/bin/pkg-config
checking pkg-config is at least version 0.9.0... yes
checking for gcc option to enable large file support... none needed
checking shared library type for unknown-Linux... plain .so
checking whether the compiler accepts -fvisibility=hidden... yes
checking build system type... x86_64-pc-linux-gnu
checking host system type... Invalid configuration `unknown-Linux': machine `unknown-unknown' not recognized
configure: error: /bin/bash ./htscodecs/config.sub unknown-Linux failed
configure: error: ./configure failed for external/htslib
```

Then you may try to set the `--host` flag according to your current
OS:
```sh
./configure --host=x86_64-linux --prefix=${HOME}/local_install
```


#### Uninstall

To remove `KIM` from your system use the following commands:

```sh
make uninstall
```

### System installation

  To download and install `KIM` globally on your system, then you have to
  follow the same procedure as described above, but simply remove the
  option passed to the configuration script and run the `make install`
  command as superuser.


Getting Started
---------------

You can test the installation using the following command.

```sh
cd ~/local_install/bin/kim --version
```

### Output example
```sh
TODO
```

About KIM
---------

### Bug Reporting

While we use an extensive set of unit tests and test coverage tools
you might still find bugs in the library. We encourage you to report
any problems with the library via the
[github issue tracking system](https://github.com/CLOCK-PhD/kim/issues)
of the project.

### Licensing

Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM
(Institut de Génétique Humaine /
Laboratoire d'Informatique, de Robotique et de Microélectronique de Montpellier /
Centre National de la Recherche Scientifique /
Université de Montpellier)

-------------------------------------------------------------------------

Ce logiciel  est un  programme informatique  permettant  d'identifier des
variations génomiques à partir de données brutes de séquençage.

Ce logiciel est régi par la licence CeCILL soumise au droit français
et respectant les principes de diffusion des logiciels libres.  Vous
pouvez utiliser, modifier et/ou redistribuer ce programme sous les
conditions de la licence CeCILL telle que diffusée par le CEA, le CNRS
et l'INRIA sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de
copie, de modification et de redistribution accordés par cette
licence, il n'est offert aux utilisateurs qu'une garantie limitée.
Pour les mêmes raisons, seule une responsabilité restreinte pèse sur
l'auteur du programme, le titulaire des droits patrimoniaux et les
concédants successifs.

À cet égard l'attention de l'utilisateur est attirée sur les risques
associés au chargement, à l'utilisation, à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant
donné sa spécificité de logiciel libre, qui peut le rendre complexe à
manipuler et qui le réserve donc à des développeurs et des
professionnels avertis possédant des connaissances informatiques
approfondies.  Les utilisateurs sont donc invités à charger et tester
l'adéquation du logiciel à leurs besoins dans des conditions
permettant d'assurer la sécurité de leurs systèmes et ou de leurs
données et, plus généralement, à l'utiliser et l'exploiter dans les
mêmes conditions de sécurité.

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez
pris connaissance de la licence CeCILL, et que vous en avez accepté
les termes.

-------------------------------------------------------------------------

This software is a computer program whose purpose is to indentify genomic
from raw sequencing data.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software. You can use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and rights to copy,
modify and redistribute granted by the license, users are provided
only with a limited warranty and the software's author, the holder of
the economic rights, and the successive licensors have only limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading, using, modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean that it is complicated to manipulate, and that also
therefore means that it is reserved for developers and experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards
their requirements in conditions enabling the security of their
systems and/or data to be ensured and, more generally, to use and
operate it in the same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.

[Click here to access the full licence](LICENSE.md)

### Auteurs/Authors:

* Rémy COSTA       <remy.costa@igh.cnrs.fr>
* William RITCHIE  <william.ritchie@igh.cnrs.fr>
* Alban MANCHERON  <alban.mancheron@lirmm.fr>


### Programmeurs/Programmers:

* Rémy COSTA       <remy.costa@igh.cnrs.fr>
* Alban MANCHERON  <alban.mancheron@lirmm.fr>

### Contact:

* KIM list         <kim@lirmm.fr>

