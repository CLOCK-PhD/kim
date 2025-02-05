/******************************************************************************
*                                                                             *
*  Copyright © 2023-2025 -- IGH / LIRMM / CNRS / UM                           *
*                           (Institut de Génétique Humaine /                  *
*                           Laboratoire d'Informatique, de Robotique et de    *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique /    *
*                           Université de Montpellier)                        *
*                                                                             *
*                                                                             *
*  Auteurs/Authors:                                                           *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - William RITCHIE  <william.ritchie@igh.cnrs.fr>                         *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Programmeurs/Programmers:                                                  *
*    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              *
*    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            *
*                                                                             *
*                                                                             *
*  Contact:                                                                   *
*    - KIM list         <kim@lirmm.fr>                                        *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce logiciel  est un  programme informatique  permettant  d'identifier des  *
*  variations génomiques à partir de données brutes de séquençage.            *
*                                                                             *
*  Ce logiciel est régi par la  licence CeCILL  soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés au   chargement, à   l'utilisation, à  la modification  et/ou au  *
*  développement et   à la reproduction du  logiciel par l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs sont   donc invités   à charger   et tester  l'adéquation du  *
*  logiciel à   leurs besoins  dans des  conditions permettant  d'assurer la  *
*  sécurité de leurs systèmes et ou de leurs données et, plus  généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait que   vous puissiez accéder à cet  en-tête signifie que vous avez  *
*  pris connaissance de la licence CeCILL,   et que vous en avez accepté les  *
*  termes.                                                                    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This software is a computer program whose purpose is to indentify genomic  *
*  from raw sequencing data.                                                  *
*                                                                             *
*  This software is governed by the CeCILL license under French law and       *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL     *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL license and that you accept its terms.             *
*                                                                             *
******************************************************************************/

#include <sort_helper.h>

#include <string>
#include <iostream>
#include <iterator>
#include <libgen.h>

using namespace std;

int main(int argc, char **argv) {

  vector<char> ref;
  vector<string> to_sort;
  vector<bool> pair_length;

  if (argc < 2) {
    const string example = " This is a Demo for the SortHelper template";
    cerr << endl
         << "usage: " << basename(argv[0]) << " <string_1> <string_2> ... <string_n>" << endl
         << endl
         << "This test/demonstration program of the sort_helper template uses the"
         << " first character of each string as a reference to sort the whole set of"
         << " strings given on the command line." << endl
         << endl
         << "Let run the command with the following arguments:" << endl
         << endl
         << argv[0] << example << endl
         << endl;
    size_t start;
    size_t end = 0;
    while ((start = example.find_first_not_of(" ", end)) != string::npos) {
      end = example.find(" ", start);
      string w = example.substr(start, end - start);
      to_sort.push_back(w);
    }
  }

  for (int i = 1; i < argc; ++i) {
    to_sort.push_back(argv[i]);
  }
  for (size_t i = 0; i < to_sort.size(); ++i) {
    ref.push_back(to_sort[i][0]);
    pair_length.push_back(!(to_sort[i].size() & 1));
  }

  cout << "Entering the demo program with the following parameters:" << endl;
  cout << "- ref is:\t\t";
  copy(ref.begin(), ref.end(), ostream_iterator<char>(cout, "\t"));
  cout << endl;
  cout << "- to_sort is:\t\t";
  copy(to_sort.begin(), to_sort.end(), ostream_iterator<string>(cout, "\t"));
  cout << endl;
  cout << "- pair_length is:\t";
  copy(pair_length.begin(), pair_length.end(), ostream_iterator<bool>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;

  cout << "SortHelper initialization:" << endl;
  kim::SortHelper<char> helper(ref);
  cout << "- ref is:\t\t";
  copy(helper.reference().begin(), helper.reference().end(), ostream_iterator<char>(cout, "\t"));
  cout << endl;
  cout << "- permutation is:\t";
  copy(helper.permutation().begin(), helper.permutation().end(), ostream_iterator<size_t>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;


  cout << "Sorting the to_sort and pair_length arrays with SortHelper:" << endl;
  helper.sort<string>(to_sort);
  helper.sort<bool, vector<bool>, vector<bool>::reference>(pair_length);
  cout << "- to_sort is now:\t";
  copy(to_sort.begin(), to_sort.end(), ostream_iterator<string>(cout, "\t"));
  cout << endl;
  cout << "- pair_length is now:\t";
  copy(pair_length.begin(), pair_length.end(), ostream_iterator<bool>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;
  struct first_letter_cmp {
    static bool cmp(const string &s1, const string &s2) {
      return s1[0] < s2[0];
    }
  };
  if (!is_sorted(to_sort.begin(), to_sort.end(), first_letter_cmp::cmp)) {
    cerr << "The to_sort array is not sorted while it should be!" << endl;
    return 1;
  }

  cout << "Sorting ref with SortHelper:" << endl;
  helper.sort<char>(ref);
  cout << "- ref is now :\t\t";
  copy(ref.begin(), ref.end(), ostream_iterator<char>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;
  if (!is_sorted(ref.begin(), ref.end())) {
    cerr << "The ref array is not sorted while it should be!" << endl;
    return 1;
  }

  cout << "Sorting the ref implies to reset the SortHelper:" << endl;
  helper.reset();
  cout << "- ref is:\t\t";
  copy(helper.reference().begin(), helper.reference().end(), ostream_iterator<char>(cout, "\t"));
  cout << endl;
  cout << "- permutation is:\t";
  copy(helper.permutation().begin(), helper.permutation().end(), ostream_iterator<size_t>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;
  if (!is_sorted(helper.permutation().begin(), helper.permutation().end())) {
    cerr << "The permutation array is not sorted while it should be!" << endl;
    return 1;
  }

  cout << "Sorting to_sort and pair_length with the SortHelper having the identity permutation works:" << endl;
  helper.sort<string>(to_sort);
  cout << "- to_sort is now:\t";
  copy(to_sort.begin(), to_sort.end(), ostream_iterator<string>(cout, "\t"));
  cout << endl;
  cout << "- pair_length is now:\t";
  copy(pair_length.begin(), pair_length.end(), ostream_iterator<bool>(cout, "\t"));
  cout << endl;
  cout << "==========================================" << endl;
  if (!is_sorted(to_sort.begin(), to_sort.end(), first_letter_cmp::cmp)) {
    cerr << "The to_sort array is not sorted while it still should be!" << endl;
    return 1;
  }

  cout << "That's All, Folk!!!" << endl;
  return 0;
}
