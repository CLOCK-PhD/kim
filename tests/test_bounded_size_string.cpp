/******************************************************************************
*                                                                             *
*  Copyright © 2023      -- IGH / LIRMM / CNRS / UM                           *
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

#include <bounded_size_string.h>

#include <string>
#include <iostream>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

int main() {

  cout << "Testing BoundedSizeString" << endl;
  cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getMaximalSize() == 0);
  assert(BoundedSizeString::getNbInstances() == 0);

  cout << "Setting BoundedSizeString maximal size to 3" << endl;
  cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  BoundedSizeString::setMaximalSize(3);
  assert(BoundedSizeString::getMaximalSize() == 3);
  assert(BoundedSizeString::getNbInstances() == 0);

  {
    cout << "Entering in a scope and creating a BoundedSizeString instance" << endl;
    BoundedSizeString s;
    cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
    cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
    cout << "The BoundedSizeString length is " << s.length() << " and its empty status is " << s.empty() << endl;
    assert(BoundedSizeString::getMaximalSize() == 3);
    assert(BoundedSizeString::getNbInstances() == 1);
    assert(s.length() == 0);
    assert(s.empty());

    cout << "Trying to set BoundedSizeString maximal size to 10" << endl;
    BoundedSizeString::setMaximalSize(10);
    cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
    cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
    assert(BoundedSizeString::getMaximalSize() == 3);
    assert(BoundedSizeString::getNbInstances() == 1);
    cout << "Leaving the scope" << endl;
  }
  cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getMaximalSize() == 3);
  assert(BoundedSizeString::getNbInstances() == 0);

  {
    cout << "Entering in a scope and set BoundedSizeString maximal size to 10" << endl;
    BoundedSizeString::setMaximalSize(10);
    cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
    cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
    assert(BoundedSizeString::getMaximalSize() == 10);
    assert(BoundedSizeString::getNbInstances() == 0);
    cout << "Creating a BoundedSizeString instance" << endl;
    BoundedSizeString s;
    cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
    assert(BoundedSizeString::getNbInstances() == 1);
    cout << "The BoundedSizeString length is " << s.length() << " and its empty status is " << s.empty() << endl;
    assert(s.length() == 0);
    assert(s.empty());

    cout << "Trying to set BoundedSizeString maximal size to 3" << endl;
    BoundedSizeString::setMaximalSize(3);
    cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
    cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
    assert(BoundedSizeString::getMaximalSize() == 10);
    assert(BoundedSizeString::getNbInstances() == 1);
    cout << "Leaving the scope" << endl;
  }
  cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;

  cout << "Set BoundedSizeString maximal size to 5" << endl;
  BoundedSizeString::setMaximalSize(5);
  cout << "The maximal size of BoundedSizeStrings is " << BoundedSizeString::getMaximalSize() << endl;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getMaximalSize() == 5);
  assert(BoundedSizeString::getNbInstances() == 0);

  cout << "Creating a BoundedSizeString with \"ABC\"" << endl;
  BoundedSizeString s1 = "ABC";
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 1);
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  assert(s1.length() == 3);
  assert(!s1.empty());
  cout << "s1 = '" << s1 << "'" << endl;

  cout << "Clearing s1" << endl;
  s1.clear();
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  assert(s1.length() == 0);
  assert(s1.empty());

  cout << "Creating a BoundedSizeString with \"0123456789\"" << endl;
  BoundedSizeString s2 = "0123456789";
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 2);
  cout << "The s2 BoundedSizeString is '" << s2 << "', its length is " << s2.length() << " and its empty status is " << s2.empty() << endl;
  assert(s2.length() == 5);
  assert(!s2.empty());
  for (size_t i = 0; i < s2.length(); ++i) {
    assert(s2[i] == char('0' + i));
  }
  cout << "s2 = '" << s2 << "'" << endl;

  cout << "Swapping s1 and s2 content" << endl;
  swap(s1, s2);
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 2);
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  assert(s1.length() == 5);
  assert(!s1.empty());
  cout << "The s2 BoundedSizeString is '" << s2 << "', its length is " << s2.length() << " and its empty status is " << s2.empty() << endl;
  assert(s2.length() == 0);
  assert(s2.empty());
  cout << "Now, s1 = '" << s1 << "' and s2 = '" << s2 << "'" << endl;

  cout << "Assign s1 to s2" << endl;
  s2 = s1;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 2);
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  assert(s1.length() == 5);
  assert(!s1.empty());
  cout << "The s2 BoundedSizeString is '" << s2 << "', its length is " << s2.length() << " and its empty status is " << s2.empty() << endl;
  assert(s2.length() == s1.length());
  assert(!s2.empty());
  cout << "Now, s1 = '" << s1 << "' and s2 = '" << s2 << "'" << endl;

  cout << "Testing comparison operators between s1 and s2: " << endl
       << "- s1 <  s2: " << (s1 < s2) << endl
       << "- s1 <= s2: " << (s1 <= s2) << endl
       << "- s1 == s2: " << (s1 == s2) << endl
       << "- s1 >= s2: " << (s1 >= s2) << endl
       << "- s1 > s2: " << (s1 > s2) << endl
       << "- s1 != s2: " << (s1 != s2) << endl;
  assert((s1 < s2) == false);
  assert((s1 <= s2) == true);
  assert((s1 == s2) == true);
  assert((s1 >= s2) == true);
  assert((s1 > s2) == false);
  assert((s1 != s2) == false);

  cout << "Setting s2[3] = '\\0'" << endl;
  s2[3] = '\0';
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  assert(s1.length() == 5);
  assert(!s1.empty());
  cout << "The s2 BoundedSizeString is '" << s2 << "', its length is " << s2.length() << " and its empty status is " << s2.empty() << endl;
  assert(s2.length() == 3);
  assert(!s2.empty());
  cout << "Now, s1 = '" << s1 << "' and s2 = '" << s2 << "'" << endl;

  cout << "Testing comparison operators between s1 and s2: " << endl
       << "- s1 <  s2: " << (s1 < s2) << endl
       << "- s1 <= s2: " << (s1 <= s2) << endl
       << "- s1 == s2: " << (s1 == s2) << endl
       << "- s1 >= s2: " << (s1 >= s2) << endl
       << "- s1 > s2: " << (s1 > s2) << endl
       << "- s1 != s2: " << (s1 != s2) << endl;
  assert((s1 < s2) == false);
  assert((s1 <= s2) == false);
  assert((s1 == s2) == false);
  assert((s1 >= s2) == true);
  assert((s1 > s2) == true);
  assert((s1 != s2) == true);

  cout << "Making an std::string using s1.c_str()" << endl;
  string s = s1.c_str();
  cout << "The s std string is '" << s << "'" << endl;
  assert(s == string("01234"));

  cout << "Setting the std::string to \"XYZ\"" << endl;
  s = "XYZ";
  cout << "Creating a BoundedSizeString with teh content of s = '" << s << "'" << endl;
  BoundedSizeString s3 = s;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 3);
  cout << "The s3 BoundedSizeString is '" << s3 << "', its length is " << s3.length() << " and its empty status is " << s3.empty() << endl;
  assert(s3.length() == 3);
  assert(!s3.empty());
  cout << "Testing comparison operators between s1 and s3: " << endl
       << "- s1 <  s3: " << (s1 < s3) << endl
       << "- s1 <= s3: " << (s1 <= s3) << endl
       << "- s1 == s3: " << (s1 == s3) << endl
       << "- s1 >= s3: " << (s1 >= s3) << endl
       << "- s1 > s3: " << (s1 > s3) << endl
       << "- s1 != s3: " << (s1 != s3) << endl;
  assert((s1 < s3) == true);
  assert((s1 <= s3) == true);
  assert((s1 == s3) == false);
  assert((s1 >= s3) == false);
  assert((s1 > s3) == false);
  assert((s1 != s3) == true);

  cout << "Testing s1.compare(s2): " << s1.compare(s2) << endl;
  assert(s1.compare(s2) > 0);
  cout << "Testing s1.compare(s3): " << s1.compare(s3) << endl;
  assert(s1.compare(s3) < 0);
  cout << "Testing s2.compare(s1): " << s2.compare(s1) << endl;
  assert(s2.compare(s1) < 0);
  cout << "Testing s2.compare(s3): " << s2.compare(s3) << endl;
  assert(s2.compare(s3) < 0);
  cout << "Testing s3.compare(s1): " << s3.compare(s1) << endl;
  assert(s3.compare(s1) > 0);
  cout << "Testing s3.compare(s2): " << s3.compare(s2) << endl;
  assert(s3.compare(s2) > 0);

  cout << "Testing s1.reverse_compare(s2): " << s1.reverse_compare(s2) << endl;
  assert(s1.reverse_compare(s2) > 0);
  cout << "Testing s1.reverse_compare(s3): " << s1.reverse_compare(s3) << endl;
  assert(s1.reverse_compare(s3) > 0);
  cout << "Testing s2.reverse_compare(s1): " << s2.reverse_compare(s1) << endl;
  assert(s2.reverse_compare(s1) < 0);
  cout << "Testing s2.reverse_compare(s3): " << s2.reverse_compare(s3) << endl;
  assert(s2.reverse_compare(s3) < 0);
  cout << "Testing s3.reverse_compare(s1): " << s3.reverse_compare(s1) << endl;
  assert(s3.reverse_compare(s1) < 0);
  cout << "Testing s3.reverse_compare(s2): " << s3.reverse_compare(s2) << endl;
  assert(s3.reverse_compare(s2) > 0);

  cout << "Assign s1.c_str() to std string s" << endl;
  s = s1.c_str();
  cout << "The s std string is '" << s << "'" << endl;
  assert(s == string("01234"));

  cout << "Assign std string s to s2" << endl;
  s2 = s;
  cout << "There is " << BoundedSizeString::getNbInstances() << " living instances at now." << endl;
  assert(BoundedSizeString::getNbInstances() == 3);
  cout << "The s1 BoundedSizeString is '" << s1 << "', its length is " << s1.length() << " and its empty status is " << s1.empty() << endl;
  cout << "The s2 BoundedSizeString is '" << s2 << "', its length is " << s2.length() << " and its empty status is " << s2.empty() << endl;
  cout << "The s3 BoundedSizeString is '" << s3 << "', its length is " << s3.length() << " and its empty status is " << s3.empty() << endl;
  assert(s1 == s2);
  assert(s1.compare(s2) == 0);
  assert(s1.reverse_compare(s2) == 0);
  assert(s2.compare(s1) == 0);
  assert(s2.reverse_compare(s1) == 0);

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
