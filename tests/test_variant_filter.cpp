/******************************************************************************
*                                                                             *
*  Copyright © 2024-2025 -- IGH / LIRMM / CNRS / UM                           *
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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
#include <vcfpp.h>
#pragma GCC diagnostic pop

#include <file_reader.h>
#include <kim_exception.h>
#include <variant_filter_driver.h>

#include <iostream>
#include <set>
#include <string>
#ifdef NDEBUG
#  undef NDEBUG
#endif
#include <cassert>

using namespace std;
using namespace kim;

void testFilter(const string &fname, const string &filter, set<string> to_keep) {

  cout << "Testing VariantFilter on file '" << fname << "'"
       << " using filter '" << filter << "'" << endl;
  vcfpp::BcfReader vcf(fname);
  // const vcfpp::BcfHeader &hdr = vcf.getHeader();
  vcfpp::BcfRecord variant(vcf.header); // construct a variant record
  while (vcf.getNextVariant(variant)) {
    string v_str = variant.asString();
    v_str.erase(v_str.find_last_not_of('\n') + 1);
    cout << "  Variant: '" << v_str << "': ";
    VariantFilterDriver drv(variant);
    drv.debug_scanning = false;
    drv.debug_parsing = false;
    set<string>::const_iterator it = to_keep.find(variant.ID());
    if (drv.apply(filter)) {
      cout << "\033[32;1mPASS\033[0m" << endl;
      if (it != to_keep.end()) {
        to_keep.erase(it);
      } else {
        Exception e;
        e << "Variant '" << variant.ID() << "'"
          << " should not have passed filter '" << filter << "'"
          << " but it does.";
        throw e;
      }
    } else {
      cout << "\033[31;1mSKIP\033[0m" << endl;
      if (it != to_keep.end()) {
        Exception e;
        e << "Variant '" << variant.ID() << "'"
          << " should have passed filter '" << filter << "'"
          << " but it doesn't.";
        throw e;
      }
    }
  }
  if (!to_keep.empty()) {
    Exception e;
    e << "All expected variants seems not to have been processed:";
    for (const auto &v: to_keep) {
      e << "\n- '" << v << "'";
    }
    throw e;
  }
}


struct FilterAndResult {
  string filter;
  set<string> to_keep;
};

int main() {

  const FilterAndResult filter_and_res[] = {
    { "Chrom~'.*'",
      {
        "rs1", "rs2", "rs3", "rs4", "rs5", "rs6", "rs7", "rs8", "rs9", "rs10",
        "rs11", "rs12", "rs13", "rs14", "rs15", "rs16", "rs17", "rs18", "rs19", "rs20",
        "rs21", "rs22", "rs23", "rs24", "rs25", "rs26", "rs27", "rs28", "rs29", "rs30",
        "rs31", "rs32", "rs33", "rs34", "rs35", "rs36", "rs37", "rs38", "rs39", "rs40"
      }
    },
    { "Chrom=\"fake\"", {} },
    { "Chrom='B'",
      { "rs11", "rs12", "rs13", "rs14", "rs15", "rs16", "rs17", "rs18", "rs19", "rs20" }
    },
    { "Pos>123",
      {
        "rs5", "rs6", "rs7", "rs8", "rs9", "rs10",
        "rs17", "rs18", "rs19", "rs20",
        "rs27", "rs28", "rs29", "rs30",
        "rs37", "rs38", "rs39", "rs40"
      }
    },
    { "Pos=48", { "rs33" } },
    { "ID~\"rs1.\"",
      { "rs10", "rs11", "rs12", "rs13", "rs14", "rs15", "rs16", "rs17", "rs18", "rs19" }
    },
    { "snp=1",
      {
        "rs1", "rs6", "rs7", "rs10", "rs11", "rs13", "rs17",
        "rs21", "rs23", "rs26", "rs30", "rs33", "rs35", "rs38"
      }
    },
    { "info:RS>=35",
      { "rs35", "rs36", "rs37", "rs38", "rs39", "rs40" }
    },
    { "info:VC='INS'",
      { "rs3", "rs4", "rs5", "rs9", "rs14", "rs20", "rs22", "rs27" }
    },
    { "info:R5=true", { "rs4" } },
    { "snp=false",
      {
        "rs2", "rs3", "rs4", "rs5", "rs8", "rs9",
        "rs12", "rs14", "rs15", "rs16", "rs18", "rs19", "rs20",
        "rs22", "rs24", "rs25", "rs27", "rs28", "rs29", "rs31", "rs32",
        "rs34", "rs36", "rs37", "rs39", "rs40"
      }
    },
    { "msnp=true",
      { "rs2", "rs8", "rs15", "rs25", "rs28", "rs31", "rs34" }
    },
    { "InDel=true",
      {
        "rs3", "rs4", "rs5", "rs9", "rs12", "rs14", "rs16", "rs18", "rs19", "rs20",
        "rs22", "rs24", "rs27", "rs29", "rs32", "rs36", "rs37", "rs39", "rs40"
      }
    },
    { "SV=true", {} },
    { "SNP=1 or MSNP=1",
      {
        "rs1", "rs2", "rs6", "rs7", "rs8", "rs10", "rs11", "rs13", "rs15", "rs17",
        "rs21", "rs23", "rs25", "rs26", "rs28", "rs30",
        "rs31", "rs33", "rs34", "rs35", "rs38"
      }
    },
    { "Chrom='D' and ( pos < 50 or pos > 100)",
      { "rs31", "rs32", "rs33", "rs36", "rs37", "rs38", "rs39", "rs40" }
    },
    { "not Chrom~'[A-C]' and not(pos >= 50 and pos <= 100)",
      { "rs31", "rs32", "rs33", "rs36", "rs37", "rs38", "rs39", "rs40" }
    },
    { "not (Chrom~'[A-C]' or pos >= 50 and pos <= 100)",
      { "rs31", "rs32", "rs33", "rs36", "rs37", "rs38", "rs39", "rs40" }
    },
    { "INFO:VC='INS' and INFO:R5=0 or not INFO:INT=false",
      { "rs1", "rs3", "rs5", "rs9", "rs14", "rs20", "rs22", "rs27" }
    }
  };

  const FilterAndResult bad_filters[] = {
    { "foo", {}},
    { "and Chrom", {}},
    { "Chrom=A and ", {}},
    { "Chrom=A and (true)", {}},
    { "Chrom>=32", {}},
    { "pos>=32s", {}},
    { "(pos>=32", {}},
    { "pos>=(32)", {}},
    { "pos>=32)", {}},
    { "info=foo", {}},
    { "info:VC=foo", {}},
    { "info:\nVC=\n'INS'\nand pos='32'\nor\npos<100", {}},
    { "info:\nVC=\n'INS'\nand pos=3\n22222222\n222\nor\npos<100", {}}
  };


  const char *dirs[] = {
    "./",
    PACKAGE_DATADIR "/",
    /* The following directories are mostly for development purpose */
    "resources/",
    "../resources/",
    "../",
    SRCDIR "/../resources/",
    /* The last array entry must be NULL */
    NULL
  };

  for (const char *f: {"test-variants.vcf", "test-variants.vcf.gz"}) {
    const string fname = FileReader::findFile(f, dirs);

    for (const auto &fr: filter_and_res) {
      testFilter(fname, fr.filter, fr.to_keep);
    }

    for (const auto &fr: bad_filters) {
      try {
        testFilter(fname, fr.filter, fr.to_keep);
      } catch (const VariantFilterDriverException &e) {
        cout << "The following exception is an expected result:" << endl;
        cout << e.what() << endl;
      }
    }
  }

  cout << "That's All, Folk!!!" << endl;
  return 0;

}
