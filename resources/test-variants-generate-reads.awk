#!/usr/bin/awk -f
###############################################################################
#                                                                             #
#  Copyright © 2025      -- IGH / LIRMM / CNRS / UM                           #
#                           (Institut de Génétique Humaine /                  #
#                           Laboratoire d'Informatique, de Robotique et de    #
#                           Microélectronique de Montpellier /                #
#                           Centre National de la Recherche Scientifique /    #
#                           Université de Montpellier)                        #
#                                                                             #
#                                                                             #
#  Auteurs/Authors:                                                           #
#    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              #
#    - William RITCHIE  <william.ritchie@igh.cnrs.fr>                         #
#    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            #
#                                                                             #
#                                                                             #
#  Programmeurs/Programmers:                                                  #
#    - Rémy COSTA       <remy.costa@igh.cnrs.fr>                              #
#    - Alban MANCHERON  <alban.mancheron@lirmm.fr>                            #
#                                                                             #
#                                                                             #
#  Contact:                                                                   #
#    - KIM list         <kim@lirmm.fr>                                        #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce logiciel  est un  programme informatique  permettant  d'identifier des  #
#  variations génomiques à partir de données brutes de séquençage.            #
#                                                                             #
#  Ce logiciel est régi par la  licence CeCILL  soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL telle que diffusée par  le CEA,  le CNRS et l'INRIA sur  #
#  le site "http://www.cecill.info".                                          #
#                                                                             #
#  En contrepartie de l'accessibilité au code source et des droits de copie,  #
#  de modification et de redistribution accordés par cette licence, il n'est  #
#  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  #
#  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  #
#  titulaire des droits patrimoniaux et les concédants successifs.            #
#                                                                             #
#  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  #
#  associés au   chargement, à   l'utilisation, à  la modification  et/ou au  #
#  développement et   à la reproduction du  logiciel par l'utilisateur étant  #
#  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  #
#  manipuler et qui le réserve donc à des développeurs et des professionnels  #
#  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  #
#  utilisateurs sont   donc invités   à charger   et tester  l'adéquation du  #
#  logiciel à   leurs besoins  dans des  conditions permettant  d'assurer la  #
#  sécurité de leurs systèmes et ou de leurs données et, plus  généralement,  #
#  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         #
#                                                                             #
#  Le fait que   vous puissiez accéder à cet  en-tête signifie que vous avez  #
#  pris connaissance de la licence CeCILL,   et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This software is a computer program whose purpose is to indentify genomic  #
#  from raw sequencing data.                                                  #
#                                                                             #
#  This software is governed by the CeCILL license under French law and       #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL     #
#  license as circulated by CEA, CNRS and INRIA at the following URL          #
#  "http://www.cecill.info".                                                  #
#                                                                             #
#  As a counterpart to the access to the source code and rights to copy,      #
#  modify and redistribute granted by the license, users are provided only    #
#  with a limited warranty and the software's author, the holder of the       #
#  economic rights, and the successive licensors have only limited            #
#  liability.                                                                 #
#                                                                             #
#  In this respect, the user's attention is drawn to the risks associated     #
#  with loading, using, modifying and/or developing or reproducing the        #
#  software by the user in light of its specific status of free software,     #
#  that may mean that it is complicated to manipulate, and that also          #
#  therefore means that it is reserved for developers and experienced         #
#  professionals having in-depth computer knowledge. Users are therefore      #
#  encouraged to load and test the software's suitability as regards their    #
#  requirements in conditions enabling the security of their systems and/or   #
#  data to be ensured and, more generally, to use and operate it in the same  #
#  conditions as regards security.                                            #
#                                                                             #
#  The fact that you are presently reading this means that you have had       #
#  knowledge of the CeCILL license and that you accept its terms.             #
#                                                                             #
###############################################################################

function assert(cond, msg) {
  if (!cond) {
    print "ERROR: " msg
    exit 1
  }
}

function rnd_nucl_aux(      r) {
  r = int(rand() * 4);
  switch (r) {
  case 0: return "A"
  case 1: return "C"
  case 2: return "G"
  case 3: return "T"
  default: return "X"
  }
}

function rnd_nucl(orig,      n) {
  do {
    n = rnd_nucl_aux()
  } while (n == toupper(orig))
  return n
}

function reads_compare(i1, v1, i2, v2) {
  if (v1["rank"] < v2["rank"]) return -1
  if (v1["rank"] > v2["rank"]) return 1
  return 0
}

BEGIN {
  FS="\t"
  OFS="\t"
  assert(length(PROCINFO["argv"]) >= 3, "PROCINFO['argv'] array is expected to have at least three elements")
  PROGNAME = PROCINFO["argv"][2]
  if (ARGC != 3) {
    printf("Usage: %s <VCF> <FASTA>\n", PROGNAME);
    print("Where:")
    print("  <VCF> is the VCF file")
    print("  <FASTA> is the fasta reference file")
    print("  Order doesn't matter")
    exit 1
  }
  srand()
  min_l = 40
  max_l = 50
  min_coverage = 5
  max_coverage = 8
  variant_rate = 0.3
  unknown_variant_rate = 0.002
  subst_error_rate = 0.0005
  ins_error_rate = 0.0002
  del_error_rate = 0.0002
  total_error_rate = subst_error_rate + ins_error_rate + del_error_rate
  delete VCF["X"][0]["name"]
  file_type=""
}

BEGINFILE {
  if (FILENAME ~ /.*\.vcf/) {
    file_type = "VCF"
  } else {
    file_type = "FASTA"
  }
}

(file_type == "VCF") && (/^#/) {
  next
}

(file_type == "VCF") {
  VCF[$1][$2]["name"] = $3
  VCF[$1][$2]["ref"] = $4
  split($5, VCF[$1][$2]["alt"], ",")
  next
}

/^>/ {
  if (length(name)) {
    sequence[name]=seq
  }
  name=substr($1,2)
  seq=""
  next
}

{
  seq = seq $0
}

END {

  if (length(name)) {
    sequence[name]=seq
  }

  print "#Chr", "Pos", "Id", "Ref", "Alt" > "/dev/stderr"
  unknown_variant_id = 0

  for (name in sequence) {

    seq = sequence[name]
    n = split(seq, sseq, "")
    if (n < min_l) continue

    read = ""

    for (i in variants) {
      delete variants[i]
    }
    for (i in unknown_variants) {
      delete unknown_variants[i]
    }

    # Selection known and unknown variants for current sequence
    for (i in sseq) {
      nucl = sseq[i]
      if (i in VCF[name]) {
        assert(match(nucl, /^[A-Z]$/), "nucl at position " i " is '" nucl "' and is expected to match /^[A-Z]$/");
        if (rand() <= variant_rate) {
          variants[i]["hdr"] = "variant " VCF[name][i]["name"] " at chrom pos " i
          variants[i]["ref"] = VCF[name][i]["ref"]
          nb_alt = length(VCF[name][i]["alt"])
          nb_alt = int(nb_alt * rand()) + 1
          variants[i]["alt"] = VCF[name][i]["alt"][nb_alt]
          print name, i, VCF[name][i]["name"], variants[i]["ref"], variants[i]["alt"] > "/dev/stderr"
        }
      } else {
        if (rand() <= unknown_variant_rate) {
          unknown_variant[i]["hdr"] = "unknown variant " ++unknown_variant_id " at chrom pos " i
          unknown_variant[i]["ref"] = nucl
          unknown_variant[i]["alt"] = rnd_nucl(nucl)
          print name, i, "unknown" unknown_variant_id, unknown_variant[i]["ref"], unknown_variant[i]["alt"] > "/dev/stderr"
        }
      }
    }

    # For each position of the current sequence
    for (i = 1; i <= n; ++i) {

      # define coverage for reads starting at position i
      cov = min_coverage + int(rand() * (max_coverage - min_coverage + 1))

      # For each read starting at current position of current sequence
      for (c = 0; c < cov; ++c) {
        # define read length
        l = min_l + int(rand() * (max_l - min_l + 1))
        if (i + l > n) {
          continue
        }

        # Building current read
        read = ""
        hdr = "@" name "[" i ".." (i + l - 1) "]"

        delta = 0
        for (j = 0; j < l; ++j) {

          nucl = sseq[i + j]
          if ((i + j) in variants) {
            if (length(variants[i + j]["ref"]) == 1) {
              assert(nucl == variants[i + j]["ref"], "nucl is '" nucl "' and variants[" (i + j) "]['ref'] = '" variants[i + j]["ref"] "'");
            }
            old_j=j
            nucl = variants[i + j]["alt"]
            delta += length(nucl) - 1
            hdr = hdr " " variants[i + j]["hdr"] " (at read pos " (j + delta + 1) ")"
            j += length(variants[i + j]["ref"]) - 1
          }

          error_val = rand()
          if (error_val <= total_error_rate) {
            hdr = hdr ", error "
            if (error_val < subst_error_rate) {
              hdr = hdr "subst"
              nucl = rnd_nucl(nucl)
            } else {
              error_val -= subst_error_rate;
              if (error_val < ins_error_rate) {
                hdr = hdr "ins"
                nucl = rnd_nucl(nucl) nucl
              } else {
                hdr = hdr "del"
                nucl = ""
              }
            }
            hdr = hdr " at pos " (j + 1)
          }

          read = read nucl
        }

        ++nb_reads
        reads[nb_reads]["rank"] = rand()
        reads[nb_reads]["hdr"] = hdr
        reads[nb_reads]["seq"] = read
      }
    }
  }

  asort(reads, shuffled, "reads_compare")
  for (r in shuffled) {
    print shuffled[r]["hdr"]
    print shuffled[r]["seq"]
    print "+"
    qual = ""
    n = length(shuffled[r]["seq"])
    while (n--) {
      qual = sprintf("%s%c", qual, int(rand() * 93) + 33)
    }
    print qual
  }

}
