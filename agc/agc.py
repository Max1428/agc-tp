#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Kermarrec Maxime"
__copyright__ = "Universite de Paris"
__credits__ = ["Kermarrec Maxime"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Kermarrec Maxime"
__email__ = "maximekermarrec14@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file, minseqlen):
    if str(amplicon_file)[len(amplicon_file)-3:len(amplicon_file)] == ".gz":
        with gzip.open(amplicon_file, "rb") as filin:
            for ligne in filin:
                if ligne.startswith(b">"):
                    pass
                else : 
                    seq = b"" + ligne.strip()
                    yield str(seq.decode('ascii'))
                #if len(str(filin.decode('UTF-8').strip())) >= minseqlen:
                    #yield filin.decode('UTF-8').strip()
    else :
        with open(amplicon_file, "r") as filin:
            for ligne in filin:
                if ligne.startswith(">"):
                    pass
                else :
                    if len(ligne.strip()) >= minseqlen:
                        yield str(ligne.strip())

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    OTU = {}
    true_OTU = {}
    for seq in read_fasta(amplicon_file,minseqlen):
        seq = str(seq)
        if seq in OTU:
            OTU[seq] = OTU[seq] + 1
        else:
            OTU[seq] = 1

    for seq in OTU:
        if OTU[seq] >= mincount:
            true_OTU[seq] = OTU[seq]
    true_OTU = sorted(true_OTU.items(), key = lambda t: t[1], reverse = True)
    return true_OTU

def read_16S(fichier):
    with open(fichier, "r") as filin:
        i = -1
        seq = []
        for ligne in filin:
            if ligne.startswith(">"):
                i  += 1
                seq.append('')
                pass
            else:
                seq[i] = str(seq[i]) + str(ligne).strip()
        return seq

def get_chunks(sequence, chunk_size):
    chunk = []
    chunk.append(sequence[0, chunk_size])
    chunk.append(sequence[chunk_size, (2*chunk_size)])
    chunk.append(sequence[(2*chunk_size), (3*chunk_size)])
    chunk.append(sequence[(3*chunk_size), (4*chunk_size)])
    return chunk

def cut_kmer(sequence, kmer_size):
    j=0
    for i in range(kmer_size,len(sequence), kmer_size):
        yield sequence[j:i]
        j = i


def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    return

#def search_mates(kmer_dict, sequence, kmer_size):
#def get_identity(alignment_list):
#def detect_chimera(perc_identity_matrix):
#def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
#def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
#def write_OTU(OTU_list, output_file):

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    OTU = dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount)
    print(OTU)
    ADN_16_S = read_16S("../data/mock_16S.fasta")
    for seq in OTU:
        k_mer = get_chunks(seq, args.chunk_size)


if __name__ == '__main__':
    main()