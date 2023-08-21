#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 23:44:55 2020

@author: Rachel Anderson
"""

import argparse
import os
import pysam
import numpy as np
import sys
import logging


def main():
    parser = argparse.ArgumentParser(description='Mate-aware trimming of repeats from reads')
    parser.add_argument('--in1', type=str, default="", required=True,
                        help='read1 file (gzipped)')
    parser.add_argument('--in2', type=str, default="", required=False,
                        help='read2 file (gzipped)')
    parser.add_argument('--out', type=str, default="", required=True,
                        help='output file')
    parser.add_argument('--sra', type=str, default="", required=True,
                        help='SRA of the study, to be embedded in the SAM file')
    parser.add_argument('--disease', type=str, default="", required=True,
                        help='Disease associated gene or CTRL, to be embedded in the SAM file')
    parser.add_argument('--minRepeats', type=int, default="4", required=True,
                        help='Minimum # of repeats to strip from reads')
    parser.add_argument('--log', type=str, default="", required=False,
                        help='File to save log information')
    parser.add_argument('--keepAllReads', action='store_true')

    args = parser.parse_args()

    find3PrimeEnds(args)


def find3PrimeEnds(_args):
    #_filePath, _whichChr, _chrStart, _chrStop, _fasta, _identifier=""

    _fileIn = pysam.AlignmentFile(_args.filePath, "rb", check_sq=False)
    _countsPerBase = np.zeros(_args.chrStop-_args.chrStart, dtype=int)

    _totalReads = 0
    for _read in _fileIn.fetch(_args.chr, _args.chrStart, _args.chrStop):
        if _args.ignoreSecondary:
            if (_read.flag & 256) != 0:
                #Skip reads which are not primary alignment
                continue
        _repeatLength = max([_read.get_tag('tL'), _read.get_tag('aL')])
        if _repeatLength <= _args.minRepeatLength:
            # Skip reads with no CAG/CTG repeats or fewer repeats than desired
            continue

        if _args.strand == "+":
            _readPositionMax = _read.get_reference_positions()[-1]
            #_readPositionMax = _readPositionMax[len(_readPositionMax)-1]
            if _readPositionMax > _args.chrStop:
                continue

            _countsPerBase[len(_countsPerBase) - (_args.chrStop-_readPositionMax)] += 1
        else:  #strand is -
            _readPositionMax = _read.get_reference_positions()[0]
            #_readPositionMax = _readPositionMax[0]
            if _readPositionMax < _args.chrStart:
                continue

            _countsPerBase[len(_countsPerBase) - (_readPositionMax-_args.chrStart)] += 1

    with open(_args.outFile, "wt") as _fileOut:
        _fileOut.write("base\t{}\n".format(_fileIn))
        for _i in range(0, len(_countsPerBase)):
            _fileOut.write("{},{}\n".format(_i + 1, _countsPerBase[_i]))

    print("Looked at", _totalReads)
    return _countsPerBase, _totalReads


if __name__ == "__main__":
    main()
