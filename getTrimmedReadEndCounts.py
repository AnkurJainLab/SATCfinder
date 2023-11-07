#!/usr/bin/env python
'''Find trimmed ends in aligned BAM file'''
__command = 'ends'

import argparse
import pysam
import numpy as np
import re
import logging


def addArgs(_parser):

    _parser.add_argument('--inBAM', type=str, default="", required=True, metavar='x',
                         help='Input BAM file to search for trimmed ends. The BAM file should have been processed by '
                              'SATCfinder (e.g. it should have SAM attributes tL/aL to indicate the number of repeats'
                              'trimmed), and must be indexed.')
    _parser.add_argument('--outTSV', type=str, default="", required=True, metavar='x',
                         help='Output TSV file.')
    _parser.add_argument('--region', type=str, default="", required=True, metavar='x',
                         help='Chromosome region to process, for example: chr4:3073876-3075876')
    _parser.add_argument('--strand', type=str, default="", required=True, metavar='x',
                     help='Strand of target region')
    _parser.add_argument('--minRepeats', type=int, default=3, required=False, metavar='x',
                         help='Minimum # of repeats required to output. Default 3')
    _parser.add_argument('--ignoreSecondary', action='store_true',
                         help='If present, ignore secondary alignments of multimapping reads, but include primary '
                              'alignments.')
    _parser.add_argument('--ignoreMultimapping', action='store_true', help='If present, ignore all multimapping reads.')
    return _parser


def find3PrimeEnds(_args):
    _processRegion = re.match('^(.+):(.+)-(.+)', _args.region)
    _chromosome, _start, _end = _processRegion.groups()
    _start = int(_start)
    _end = int(_end)
    if _end == 0 | _start >= _end:
        logging.error("Invalid start/end coordinates")
        return
    _countsPerBase = np.zeros(_end - _start, dtype=int)
    _totalReads = 0
    with pysam.AlignmentFile(_args.inBAM, "rb", check_sq=False) as _fileIn:
        for _read in _fileIn.fetch(_chromosome, _start, _end):
            _totalReads += 1

            if _args.ignoreSecondary:
                if (_read.flag & 256) != 0:
                    # Skip reads which are not primary alignment
                    continue
            if _args.ignoreMultimapping:
                if int(_read.get_tag('NH')) > 1:
                    # Skip reads with multiple alignments
                    continue

            _repeatLength = max([_read.get_tag('tL'), _read.get_tag('aL')])
            if _repeatLength < _args.minRepeats:
                # Skip reads with no repeats or fewer repeats than desired
                continue

            if _args.strand == "+":
                _readPositionMax = _read.get_reference_positions()[-1]
                if _readPositionMax > _end:
                    continue
                _countsPerBase[len(_countsPerBase) - (_end - _readPositionMax)] += 1
            else:  # strand is -
                _readPositionMax = _read.get_reference_positions()[0]
                if _readPositionMax < _start:
                    continue
                _countsPerBase[len(_countsPerBase) - (_readPositionMax - _start)] += 1

    with open(_args.outTSV, "wt") as _fileOut:
        _fileOut.write(f"base\t{_args.inBAM}_{_args.region}\n")
        for _i in range(0, len(_countsPerBase)):
            _fileOut.write("{}\t{}\n".format(_i + 1 + _start, _countsPerBase[_i]))

    logging.info(f"Done! Looked at {_totalReads} reads.")


def main(_args):
    logging.info(f'Starting {__command} module')
    find3PrimeEnds(_args)