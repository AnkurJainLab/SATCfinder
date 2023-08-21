#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 23:44:55 2020

@author: Rachel Anderson
"""
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S', level=logging.INFO)

import argparse
import copy
import gzip
import re
import sys
import time
from itertools import zip_longest
import numpy as np
import rachelcore as rc



def main(argv):
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
    parser.add_argument('--UMIlength', type=int, default=0, required=False,
                        help='Length of UMI (barcode) present, assumed to be first bases in read1')
    parser.add_argument('--repeatSequence', type=str, default="CAG", required=True,
                        help='Forward repeat to strip from reads. Reverse complement is also searched. '
                             'Degenerate IUPAC bases [RYSWKMBDHVN] are accepted.')
    parser.add_argument('--log', type=str, default="", required=False,
                        help='File to save log information')
    parser.add_argument('--keepAllReads', action='store_true')

    args = parser.parse_args()

    removeRepeatsAndConvertFASTQToSAM(args)


# args.sra, args.disease, args.in1, args.in2, args.out, args.minRepeats, args.UMIlength,
# args.repeatSequence, args.log, args.keepAllReads)


def removeRepeatsAndConvertFASTQToSAM(_args):
    # _SRA, _disease, _filePathIn1, _filePathIn2, _filePathOut, _minRepeats, _UMIlength,
    # _args.repeatForward, _logPath, _keepAllReads):
    _tStart = time.time()

    _readBlank = {'read': "",
                  'phred': "",
                  'ID': "",
                  'repeatFwdSequence': "*",
                  'repeatFwdLength': 0,
                  'repeatFwdPhred': "*",
                  'repeatRevSequence': "*",
                  'repeatRevLength': 0,
                  'repeatRevPhred': "*",
                  'UMI': "*",
                  'UMIphred': "*"}
    _read2 = copy.deepcopy(_readBlank)  # In case read2 is not provided, initialize as blank for calculations later
    _repeatDistribution = np.zeros(100)
    _i = 0
    _n = 1
    _data = {}
    _short = 0
    _long = 0
    # Build the repeats, allowing for degenerate IUPAC bases ## TODO replace with commonFunctions
    _degenerateBases = {'R': '[AG]',
                        'Y': '[CT]',
                        'S': '[GC]',
                        'W': '[AT]',
                        'K': '[GT]',
                        'M': '[AC]',
                        'B': '[CGT]',
                        'D': '[AGT]',
                        'H': '[ACT]',
                        'V': '[ACG]',
                        'N': '[ATGC]'}

    _editedRepeat = ""
    for _base in range(0, len(_args.repeatForward)):
        if _args.repeatForward[_base] in _degenerateBases:
            _editedRepeat = _editedRepeat + _degenerateBases[_args.repeatForward[_base]]
        else:
            _editedRepeat = _editedRepeat + _args.repeatForward[_base]
    _repeatForward = _editedRepeat
    _repeatReverse = rc.getReverseComplement(_repeatForward).replace('[', '~').replace(']', '[').replace('~', ']')
    _findFwdRepeat = re.compile(r'(' + _repeatForward + '){' + str(_args.minRepeats) + ',}')
    _findRevRepeat = re.compile(r'(' + _repeatReverse + '){' + str(_args.minRepeats) + ',}')

    _position = 0

    # Open both input files (paired end) and create the output sam file
    # using popen because it is slightly faster than built-in gzip. Not a lot, but..
    _fileOut = gzip.open(_args.filePathOut, 'wt', compresslevel=2)
    _fileOut.write("@HD\tVN:1.6\tSO:queryname\n")  # Write SAM header
    _fileIn1 = gzip.open(_args.filePathIn1)
    if _args.filePathIn2:
        _fileIn2 = gzip.open(_args.filePathIn2)
    else:
        # This is a hack to allow us to iterate over either 1 or 2 files with the same code
        _fileIn2 = [None]
    for _in1, _in2 in zip_longest(_fileIn1, _fileIn2):
        # Expected format is roughly:
        # @SRR8393695.1 1/1
        # CTGGCGGCCGCGGGGACCAGCCGCGCTTTCAGCAGCACCACGGCCAGGCCGAGAAGCAGGGTGCAGGGGACACGCCGGCAGAGCCTCGCCATGGCCTAGAG
        # +
        # AAAFFJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJAFFJJJJJJJJJJJJJJJJJJJJJJ

        # Logic here:
        #   Line 1: save the ID; skip the rest.
        #   Line 2: save the read seq; skip the rest.
        #   Line 3: skip everything entirely.
        #   Line 4: save phred.
        #   Process all four lines into SAM format.
        # Repeat.
        if _position == 0:
            _read1 = copy.deepcopy(_readBlank)
            _read1['ID'] = _in1.strip().replace('@', '').replace(' ', '_')

            if _in2:
                _read2 = copy.deepcopy(_readBlank)
                _read2['ID'] = _in2.strip().replace('@', '').replace(' ', '_')
            _position += 1
            continue

        if _position == 1:
            _read1['read'] = _in1.strip()
            if _in2:
                _read2['read'] = _in2.strip()
            _position += 1
            continue

        if _position == 2:
            _position += 1
            continue

        if _position == 3:
            _read1['phred'] = _in1.strip()
            if _in2:
                _read2['phred'] = _in2.strip()
            _position = 0
            _n += 1

        if _n % 1000000 == 0:
            logging.info(f"{_n // 1000000}M reads trimmed in {round(time.time() - _tStart, 2)} seconds")
            writeRecordsToDisk(_fileOut, _data, _args)
            del _data
            _data = {}
        # Sometimes it fails silently (with empty lines) rather than stopping at end of file.
        # Not sure why. Hack to avoid this
        if _read1['ID'] == "":
            break

        if _args.UMIlength > 0:
            # Trim and save UMI from both read and phred
            _read1['UMI'] = _read1['read'][0:_args.UMIlength]
            _read1['UMIphred'] = _read1['phred'][0:_args.UMIlength]
            _read1['read'] = _read1['read'][_args.UMIlength:]
            _read1['phred'] = _read1['phred'][_args.UMIlength:]

        # Assume we throw out reads, only keep if we find a repeat of sufficient length
        _hit = False

        _fwdRepeatFile1 = re.search(_findFwdRepeat, _read1['read'])
        if _fwdRepeatFile1 is not None:
            _hit = True
            _start = _fwdRepeatFile1.start(0)
            _read1['repeatFwdSequence'] = _read1['read'][_start:]
            _read1['repeatFwdPhred'] = _read1['phred'][_start:]
            _read1['repeatFwdLength'] = (_fwdRepeatFile1.end(0) - _start) // 3
            _read1['read'] = _read1['read'][:_start]
            _read1['phred'] = _read1['phred'][:_start]

        _revRepeatFile1 = re.finditer(_findRevRepeat, _read1['read'])
        _result = None
        for _result in _revRepeatFile1:
            pass

        if _result is not None:
            _hit = True
            _end = _result.end(0)
            _read1['repeatRevSequence'] = _read1['read'][:_end]
            _read1['repeatRevPhred'] = _read1['phred'][:_end]
            _read1['repeatRevLength'] = (_end - _result.start(0)) // 3
            _read1['read'] = _read1['read'][_end:]
            _read1['phred'] = _read1['phred'][_end:]

        if _in2:
            _fwdRepeatFile2 = re.search(_findFwdRepeat, _read2['read'])
            if _fwdRepeatFile2 is not None:
                _hit = True
                _start = _fwdRepeatFile2.start(0)
                _read2['repeatFwdSequence'] = _read2['read'][_start:]
                _read2['repeatFwdPhred'] = _read2['phred'][_start:]
                _read2['repeatFwdLength'] = (_fwdRepeatFile2.end(0) - _start) // 3
                _read2['read'] = _read2['read'][:_start]
                _read2['phred'] = _read2['phred'][:_start]

            _revRepeatFile2 = re.finditer(_findRevRepeat, _read2['read'])
            _result = None
            for _result in _revRepeatFile2:
                pass

            if _result is not None:
                _hit = True
                _end = _result.end(0)
                _read2['repeatRevSequence'] = _read2['read'][:_end]
                _read2['repeatRevPhred'] = _read2['phred'][:_end]
                _read2['repeatRevLength'] = (_end - _result.start(0)) // 3
                _read2['read'] = _read2['read'][_end:]
                _read2['phred'] = _read2['phred'][_end:]

        # Save repeat length distribution (optional statistics)
        _repeatDistribution[max([_read1['repeatRevLength'], _read1['repeatFwdLength'], _read2['repeatRevLength'],
                                 _read2['repeatFwdLength']])] += 1

        # if we got a match, save the read(s)
        if _hit | _args.keepAllReads:
            if len(_read1['read']) != len(_read1['phred']):
                rc.msg("Fatal error: Mismatched lengths!", _read1)
                break
            if _in2 and (len(_read2['read']) != len(_read2['phred'])):
                rc.msg("Fatal error: Mismatched lengths!", _read2)
                break
            # 0-length reads kill STAR. Also, short reads aren't mappable... so let's omit anything shorter than 8 bases
            if not _in2:
                if len(_read1['read']) >= 8:
                    _data[_i] = _read1
                    _i += 1
                    _long += 1
                else:
                    _short += 1
            elif len(_read1['read']) >= 8 and len(_read2['read']) >= 8:
                _data[_i] = _read1
                _data[_i + 1] = _read2
                _i += 2
                _long += 1
            else:
                _short += 1
        else:
            # We didn't find anything. Continue!
            continue

    # Write any remaining SAM entries
    writeRecordsToDisk(_fileOut, _data, _args)
    _fileOut.close()

    _tEnd = time.time()
    _tDiff = _tEnd - _tStart
    logging.info(f"Processed {str(_n)} reads in {str(round(_tDiff, 2))} seconds ({_n / _tDiff // 1}/sec)")
    logging.info(f"{_long} reads ({_long / _n * 100}%) had repeats and were successfully trimmed "
                 f"and converted to SAM format")
    logging.info(f"{_short} reads ({_short / _n * 100}%) had repeats, were discarded because the non-repetitive region "
                 f"was shorter than --minReadLength")
    logging.info(_repeatDistribution)
    if _args.logPath != "":
        with open(_args.logPath, "wt") as _logFileOut:
            for _i in range(0, len(_repeatDistribution)):
                _logFileOut.write("{},{}\n".format(_i, _repeatDistribution[_i]))
            _logFileOut.write("Run time (s):,{}\n".format(_tDiff // 1))
            _logFileOut.write("Total reads:,{}\n".format(_n))
            _logFileOut.write("Reads too short,{}\n".format(_short))


def writeRecordsToDisk(_fileOut, _data, _args):
    for _line in _data.values():
        _output = [_line['ID'],
                   "0",
                   '*',
                   '0',
                   '0',
                   '*',
                   '*',
                   '0',
                   '0',
                   _line['read'],
                   _line['phred'],
                   'aS:Z:' + _line['repeatFwdSequence'],
                   'aL:i:' + str(_line['repeatFwdLength']),
                   'aP:Z:' + _line['repeatFwdPhred'],
                   'tS:Z:' + _line['repeatRevSequence'],
                   'tL:i:' + str(_line['repeatRevLength']),
                   'tP:Z:' + _line['repeatRevPhred'],
                   'rS:Z:' + _args.SRA,
                   'rD:Z:' + _args.disease,
                   'RX:Z:' + _line['UMI'],
                   'pX:Z:' + _line['UMIphred']]

        _fileOut.write('\t'.join(_output) + "\n")


if __name__ == "__main__":
    main(sys.argv[1:])
