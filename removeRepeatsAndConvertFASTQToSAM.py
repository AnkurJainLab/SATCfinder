#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 23:44:55 2020

@author: Rachel Anderson
"""
import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S', level=logging.INFO)
import commonFunctions as CF
import argparse
import copy
import gzip
import re
import time
from itertools import zip_longest
import numpy as np


def main():

    parser = argparse.ArgumentParser(description='Mate-aware trimming of repeats from reads')
    parser.add_argument('--inFASTQ1', type=str, default="", required=True,
                        help='path to read1 fastq file (gzipped)')
    parser.add_argument('--inFASTQ2', type=str, default="", required=False,
                        help='path to read2 fastq file (gzipped)')
    parser.add_argument('--outSAM', type=str, default="", required=True,
                        help='path to output SAM file (gzipped)')
    parser.add_argument('--note', type=str, default="", required=False,
                        help='Optional note (e.g. dataset ID) to save in read tags. Will be saved with every read, so'
                             'Should be short and only alphanumeric characters.')
    parser.add_argument('--minRepeats', type=int, default=3, required=True,
                        help='Minimum # of repeats to strip from reads')
    parser.add_argument('--UMIlength', type=int, default=0, required=False,
                        help='Length of UMI (barcode) present, assumed to be first bases in read1')
    parser.add_argument('--minReadLengthAfterTrim', type=int, default=15, required=False,
                        help='Minimum length of read after trimming needed to output')
    parser.add_argument('--repeatSequence', type=str, default="CAG", required=True,
                        help='Forward repeat to strip from reads. Reverse complement is also searched. '
                             'Degenerate IUPAC bases [RYSWKMBDHVN] are accepted.')
    parser.add_argument('--outLog', type=str, default="", required=False,
                        help='path to save log information file')
    parser.add_argument('--keepAllReads', action='store_true', help='If present, output all reads to SAM file,'
                                                                    'even those without repeats.')

    args = parser.parse_args()

    removeRepeatsAndConvertFASTQToSAM(args)


def removeRepeatsAndConvertFASTQToSAM(_args):
    logging.info(f"Starting run with options {_args}")
    _tStart = time.time()
    # Structure for every read:
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
    _numReads = 0
    _data = {}
    _readOutcomes = {'tooShort': 0, 'noRepeat': 0, 'good': 0}
    _position = 0
    # Build the repeats, allowing for degenerate IUPAC bases
    (_repeatForward, _repeatReverse) = CF.parseDegenerateSequence(_args.repeatSequence)
    _findFwdRepeat = re.compile(r'(' + _repeatForward + '){' + str(_args.minRepeats) + ',}')
    _findRevRepeat = re.compile(r'(' + _repeatReverse + '){' + str(_args.minRepeats) + ',}')

    # Open both input files (paired end) and create the output sam file
    # using popen because it is slightly faster than built-in gzip. Not a lot, but..
    _fileOut = gzip.open(_args.outSAM, 'wt', compresslevel=2)
    # Write SAM header
    _fileOut.write("@HD\t"
                   "VN:1.6\tSO:queryname\n")
    _fileOut.write("@PG\t"
                   "ID:SATCfinder\t"
                   "PN:SATCfinder\t"
                   f"VN:{CF.VERSION}\t"
                   "CL:removeRepeatsAndConvertFASTQToSAM " +
                   "\t".join([f"--{key}\t{value}" for key, value in vars(_args).items()]) + "\n")

    _fileIn1 = gzip.open(_args.inFASTQ1, "rt")
    if _args.inFASTQ2:
        _fileIn2 = gzip.open(_args.inFASTQ2, "rt")
    else:
        # This is a hack to allow us to iterate over either 1 or 2 files with the same code
        _fileIn2 = [None]
    for _lineIn1, _lineIn2 in zip_longest(_fileIn1, _fileIn2):
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
            _read1['ID'] = _lineIn1.strip().replace('@', '').replace(' ', '_')  # Handle spaces in names

            if _lineIn2:
                _read2 = copy.deepcopy(_readBlank)
                _read2['ID'] = _lineIn2.strip().replace('@', '').replace(' ', '_')
            _position += 1
            continue

        if _position == 1:
            _read1['read'] = _lineIn1.strip()
            if _lineIn2:
                _read2['read'] = _lineIn2.strip()
            _position += 1
            continue

        if _position == 2:
            _position += 1
            continue

        if _position == 3:
            _read1['phred'] = _lineIn1.strip()
            if _lineIn2:
                _read2['phred'] = _lineIn2.strip()
            _position = 0
            _numReads += 1

        # Let's not store too much in memory
        if _numReads % 1000000 == 0:
            logging.info(f"{_numReads // 1000000}M reads trimmed in {round(time.time() - _tStart, 0)} seconds")
            writeRecordsToDisk(_fileOut, _data, _args)
            del _data
            _data = {}
        # Sometimes open fails silently (with empty lines) rather than stopping at end of file.
        # Check for empty line and abort if found
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
        _readNum = 1
        for _whichRead in [_read1, _read2]:
            if (_readNum == 2) and not _lineIn2:
                continue
            _fwdRepeatMatches = re.search(_findFwdRepeat, _whichRead['read'])
            if _fwdRepeatMatches is not None:
                _hit = True
                _start = _fwdRepeatMatches.start(0)
                _whichRead['repeatFwdSequence'] = _whichRead['read'][_start:]
                _whichRead['repeatFwdPhred'] = _whichRead['phred'][_start:]
                _whichRead['repeatFwdLength'] = (_fwdRepeatMatches.end(0) - _start) // 3
                _whichRead['read'] = _whichRead['read'][:_start]
                _whichRead['phred'] = _whichRead['phred'][:_start]

            # Skip to last instance of repeat in read. Fix for interrupted repeats
            _revRepeatMatches = re.finditer(_findRevRepeat, _whichRead['read'])
            _result = None
            for _result in _revRepeatMatches:
                pass

            if _result is not None:
                _hit = True
                _end = _result.end(0)
                _whichRead['repeatRevSequence'] = _whichRead['read'][:_end]
                _whichRead['repeatRevPhred'] = _whichRead['phred'][:_end]
                _whichRead['repeatRevLength'] = (_end - _result.start(0)) // 3
                _whichRead['read'] = _whichRead['read'][_end:]
                _whichRead['phred'] = _whichRead['phred'][_end:]
            _readNum = 2

        # Save repeat length distribution (optional statistics)
        _repeatDistribution[max([_read1['repeatRevLength'], _read1['repeatFwdLength'], _read2['repeatRevLength'],
                                 _read2['repeatFwdLength']])] += 1

        # if we got a match, save the read(s)
        if _hit | _args.keepAllReads:
            if len(_read1['read']) != len(_read1['phred']):
                logging.error("Fatal error: Read has mismatched sequence and phred lengths", _read1)
                break
            if _lineIn2 and (len(_read2['read']) != len(_read2['phred'])):
                logging.error("Fatal error: Read has mismatched sequence and phred lengths", _read2)
                break
            # 0-length reads kill STAR. Also, short reads aren't mappable
            if not _lineIn2:
                if len(_read1['read']) >= _args.minReadLengthAfterTrim:
                    _data[_i] = _read1
                    _i += 1
                    _readOutcomes['good'] += 1
                else:
                    _readOutcomes['tooShort'] += 1
            elif len(_read1['read']) >= _args.minReadLengthAfterTrim and \
                    len(_read2['read']) >= _args.minReadLengthAfterTrim:
                _data[_i] = _read1
                _data[_i + 1] = _read2
                _i += 2
                _readOutcomes['good'] += 1
            else:
                _readOutcomes['tooShort'] += 1
        else:
            _readOutcomes['noRepeat'] += 1
            # We didn't find anything. Continue!
            continue

    # Write any remaining SAM entries
    writeRecordsToDisk(_fileOut, _data, _args)
    _fileOut.close()

    _tEnd = time.time()
    _tDiff = _tEnd - _tStart
    logging.info(f"Processed {str(_numReads)} reads in {round(_tDiff, 2)} seconds ({_numReads / _tDiff // 1}/sec)")
    logging.info(f"{_readOutcomes['good']} reads ({round(_readOutcomes['good'] / _numReads * 100,2)}%) "
                 f"had repeats and were successfully trimmed and converted to SAM format")
    logging.info(f"{_readOutcomes['tooShort']} reads ({round(_readOutcomes['tooShort'] / _numReads * 100, 2)}%)"
                 f" had repeats, but were discarded because the non-repetitive region "
                 f"was shorter than --minReadLength")
    logging.info(f"{_readOutcomes['noRepeat']} reads ({round(_readOutcomes['noRepeat'] / _numReads * 100, 2)}%)"
                 f" had no repeats matching >={_args.minRepeats}x{_repeatForward}/{_repeatReverse}")

    if _args.outLog != "":
        with open(_args.outLog, "wt") as _logFileOut:
            _logFileOut.write(f"Run time (s):,{_tDiff // 1}\n")
            _logFileOut.write(f"Total reads:,{_numReads}\n")
            _logFileOut.write(f"Reads too short,{_readOutcomes['tooShort']}\n")
            _logFileOut.write(f"Reads without repeats,{_readOutcomes['noRepeat']}\n")
            for _i in range(0, len(_repeatDistribution)):
                _logFileOut.write(f"{_i},{_repeatDistribution[_i]}\n")


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
                   'rS:Z:' + _args.note,
                   'RX:Z:' + _line['UMI'],
                   'pX:Z:' + _line['UMIphred']]

        _fileOut.write('\t'.join(_output) + "\n")


if __name__ == "__main__":
    main()
