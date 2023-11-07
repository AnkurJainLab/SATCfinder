#!/usr/bin/env python
'''Stand-alone module for mate-aware trimming of repeats from reads. Does not perform any other steps of pipeline
such as removing adapters or aligning to genome.'''
__command = 'trim'

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S', level=logging.INFO)
import commonFunctions as CF
import copy
import gzip
import re
import time
from itertools import zip_longest
import numpy as np


def addArgs(_parser):
    _parser.add_argument('--inFASTQ1', type=str, default="", required=True, metavar='x',
                         help='Path to input read1 fastq file (gzipped)')
    _parser.add_argument('--inFASTQ2', type=str, default="", required=False, metavar='x',
                         help='(optional) Path to input read2 fastq file (gzipped)')
    _parser.add_argument('--outSAM', type=str, default="", required=True, metavar='x',
                         help='Path to output SAM file (gzipped)')
    _parser.add_argument('--note', type=str, default="", required=False, metavar='x',
                         help='(optional) Note (e.g. dataset ID) to save in read tags. Will be saved with every read, so'
                             'should be short and only contain alphanumeric characters.')
    _parser.add_argument('--minRepeats', type=int, default=3, required=True, metavar='x',
                         help='Minimum # of repeats to strip from reads')
    _parser.add_argument('--UMIlength', type=int, default=0, required=False, metavar='x',
                         help='Length of UMI (barcode) present, assumed to be first bases in read1')
    _parser.add_argument('--minTrimmedLength', type=int, default=15, required=False, metavar='x',
                         help='Minimum length of read after trimming needed to output')
    _parser.add_argument('--repeatSequence', type=str, default="CAG", required=True, metavar='x',
                         help='Forward repeat to strip from reads. Reverse complement is also searched. '
                             'Degenerate IUPAC bases [RYSWKMBDHVN] are accepted.')
    _parser.add_argument('--outLog', type=str, default="", required=False, metavar='x',
                         help='(optional) Path to save trimming log file')
    _parser.add_argument('--keepAllReads', action='store_true', help='(flag) If present, output all reads to SAM file,'
                                                                    'even those without repeats.')

    return _parser


def trimRepeats(_args):
    logging.info(f"Starting run with options {_args}")
    _tStart = time.time()

    _repeatDistribution = np.zeros(100)
    _i = 0
    _numReads = 0
    _processedReads = {}
    _readOutcomes = {'tooShort': 0, 'noRepeat': 0, 'good': 0}
    _position = 0

    # Build the repeats, allowing for degenerate IUPAC bases
    (_repeatForward, _repeatReverse) = CF.parseDegenerateSequence(_args.repeatSequence)
    _findFwdRepeat = re.compile(r'(' + _repeatForward + '){' + str(_args.minRepeats) + ',}')
    _findRevRepeat = re.compile(r'(' + _repeatReverse + '){' + str(_args.minRepeats) + ',}')

    # Open both input files (paired end) and create the output sam file
    _fileOut = gzip.open(_args.outSAM, 'wt', compresslevel=2)
    # Write SAM header
    _fileOut.write("@HD\t"
                   "VN:1.6\n")
    _fileOut.write("@PG\t"
                   "ID:SATCfinder\t"
                   "PN:SATCfinder\t"
                   f"VN:{CF.__version}\t"
                   "CL:SATCfinder " +
                   " ".join([f"{str(key)}={str(value)}" for key, value in vars(_args).items()][:-1]) + "\n")

    _fileIn1 = gzip.open(_args.inFASTQ1, "rt")
    if _args.inFASTQ2:
        _fileIn2 = gzip.open(_args.inFASTQ2, "rt")
    else:
        # This allows us to iterate over either 1 or 2 files with the same code
        _fileIn2 = [None]
    for _lineIn1, _lineIn2 in zip_longest(_fileIn1, _fileIn2):
        # Phase 1: read FASTQ in groups of 4 lines, for both reads
        if _position == 0:
            _lines = [["", ""],
                      ["", ""],
                      ["", ""],
                      ["", ""]]
        _lines[_position][0] = _lineIn1
        _lines[_position][1] = _lineIn2

        if _position == 3:
            _position = 0
            _numReads += 1
        else:
            _position += 1
            continue

        # Phase 2: process reads
        _read1, _read2 = createFASTQreads(_lines)

        # Let's not store too much in memory
        if _numReads % 1000000 == 0:
            logging.info(f"{_numReads // 1000000}M reads trimmed in {round(time.time() - _tStart, 0)} seconds")
            writeRecordsToDisk(_fileOut, _processedReads, _args)
            del _processedReads
            _processedReads = {}

        # Sometimes open(file) fails silently with empty lines rather than stopping at end of file.
        # Check for empty line and abort if found
        if _read1['ID'] == "":
            break

        _read1, _read2, _maxRepeatLength = processReads(_args, _read1, _read2, _findFwdRepeat, _findRevRepeat)

        # Save repeat length distribution (optional statistics)
        _repeatDistribution[_maxRepeatLength] += 1

        # if we got a match, save the read(s)
        if not (_maxRepeatLength >= _args.minRepeats | _args.keepAllReads):
            _readOutcomes['noRepeat'] += 1
            # We didn't find anything. Continue!
            continue

        # Error checking (not mappable by STAR)
        if len(_read1['read']) != len(_read1['phred']):
            logging.error("Fatal error: Read 1 has mismatched sequence and phred lengths")
            logging.error(_read1)
            break
        if _lineIn2 and (len(_read2['read']) != len(_read2['phred'])):
            logging.error("Fatal error: Read 2 has mismatched sequence and phred lengths")
            logging.error(_read2)
            break

        # 0-length reads kill STAR and short reads aren't mappable
        if not _lineIn2:
            if len(_read1['read']) >= _args.minTrimmedLength:
                _processedReads[_i] = _read1
                _i += 1
                _readOutcomes['good'] += 1
            else:
                _readOutcomes['tooShort'] += 1
        elif len(_read1['read']) >= _args.minTrimmedLength and \
                len(_read2['read']) >= _args.minTrimmedLength:
            _processedReads[_i] = _read1
            _processedReads[_i + 1] = _read2
            _i += 2
            _readOutcomes['good'] += 1
        else:
            _readOutcomes['tooShort'] += 1

    # Write any remaining SAM entries
    writeRecordsToDisk(_fileOut, _processedReads, _args)
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
            _logFileOut.write(f"Run time (s):\t{_tDiff // 1}\n")
            _logFileOut.write(f"Total reads:\t{_numReads}\n")
            _logFileOut.write(f"Reads too short:\t{_readOutcomes['tooShort']}\n")
            _logFileOut.write(f"Reads without repeats:\t{_readOutcomes['noRepeat']}\n")
            _logFileOut.write(f"\n")
            _logFileOut.write(f"Repeat Length\t# reads\n")
            for _i in range(0, len(_repeatDistribution)):
                _logFileOut.write(f"{_i},{_repeatDistribution[_i]}\n")

    return 1

def createFASTQreads(_lines):
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

    _read1 = copy.deepcopy(_readBlank)
    _read2 = copy.deepcopy(_readBlank)
    # Expected format is four lines for each read and mate, roughly:
    # @SRR8393695.1 1/1
    # CTGGCGGCCGCGGGGACCAGCCGCGCTTTCAGCAGCACCACGGCCAGGCCGAGAAGCAGGGTGCAGGGGACACGCCGGCAGAGCCTCGCCATGGCCTAGAG
    # +
    # AAAFFJJJJJJJJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJFJJJJJJJJJJJJJJJJJJJJJJJJAFFJJJJJJJJJJJJJJJJJJJJJJ
    _numRead = 0
    for _whichRead in [_read1, _read2]:
        if (_numRead == 1) and (_lines[0][1] is None):
            # No read 2, skip
            continue
        #   Line 1: save the ID; skip the rest.
        _whichRead['ID'] = _lines[0][_numRead].strip().replace('@', '').replace(' ', '_')  # Handle spaces in names
        #   Line 2: save the read seq; skip the rest.
        _whichRead['read'] = _lines[1][_numRead].strip()
        #   Line 3: skip everything entirely.
        #   Line 4: save phred.
        _whichRead['phred'] = _lines[3][_numRead].strip()
        _numRead = 1

    return _read1, _read2


def processReads(_args, _read1, _read2, _findFwdRepeat, _findRevRepeat):
    # Do some string math to trim the repeats.
    _readNum = 1
    for _whichRead in [_read1, _read2]:
        if (_readNum == 2) and _whichRead['ID'] == '':
            continue

        if _args.UMIlength > 0:
            # Trim and save UMI from both read and phred
            _whichRead['UMI'] = _whichRead['read'][0:_args.UMIlength]
            _whichRead['UMIphred'] = _whichRead['phred'][0:_args.UMIlength]
            _whichRead['read'] = _whichRead['read'][_args.UMIlength:]
            _whichRead['phred'] = _whichRead['phred'][_args.UMIlength:]

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

    return _read1, _read2, max([_read1['repeatRevLength'],
                                _read1['repeatFwdLength'],
                                _read2['repeatRevLength'],
                                _read2['repeatFwdLength']])


def writeRecordsToDisk(_file, _data, _args):
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

        _file.write('\t'.join(_output) + "\n")

def main(_args):
    logging.info(f'Starting {__command} module')
    trimRepeats(_args)