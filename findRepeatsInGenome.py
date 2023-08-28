#!/usr/bin/env python

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S',
                        level=logging.INFO)
import pandas as pd
import re
import time
import copy  # for deepcopy
from gtfparse import read_gtf

import sys
from Bio import SeqIO
import argparse
from commonFunctions import *


def main(argv):
    parser = argparse.ArgumentParser(description='Find desired repeats on genome and add them to annotation GTF file')
    parser.add_argument('--genomeFASTA', type=str, default="", required=True,
                        help='Path to genome FASTA file')
    parser.add_argument('--gtfFile', type=str, default="", required=True,
                        help='Path to genome annotation file (GTF format)')
    parser.add_argument('--outFile', type=str, default="", required=True,
                        help='File path to save output dataframe')
    parser.add_argument('--minRepeats', type=int, default="3", required=True,
                        help='Minimum # of times sequence must be repeated without interruption')
    parser.add_argument('--repeatSequence', type=str, default="CAG", required=True,
                        help='Forward repeat to find in genome. Reverse complement is also searched. '
                             'Degenerate IUPAC bases [RYSWKMBDHVN] are accepted but may significantly increase runtime')

    args = parser.parse_args()
    findRepeatsInGenome(args)


def parseGTF(_dataframe):
    _tStart = time.time()
    _feature = {'chr': '', 'name': '', 'id': '', 'type': '', 'featureStart': 0, 'featureEnd': 0,
                'strand': '', 'repeatStart': 0, 'repeatLength': 0, 'exonNumber': -1, 'transcript': ''}
    _numFeatures = 0

    _featureList = {}

    # some GTF files (c elegans) use gene_id, others (homo sapiens) have gene name
    _hasGeneName = False
    if "gene_name" in _dataframe.columns:
        _hasGeneName = True

    # Group by all items with same gene_id to get exons, start/stop codons, etc
    for _i, _GOI in _dataframe.groupby("gene_id"):
        if len(_GOI[_GOI['feature'] == 'gene']) == 0:
            # Skip genes which are missing a gene annotation somehow
            continue

        _whichGene = _GOI[_GOI['feature'] == 'gene'].iloc[0]
        # Add an entry for the gene. We will store the data in dicts temporarily
        # because repeatedly adding to a dataframe is slow.
        _f = copy.deepcopy(_feature)
        _f['chr'] = _whichGene['seqname']
        _f['id'] = _whichGene['gene_id']
        if _hasGeneName:
            _f['name'] = _whichGene['gene_name']
        _f['type'] = "gene"
        _f['featureStart'] = int(_whichGene['start'])
        _f['featureEnd'] = int(_whichGene['end'])
        _f['strand'] = _whichGene['strand']
        _featureList[_numFeatures] = copy.deepcopy(_f)
        _numFeatures += 1

        # get transcripts associated with gene
        _geneTranscripts = _GOI[_GOI['feature'] == 'transcript'].copy()
        #If there is at least one annotated transcript, take the transcript with the lowest Consensus CDS (CCDS). 
        #Otherwise, just take the longest transcript
        if(len(_geneTranscripts) > 0):
            _geneTranscripts['length'] = _geneTranscripts['end'] - _geneTranscripts['start']
            _geneTranscriptsCCDS = _geneTranscripts[~_geneTranscripts['ccds_id'].eq('')]
            if(len(_geneTranscriptsCCDS) > 0):
                _geneTranscripts = _geneTranscriptsCCDS.sort_values('ccds_id', ascending=True)
            else: 
                _geneTranscripts = _geneTranscripts.sort_values('length', ascending=False)
            _whichTranscript = _geneTranscripts.iloc[0]['transcript_id']

            # Get the exons associated with the longest transcript
            _exons = _GOI[(_GOI['feature'] == 'exon') & (_GOI['transcript_id'] == _whichTranscript)]

            # Identify start and stop codons if they exist - would be helpful to figure out UTRs later
            # Also, assume only one start/stop codon exist per transcript
            _startCodon = _GOI[
                (_GOI['feature'].str.contains('start_codon')) & (_GOI['transcript_id'] == _whichTranscript)]
            _stopCodon = _GOI[
                (_GOI['feature'].str.contains('stop_codon')) & (_GOI['transcript_id'] == _whichTranscript)]
            if len(_startCodon) > 0:
                _startCodon = _startCodon.iloc[0]
                _f['type'] = "startCodon"
                _f['featureStart'] = _startCodon['start']
                _f['featureEnd'] = _startCodon['end']
                _f['transcript'] = _whichTranscript
                _featureList[_numFeatures] = copy.deepcopy(_f)
                _numFeatures += 1
            if len(_stopCodon) > 0:
                _stopCodon = _stopCodon.iloc[0]
                _f['type'] = "stopCodon"
                _f['featureStart'] = _stopCodon['start']
                _f['featureEnd'] = _stopCodon['end']
                _f['transcript'] = _whichTranscript
                _featureList[_numFeatures] = copy.deepcopy(_f)
                _numFeatures += 1

            # Now, add the exons
            for _j, _exon in _exons.iterrows():
                _f['type'] = "exon"
                _f['repeatStart'] = -1
                _f['repeatLength'] = -1
                _f['exonNumber'] = _exon['exon_number']
                _f['featureStart'] = _exon['start']
                _f['featureEnd'] = _exon['end']
                _f['transcript'] = _whichTranscript
                _featureList[_numFeatures] = copy.deepcopy(_f)
                _numFeatures += 1
    # Store that all in a dataframe.
    _features = pd.DataFrame.from_dict(_featureList, "index")
    _featureList = {}

    # Now, go back to each gene and generate features for the introns
    for _i, _GOI in _features.groupby("id"):
        _exons = _GOI[_GOI['type'] == 'exon'].sort_values('featureStart', ascending=True)

        if len(_exons) > 0:
            _firstExonStart = _exons['featureStart'].iloc[0]
            _lastExonEnd = _exons['featureEnd'].iloc[len(_exons) - 1]
            _gene = _GOI[_GOI['type'] == 'gene'].iloc[0]

            _f = copy.deepcopy(_feature)
            _f['chr'] = _gene['chr']
            _f['id'] = _gene['id']
            _f['name'] = _gene['name']
            _f['strand'] = _gene['strand']
            _f['type'] = 'intron'
            _f['repeatStart'] = -1
            _f['repeatLength'] = -1
            _f['transcript'] = _gene['transcript']

            # Add introns
            _previousExonEnd = -1
            _intronNum = 1
            for _j, _exon in _exons.iterrows():
                if _previousExonEnd == -1:
                    _previousExonEnd = _exon['featureEnd']
                else:
                    _f['featureStart'] = _previousExonEnd + 1
                    _f['featureEnd'] = _exon['featureStart'] - 1
                    _f['exonNumber'] = _intronNum
                    _intronNum += 1
                    _featureList[_numFeatures] = copy.deepcopy(_f)
                    _numFeatures += 1
                    _previousExonEnd = _exon['featureEnd']
    _newFeatures = pd.DataFrame.from_dict(_featureList, "index")
    _features = _features.append(_newFeatures)

    _tEnd = time.time()
    _tDiff = _tEnd - _tStart
    logging.info(f"Parsed {_numFeatures} features in {round(_tDiff / 60, 1)} minutes")
    return _features


def findAllRepeats(_args):
    _tStart = time.time()
    _CHROMOSOME = {}

    for _whichChromosome in SeqIO.parse(_args.genomeFASTA, "fasta"):
        logging.info(f"Parsing {_whichChromosome.id} / {len(_whichChromosome)} bp")
        _CHROMOSOME[_whichChromosome.id] = {}
        _CHROMOSOME[_whichChromosome.id]['ID'] = _whichChromosome.id
        _CHROMOSOME[_whichChromosome.id]['sequence'] = str(_whichChromosome.seq).upper()
        _CHROMOSOME[_whichChromosome.id]['length'] = len(_whichChromosome)

    _feature = {'chr': '', 'name': '', 'id': '', 'type': '', 'featureStart': -1, 'featureEnd': -1,
                'strand': '', 'repeatStart': 0, 'repeatLength': 0, 'exonNumber': -1, 'transcript': ''}
    _numFeatures = 0
    _featureList = {}
    _repeatFwd, _repeatRev = parseDegenerateSequence(_args.repeatSequence)
    _findRepeat = {'+': re.compile(r'(?:' + _args.repeatSequence + '){' + str(_args.minRepeats) + ',}'),
                   '-': re.compile(r'(?:' + _repeatRev + '){' + str(_args.minRepeats) + ',}')}

    for _strand in ['+', '-']:
        for _whichChr in _CHROMOSOME:
            _results = _findRepeat[_strand].finditer(_CHROMOSOME[_whichChr]['sequence'])
            for _match in _results:
                _f = copy.deepcopy(_feature)
                _f['chr'] = _whichChr
                _f['type'] = "intergene"
                _f['strand'] = _strand
                _f['repeatStart'] = _match.start()
                _f['repeatLength'] = len(_match[0]) // 3

                _featureList[_numFeatures] = copy.deepcopy(_f)
                _numFeatures += 1

    _featureList = pd.DataFrame.from_dict(_featureList, "index")

    _tEnd = time.time()
    _tDiff = _tEnd - _tStart
    logging.info("Found {} repeats in {} minutes".format(_numFeatures, round(_tDiff / 60, 1)))
    return _featureList


def addRepeatsToFeatures(_geneFeatures, _allRepeats):
    _tStart = time.time()
    _toAdd = {}
    _test = 0
    _mergedRepeats = copy.deepcopy(_geneFeatures)

    for _i, _row in _allRepeats.iterrows():
        _inGene = _geneFeatures.query("(chr == '{}') & \
                                       ((featureStart <= {} <= featureEnd) | (featureStart <= {} <= featureEnd))\
                                       ".format(_row['chr'],
                                                _row['repeatStart'],
                                                _row['repeatStart'] + 3 * _row['repeatLength']))
        if len(_inGene) == 0:
            # It's intergenic. Add a new feature for it
            _toAdd[len(_toAdd)] = copy.deepcopy(_row)
        else:
            # It's in a gene region
            # is it on the right strand as the GOI? (e.g. +strand CAG with +strand gene)
            for _ignore, _GOI in _inGene.iterrows():  # account for overlapping genes
                if _GOI['strand'] == _row['strand']:
                    _row['type'] = _GOI['type']
                    _row['name'] = _GOI['name']
                    _row['id'] = str(_GOI['id'])
                    _row['featureStart'] = int(_GOI['featureStart'])
                    _row['featureEnd'] = int(_GOI['featureEnd'])
                    _row['exonNumber'] = _GOI['exonNumber']
                    _toAdd[len(_toAdd)] = copy.deepcopy(_row)

    _mergedRepeats = _mergedRepeats.append(pd.DataFrame.from_dict(_toAdd, "index"))
    _tDiff = time.time() - _tStart
    logging.info("Added repeats to features")
    return _mergedRepeats


def findRepeatsInGenome(_args):
    logging.info(f"Starting run. Will find all instances of {_args.minRepeats}x({_args.repeatSequence})")
    logging.info(f"Parsing features in gtf file")
    _gtf = read_gtf(_args.gtfFile)
    _features = parseGTF(_gtf)

    logging.info(f"Finding repeats")
    _allRepeats = findAllRepeats(_args)
    logging.info(f"Assigning repeats to features")
    _mergedFeatures = addRepeatsToFeatures(_features, _allRepeats)

    _mergedFeatures.to_csv(_args.outFile)
    logging.info(f"Done!")


if __name__ == "__main__":
    main(sys.argv[1:])
