import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S', level=logging.INFO)

VERSION = 0.1

def getReverseComplement(sequence):
    # First, create a list of the complementary bases
    comp = [getComplementBase(base) for base in sequence]
    # Now implode the list and reverse
    return ("".join(comp))[::-1]


def getComplement(sequence):
    # First, create a list of the complementary bases
    comp = [getComplementBase(base) for base in sequence]
    # Now implode the list
    return "".join(comp)


def getComplementBase(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'G':
        return 'C'
    elif base == 'C':
        return 'G'
    else:
        return base


def parseDegenerateSequence(_sequence):
    # Build the repeats, allowing for degenerate IUPAC bases
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

    _editedSequence = ""
    for _base in range(0, len(_sequence)):
        if _sequence[_base] in ['A', 'T', 'G', 'C', 'U']:
            _editedSequence = _editedSequence + _sequence[_base]
        elif _sequence[_base] in _degenerateBases:
            _editedSequence = _editedSequence + _degenerateBases[_sequence[_base]]
        else:
            _editedSequence = _editedSequence + _sequence[_base]

            logging.warning(f"Base {_sequence[_base]} not in nucleotides [ATGCU]"
                            f" or IUPAC degenerate bases [RYSWKMBDHVN]")
    _sequenceForward = _editedSequence
    _sequenceReverse = getReverseComplement(_sequenceForward)
    return _sequenceForward, _sequenceReverse
