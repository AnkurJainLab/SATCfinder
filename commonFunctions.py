
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
        if _sequence[_base] in _degenerateBases:
            _editedSequence = _editedSequence + _degenerateBases[_sequence[_base]]
        else:
            _editedSequence = _editedSequence + _sequence[_base]
    _sequenceForward = _editedSequence
    _sequenceReverse = getReverseComplement(_sequenceForward).replace('[', '~').replace(']', '[').replace('~', ']')
    return _sequenceForward, _sequenceReverse
