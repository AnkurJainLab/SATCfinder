#!/usr/bin/env python3

import logging
logging.basicConfig(format='%(asctime)s %(levelname)s : %(message)s', datefmt='%d-%m-%y %H:%M:%S', level=logging.INFO)
import argparse

import trimRepeats
import findRepeatsInGenome
import getTrimmedReadEndCounts


def main():
    # setting up to parse arguments. Nothing fancy happens here
    _parser = argparse.ArgumentParser(prog='SATCfinder', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    _subparsers = _parser.add_subparsers(dest='command', title='Choose a SATCfinder module')

    # define modules, each of which has a __command, docstring, and addArgs() function
    _modules = [findRepeatsInGenome,
                trimRepeats,
                getTrimmedReadEndCounts]
    for _module in _modules:
        # create new subparser for module
        _whichParser = _subparsers.add_parser(_module.__command, help=_module.__doc__)
        # add args from module
        _module.addArgs(_whichParser)
        # this allows the module to automatically run
        _whichParser.set_defaults(func=_module.main)

    _args = _parser.parse_args()
    # if no args passed, assume user wants help
    if(_args.command == None):
        _parser.print_help()
        return
    _args.func(_args)

if __name__ == '__main__':
    main()