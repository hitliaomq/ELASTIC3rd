# -*- coding: utf-8 -*-
# This is the script run elastic3rd in the terminal
# Simply run     elastic3rd
import argparse
from elastic3rd import __version__
import elastic3rd.elastic

def run(args):
    INPUT = args.INPUT
    elastic3rd.elastic.elastic3(INPUT=INPUT)

def run_e3rd():
    print("ELASTIC3RD version: " + __version__)
    print("Copyright \u00a9 FGMS @ HIT\n")

    parser = argparse.ArgumentParser(description='Run elastic3rd.')
    
    subparsers = parser.add_subparsers()

    #SUB-PROCESS: run
    prun = subparsers.add_parser("run", help="Run elastic3rd.")
    prun.add_argument("-i", "--input", dest="INPUT", type=str, default="INPUT",
                      help="The input file for Elastic3rd, Default: 'INPUT'\n")
    prun.set_defaults(func=run)

    args = parser.parse_args()

    try:
        a = getattr(args, "func")
    except AttributeError:
        parser.print_help()
        sys.exit(0)
    args.func(args)

if __name__ == '__main__':
    run_dfttk()