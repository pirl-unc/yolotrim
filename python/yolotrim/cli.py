# pylint: disable = no-name-in-module
from .yolotrim import trim_fastq

import argparse

parser = argparse.ArgumentParser("yolotrim", description="Heuristic adapter and poly-A trimmer for Iso-Seq reads")
parser.add_argument("--input", help="Input FASTQ file", required=True)
parser.add_argument("--output", help="Output FASTQ file", required=True)

def main():
    print("=== yolotrim ===")
    args = parser.parse_args()
    print("Args: %s" % (args, ))
    trim_fastq(args.input, args.output)
    