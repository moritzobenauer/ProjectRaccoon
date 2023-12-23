import argparse
import os

from raccoon.src.ui import start_racoon, welcome, tschau_kakao


if __name__ == "__main__":
    global SEQUENZEFILE
    global OUTPUTFILE
    global MONOMERFILE
    global EXPLICITBONDS
    global REMOVEDUPLICATES

    argParser = argparse.ArgumentParser()
    argParser.add_argument(
        "-s", "--sequence", help="Sequence File", type=str, default="seq.txt"
    )
    argParser.add_argument(
        "-o", "--output", help="Output File", type=str, default="out"
    )
    argParser.add_argument(
        "-m", "--monomers", help="Monomer File", type=str, default="monomers.dat"
    )
    argParser.add_argument(
        "-e", "--explicit", help="Explicit Bonds", type=bool, default=False
    )
    argParser.add_argument(
        "-r", "--removeduplicates", help="Remove Duplicates", type=bool, default=True
    )

    args = argParser.parse_args()

    SEQUENZEFILE = args.sequence
    OUTPUTFILE = args.outputs
    MONOMERFILE = args.monomers
    EXPLICITBONDS = args.explicit
    REMOVEDUPLICATES = args.removeduplicates

    if not os.path.isfile(SEQUENZEFILE):
        print("Sequence file does not exist!")
        exit(1)

    if not os.path.isfile(MONOMERFILE):
        print("Monomer file does not exist!")
        exit(1)

    welcome()
    start_racoon(SEQUENZEFILE, OUTPUTFILE, MONOMERFILE, EXPLICITBONDS, REMOVEDUPLICATES)
    tschau_kakao()
