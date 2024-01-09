import argparse
import os

from project_raccoon import start_racoon, welcome, tschau_kakao


def main():
    argParser = argparse.ArgumentParser()
    argParser.add_argument(
        "-s", "--sequence", help="Sequence File", type=str, default="seq.txt"
    )
    argParser.add_argument(
        "-o", "--output", help="Output File", type=str, default="out.pdb"
    )
    argParser.add_argument(
        "-m", "--monomers", help="Monomer File", type=str, default=None
    )
    argParser.add_argument(
        "-e", "--explicit", help="Explicit Bonds", type=bool, default=False
    )
    argParser.add_argument(
        "-r", "--removeduplicates", help="Remove Duplicates", type=bool, default=True
    )

    args = argParser.parse_args()

    SEQUENCEFILE = args.sequence
    OUTPUTFILE = args.output
    MONOMERFILE = args.monomers
    EXPLICITBONDS = args.explicit
    REMOVEDUPLICATES = args.removeduplicates

    if not os.path.isfile(SEQUENCEFILE):
        print("Sequence file does not exist!")
        exit(1)

    if MONOMERFILE is not None and not os.path.isfile(MONOMERFILE):
        print("Monomer file does not exist!")
        exit(1)

    welcome()
    start_racoon(
        sequence_file=SEQUENCEFILE,
        out_file=OUTPUTFILE,
        monomer_file=MONOMERFILE,
        explicitbonds=EXPLICITBONDS,
        remove_duplicates=REMOVEDUPLICATES,
    )
    tschau_kakao()


if __name__ == "__main__":
    main()
