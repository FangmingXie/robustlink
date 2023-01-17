"""
"""
import argparse
import logging

from scf import SCF_main_repeat_subsampling
import generate_metacells_rna
import correlate_metacells_mc_rna

def create_parser():
    """
    """
    parser = argparse.ArgumentParser(prog='robustlink')
    subparsers = parser.add_subparsers()
    #  choices=['scf', 'metacell', 'corr']

    # scf 
    scf = subparsers.add_parser('scf')
    scf.set_defaults(cmd='scf')
    SCF_main_repeat_subsampling.add_args(scf)

    # generat metacells rna
    metacell = subparsers.add_parser('metacell')
    metacell.set_defaults(cmd='metacell')
    generate_metacells_rna.add_args(metacell)

    # corr
    corr = subparsers.add_parser('corr')
    corr.set_defaults(cmd='corr')
    correlate_metacells_mc_rna.add_args(corr)

    return parser


if __name__ == "__main__":
    # 
    parser = create_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if args.cmd == 'scf':
        SCF_main_repeat_subsampling.main(args)

    if args.cmd == 'metacell':
        generate_metacells_rna.main(args)

    if args.cmd == 'corr':
        correlate_metacells_mc_rna.main(args)