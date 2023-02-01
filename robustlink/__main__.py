"""
"""
import argparse
import logging

from robustlink.scf import SCF_main_repeat_subsampling
from robustlink import generate_metacells_rna
from robustlink import correlate_metacells_mc_rna
from robustlink import correlate_metacells_atac_rna

def create_parser():
    """
    """
    parser = argparse.ArgumentParser(prog='python -m robustlink')
    subparsers = parser.add_subparsers()

    # scfusion 
    scfusion = subparsers.add_parser('scfusion')
    scfusion.set_defaults(cmd='scfusion')
    SCF_main_repeat_subsampling.add_args(scfusion)

    # generat metacells rna
    metacell = subparsers.add_parser('metacell')
    metacell.set_defaults(cmd='metacell')
    generate_metacells_rna.add_args(metacell)

    # corr
    corr_mc = subparsers.add_parser('corr_mc')
    corr_mc.set_defaults(cmd='corr_mc')
    correlate_metacells_mc_rna.add_args(corr_mc)

    # corr
    corr_atac = subparsers.add_parser('corr_atac')
    corr_atac.set_defaults(cmd='corr_atac')
    correlate_metacells_atac_rna.add_args(corr_atac)

    return parser


if __name__ == "__main__":
    # 
    parser = create_parser()
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)

    if not hasattr(args, 'cmd'):
        raise ValueError("Wrong input, please run: `python -m robustlink --help`")
    elif args.cmd == 'scfusion':
        SCF_main_repeat_subsampling.main(args)
    elif args.cmd == 'metacell':
        generate_metacells_rna.main(args)
    elif args.cmd == 'corr_mc':
        correlate_metacells_mc_rna.main(args)
    elif args.cmd == 'corr_atac':
        correlate_metacells_atac_rna.main(args)