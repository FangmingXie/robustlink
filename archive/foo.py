
import sys

def main():
    config_py = "config_mop_rna_mc_ka30_knn30_211115"
    sys.path.insert(0, '/Users/fangmingxie/Documents/Code/robustlink/robustlink/scf/configs')
    print(f"from {config_py} import *")
    exec(f"from {config_py} import *")
    print(outdir)

main()