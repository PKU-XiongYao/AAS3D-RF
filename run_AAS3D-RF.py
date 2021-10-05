#-*- coding:utf-8 -*-
"""
Created on Tue Jan 7 19:02:11 2020 

@author: XiongYao; zhoujb
"""

import sys, os, argparse
import PredMod

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="AAS3D-RF")
    parser.add_argument("--in", dest="AAS", required=True, help="Input amino acid subsitution file")
    parser.add_argument("--out", dest="OUTPUT", default=None, help="Output file name")
    parser.add_argument("--workPath", dest="PATH", required=True, help="Working Path (please use absolute path)")
#    parser.add_argument("--seqDir", dest="seqDir", required=True, help="Sequence directory name under workPath")
#    parser.add_argument("--pdbDir", dest="pdbDir", required=True, help="PDB directory name under workPath")
    parser.add_argument("--num", dest="NUM_THREADS", default=1, type=int, help="Number of threads (CPUs) to use, default=1")
    
    args = parser.parse_args()
    PredMod.DoPrediction(args.AAS, args.OUTPUT, args.PATH, args.NUM_THREADS)
