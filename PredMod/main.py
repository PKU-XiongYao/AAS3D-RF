#-*- coding: utf-8 -*-

"""
Created on Mon Jan 6 16:35:55 2020 

@author: XiongYao
"""


from preparation import checkStruSeqPosANDaas
from data import runSoftware
from features import addFeatures
from prediction import AAS3DrfPred


def DoPrediction(input_file, out_file, work_path, num_threads=1, uniprot_dir="UniProt_XML", file3d_dir="STRUCTURES"):
    d_AAS, d_disulfide = checkStruSeqPosANDaas(input_file, uniprot_dir, file3d_dir, work_path)
    runSoftware(d_AAS, file3d_dir, work_path, num_threads)
    addFeatures(d_AAS, d_disulfide)
    AAS3DrfPred(input_file, out_file, work_path)
