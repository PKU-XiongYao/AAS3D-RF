# -*- coding: utf-8 -*-
import os

"""
Created on Fri Jan 3 12:11:26 2020 

@author: XiongYao
"""

##*************************************##
## Set paths for software and database ##
##*************************************##



## Paths not to change !!!

# Path for AAS3D-RF predictor
AAS3Drf = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Path for PSI-BLAST 
PSIBLAST = os.path.join(AAS3Drf, "Tools/blast-2.2.29+/bin")

# Path for mkdssp in DSSP
DSSP = os.path.join(AAS3Drf, "Tools/DSSP-v2.2.1")



## Paths provided by user

# Set database for PSI-BLAST (absolute path)
BLAST_DATABASE = os.path.join(AAS3Drf, "Tools/blast-2.2.29+/UniRef90/uniref90")

# Set RING path (where is "bin" folder; absolute path)
RING = os.path.join(AAS3Drf, "Tools/RING")  

# Set Naccess path (where is "naccess" executable file; absolute path)
Naccess = os.path.join(AAS3Drf, "Tools/Naccess")

# Set FoldX path (where is "foldx" executable file; absolute path)
FoldX = os.path.join(AAS3Drf, "Tools/FoldX")



'''
## Paths for author: XiongYao
AAS3Drf = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PSIBLAST = os.path.join(AAS3Drf, "Tools/blast-2.2.29+/bin")
DSSP = os.path.join(AAS3Drf, "Tools/DSSP-v2.2.1")
BLAST_DATABASE = "/home/xiongy/database/uniprot201806/uniref90.fasta"
RING = "/home/xiongy/software/Ring"
Naccess = "/home/xiongy/software/naccess2.1.1"
FoldX = "/home/xiongy/software/FoldX-v.3"
'''
