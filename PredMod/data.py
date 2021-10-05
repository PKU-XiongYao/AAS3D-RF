#-*- coding:utf-8 -*-

"""
Created on Sat Jan 4 08:25:00 2020 

@author: XiongYao
"""

import os, sys, logging, subprocess, multiprocessing
from functools import partial
from userConfig import AAS3Drf, PSIBLAST, BLAST_DATABASE, DSSP, RING, Naccess, FoldX


##*****************##
## Build directory ##
##*****************##

for dir_name in ["UniProtXML", "UniProtPSSM", "MSAforRING", "RING", "DSSP", "Naccess", "FoldX/Repaired_PDB", "FoldX/structure_mut", "FoldX/BuildModel_OUT"]:
    build_dir = os.path.join(AAS3Drf, "Data", dir_name)
    try:
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)
    except OSError as e:
        pass


##*****************##        
## Run other tools ##
##*****************##        

def runPSIBLAST(uniprot_ac, num_threads, outscreen, errorscreen):
    psiblast_exe = os.path.join(PSIBLAST, "psiblast")
    query = os.path.join(AAS3Drf, "Data", "UniProtSeq", uniprot_ac+".fasta")
    xml_out = os.path.join(AAS3Drf, "Data", "UniProtXML", uniprot_ac+".xml")
    pssm_out = os.path.join(AAS3Drf, "Data", "UniProtPSSM", uniprot_ac+".pssm")
    comd = psiblast_exe+" -query "+query+" -db "+BLAST_DATABASE+" -evalue 0.0001 -num_iterations 3 -outfmt 5 \
           -num_threads "+str(num_threads)+" -out "+xml_out+" -out_ascii_pssm "+pssm_out
    try:
        retcode = subprocess.check_call(comd, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
    except subprocess.CalledProcessError as exc:
        print "***Error: UniProt full sequence PSI-BLAST***"
        print exc.cmd



def preMSAforRING(stru_id, num_threads, outscreen, errorscreen):
    psiblast_exe = os.path.join(PSIBLAST, "psiblast")
    query = os.path.join(AAS3Drf, "Data", "StructureSeq", stru_id+".fasta")
    out_msa = os.path.join(AAS3Drf, "Data", "MSAforRING", stru_id+".msa")
    comd = psiblast_exe+" -query "+query+" -db "+BLAST_DATABASE+" -evalue 0.0001 -num_iterations 3 -outfmt 4 \
            -num_threads "+str(num_threads)+" -out "+out_msa

    try:
        retcode = subprocess.check_call(comd, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
    except subprocess.CalledProcessError as exc:
        print "***Error: structure sequence region PSI-BLAST for RING***"
        print exc.cmd


def runRING(stru_id, file3d_dir, work_path, outscreen, errorscreen):
    ring_exe = os.path.join(RING, "bin/Ring")
    structure_file = os.path.join(work_path, file3d_dir, stru_id+".pdb")
    out_msa = os.path.join(AAS3Drf, "Data", "MSAforRING", stru_id+".msa")
    out_node = os.path.join(AAS3Drf, "Data", "RING", stru_id+".node")
    out_edge = os.path.join(AAS3Drf, "Data", "RING", stru_id+".edge")
    comd = ring_exe+" -i "+structure_file+" -p "+out_msa+" -n ca --first_edge --get_iac -t 6.5 -g 1 -N "+out_node+" -E "+out_edge
    try:
        retcode = subprocess.check_call(comd, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
    except subprocess.CalledProcessError as exc:
        print "***Error: construct residue-residue interaction network by RING***"
        print exc.cmd


def runDSSP(stru_id, file3d_dir, work_path, outscreen, errorscreen):
    dssp_exe = os.path.join(DSSP, "mkdssp")
    structure_file = os.path.join(work_path, file3d_dir, stru_id+".pdb")
    out_dssp = os.path.join(AAS3Drf, "Data", "DSSP", stru_id+".dssp")
    comd = dssp_exe+" -i "+structure_file+" -o "+out_dssp
    try:
        retcode = subprocess.check_call(comd, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
    except subprocess.CalledProcessError as exc:
        print "***Error: run DSSP***"
        print exc.cmd


def runNaccess(stru_id, file3d_dir, work_path, outscreen, errorscreen):
    naccess_exe = os.path.join(Naccess, "naccess")
    structure_file = os.path.join(work_path, file3d_dir, stru_id+".pdb")
    comd1 = naccess_exe+" "+structure_file
    
    out_path = os.path.join(AAS3Drf, "Data", "Naccess")
    comd2 = "mv "+os.path.join("./", stru_id+".[lra][os][ga]")+" "+out_path
    
    try:
        retcode1 = subprocess.check_call(comd1, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
        retcode2 = subprocess.check_call(comd2, shell=True)
    except subprocess.CalledProcessError as exc:
        print "***Error: run Naccess***"
        print exc.cmd



def runFoldXrepair(stru_id, file3d_dir, work_path, outscreen, errorscreen):
    foldx_exe = os.path.join(FoldX, "foldx")
    filename = stru_id+".pdb"
    old_path = os.path.join(work_path, file3d_dir)
    new_path = os.path.join(AAS3Drf, "Data", "FoldX", "Repaired_PDB")
    comd1 = foldx_exe+" --command=RepairPDB --pdb="+filename+" --pdb-dir="+old_path+" --output-dir="+new_path
    comd2 = "mv "+os.path.join(new_path,stru_id+"_Repair.pdb")+" "+os.path.join(new_path,filename)
    try:
        retcode1 = subprocess.check_call(comd1, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
        retcode2 = subprocess.check_call(comd2, shell=True)
    except subprocess.CalledProcessError as exc:
        print "***Error: FoldX repair PDB***"
        print exc.cmd



def preFoldXmut(stru_id, list_aas):
    out_filename = "individual_list_"+stru_id+".mut"
    out_mut = os.path.join(AAS3Drf, "Data", "FoldX", "structure_mut", out_filename)
    with open(out_mut, "w") as outfile:
        l_mut = list(set(list_aas))
        print >> outfile, "\n".join(l_mut)


def runFoldXBuildModel(stru_id, outscreen, errorscreen):
    foldx_exe = os.path.join(FoldX, "foldx")
    pdb_name = stru_id+".pdb"
    pdb_path = os.path.join(AAS3Drf, "Data", "FoldX", "Repaired_PDB")
    mut_name = "individual_list_"+stru_id+".mut"
    mut_file = os.path.join(AAS3Drf, "Data", "FoldX", "structure_mut", mut_name)
    out_path = os.path.join(AAS3Drf, "Data", "FoldX", "BuildModel_OUT")

    comd1 = "cp "+mut_file+" "+"./"
    comd2 = foldx_exe+" --command=BuildModel --pdb="+pdb_name+" --mutant-file="+mut_name+" --pdb-dir="+pdb_path+" --output-dir="+out_path
    comd3 = "rm "+mut_name
    try:
        retcode1 = subprocess.check_call(comd1, shell=True)
        retcode2 = subprocess.check_call(comd2, shell=True, stdout=open(outscreen,"a"), stderr=open(errorscreen,"a"))
        retcode3 = subprocess.check_call(comd3, shell=True)
    except subprocess.CalledProcessError as exc:
        print "***Error: FoldX BuildModel***"
        print exc.cmd

    

def creatParallelJobs(input_func, input_list, file3d_dir, work_path, num_threads, outscreen, errorscreen):
    func_partial = partial(input_func, file3d_dir=file3d_dir, work_path=work_path, outscreen=outscreen, errorscreen=errorscreen)
    p = multiprocessing.Pool(processes=num_threads)
    p.map(func_partial, input_list)
    p.close()
    p.join()


def ParallelFoldXBuildModel(input_func, input_list, num_threads, outscreen, errorscreen):
    f_partial = partial(input_func, outscreen=outscreen, errorscreen=errorscreen )
    p = multiprocessing.Pool(processes=num_threads)
    p.map(f_partial, input_list)
    p.close()
    p.join()



def runSoftware(d_AAS, file3d_dir, work_path, num_threads):
    l_uniprot = []
    l_stru = []
    for key in d_AAS:
        uniprot_ac, stru_id, gene_name, chain = key
        if uniprot_ac not in l_uniprot:
            l_uniprot.append(uniprot_ac)
        
        if stru_id not in l_stru:
            l_stru.append(stru_id)
        
        l_aas = []
        if chain == " ":
            for aas in d_AAS[key]:
                new_aas = aas[2:]
                l_aas.append(new_aas+";")
        else:
            for aas in d_AAS[key]:
                wt = aas[2]
                new_aas = wt+chain+aas[3:]
                l_aas.append(new_aas+";")

        preFoldXmut(stru_id, l_aas)

    OutFile = os.path.join(AAS3Drf, "tmp", "Out_Screen")
    ErrorFile = os.path.join(AAS3Drf, "tmp", "Error_Screen")
    ofile_out = open(OutFile, "w")
    ofile_out.close()
    ofile_error = open(ErrorFile, "w")
    ofile_error.close()

    logging.info("PSI-BLAST running...")
    for uniprot in l_uniprot:
        runPSIBLAST(uniprot, num_threads, OutFile, ErrorFile)

    logging.info("RING running...")
    for stru in l_stru:
        preMSAforRING(stru, num_threads, OutFile, ErrorFile)
    
    creatParallelJobs(runRING, l_stru, file3d_dir, work_path, num_threads, OutFile, ErrorFile)

    logging.info("DSSP running...")
    for stru in l_stru:
        runDSSP(stru, file3d_dir, work_path, OutFile, ErrorFile)

    logging.info("Naccess running...")
    for stru in l_stru:
        runNaccess(stru, file3d_dir, work_path, OutFile, ErrorFile)

    logging.info("FoldX running...")
    foldx_param_file = os.path.join(FoldX, "rotabase.txt")
    com_cp = "cp "+foldx_param_file+" ./"
    subprocess.check_call(com_cp, shell=True)

    creatParallelJobs(runFoldXrepair, l_stru, file3d_dir, work_path, num_threads, OutFile, ErrorFile)
    ParallelFoldXBuildModel(runFoldXBuildModel, l_stru, num_threads, OutFile, ErrorFile)



##*********************************##
## Aquire feature value from file  ##
##*********************************##

## Composition in percent for the complete SwissProt database
## (only Organism: Human)
_AA_Freq = {
    "A": 7.015, "C": 2.299, "E": 7.094, "D": 4.737, "G": 6.577, "F": 3.655,
    "I": 4.340, "H": 2.631, "K": 5.723, "M": 2.131, "L": 9.962, "N": 3.590,
    "Q": 4.764, "P": 6.306, "S": 8.318, "R": 5.645, "T": 5.356, "W": 1.220,
    "V": 5.973, "Y": 2.666}


## Residual Variation Intolerance Score (RVIS)
## Source: http://genic-intolerance.org
rvis_file = os.path.join(AAS3Drf, "PredMod", "preData", "RVIS_Unpublished_ExACv2_March2017.txt")
_rvis = {}
with open(rvis_file, "r") as ifile_rvis:
    title_rvis = ifile_rvis.readline()
    for line in ifile_rvis:
        line = line.rstrip("\n")
        l_info = line.split("\t")
        try:
            _rvis[l_info[0]] = float(l_info[3])
        except ValueError:
            pass

## Miyata matrice
## Ref: Miyata, T., Miyazawa, S., & Yasunaga, T. (1979). Two types of amino acid substitutions in protein evolution. J Mol Evol, 12(3), 219-236. 
miyata_file = os.path.join(AAS3Drf, "PredMod", "preData", "Miyata_AminoAcidPairDistance.txt")
miyata_title = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
_miyata = {}
with open(miyata_file, "r") as ifile_miyata:
    n = 0
    for line in ifile_miyata:
        line = line.rstrip("\n")
        l_info = line.split()
        for index in range(len(l_info)):
            _miyata[(miyata_title[index], miyata_title[n])] = l_info[index]
            _miyata[(miyata_title[n], miyata_title[index])] = l_info[index]
        n += 1


## Physicochemical properties 
## Ref for F1: Atchley, W. R., Zhao, J. P., Fernandes, A. D., & Druke, T. (2005). Solving the protein sequence metric problem. Proceedings of the National Academy of Sciences of the United States of America, 102(18), 6395-6400. doi:10.1073/pnas.0408677102
## Ref for Grantham volume: Grantham, R. (1974). Amino-Acid Difference Formula to Help Explain Protein Evolution. Science, 185(4154), 862-864. doi:DOI 10.1126/science.185.4154.862
## Ref for Kyte-Doolittle hydropathy index: Kyte, J., & Doolittle, R. F. (1982). A Simple Method for Displaying the Hydropathic Character of a Protein. Journal of Molecular Biology, 157(1), 105-132. doi:Doi 10.1016/0022-2836(82)90515-0
## Ref for Isoelectric point: Zimmerman, J. M., Eliezer, N., & Simha, R. (1968). The characterization of amino acid sequences in proteins by statistical methods. J Theor Biol, 21(2), 170-201. 
prop_file = os.path.join(AAS3Drf, "PredMod", "preData", "Physicochemical_Property.tab")
_F1 = {}
_Vol = {}
_Hydro = {}
_pI = {}
with open(prop_file, "r") as ifile_prop:
    title_prop = ifile_prop.readline()
    for line in ifile_prop:
        line = line.rstrip("\n")
        l_AA_value = line.split("\t")
        for (_prop_dic,index) in ((_F1,1), (_Vol,6), (_Hydro,7), (_pI,8)):
            _prop_dic[l_AA_value[0]] = l_AA_value[index]



## Context-Dependent Substitution Score (CDSS)
## Ref: Koshi, J. M., & Goldstein, R. A. (1995). Context-Dependent Optimal Substitution Matrices. Protein Engineering, 8(7), 641-645. 
CDSS_path = os.path.join(AAS3Drf, "PredMod", "preData", "CDSS")
CDSS_title = ["-","A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
binary_rsa = ["buried", "exposed"]
four_SS = ["beta", "helix", "turn", "coil"]
_CDSS = {}
for rsa in binary_rsa:
    for ss in four_SS:
        context = rsa+"_"+ss
        context_file = os.path.join(CDSS_path, context)
        with open(context_file, "r") as ifile_context:
            n = 0
            for line in ifile_context:
                line = line.rstrip("\n")
                l_value = line.split()
                for index in range(len(l_value)):
                    _CDSS[(rsa, ss, CDSS_title[index], CDSS_title[n])] = l_value[index]
                n += 1



## Environment-dependent Residue Contact Energy (ERCE)
## Ref: Zhang, C., & Kim, S. H. (2000). Environment-dependent residue contact energies for proteins. Proc Natl Acad Sci U S A, 97(6), 2550-2555. doi:10.1073/pnas.040573597
ERCE_path = os.path.join(AAS3Drf, "PredMod", "preData", "ERCE")
ERCE_title = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
contact_list = ["helix_helix", "helix_strand", "helix_coil", "strand_strand", "strand_coil", "coil_coil"]

SS_one = ["H","E","C"]
SS_full = ["helix","strand","coil"]
d_map_ss = dict(zip(SS_full, SS_one))
_ERCE = {}
for contact in contact_list:
    contact_file = os.path.join(ERCE_path, contact+".txt")
    row = contact.split("_")[0]
    col = contact.split("_")[1]
    with open(contact_file, "r") as ifile_contact:
        n = 0
        for line in ifile_contact:
            line = line.rstrip("\n")
            l_energy = line.split()
            for index in range(len(l_energy)):
                _ERCE[(d_map_ss[col], ERCE_title[index], d_map_ss[row], ERCE_title[n])] = l_energy[index]
                _ERCE[(d_map_ss[row], ERCE_title[n], d_map_ss[col], ERCE_title[index])] = l_energy[index]
            n += 1 


