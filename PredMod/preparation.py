#-*- coding:utf-8 -*-

"""
Created on Fri Jan 3 20:38:30 2020 

@author: XiongYao
"""



import os, sys, logging, traceback
from Bio import SwissProt, SeqIO
from Bio.PDB import PDBParser
from xml.etree import ElementTree as ET
from userConfig import AAS3Drf


##*****************##
## Build directory ##
##*****************##
for dir_name in ["UniProtSeq", "StructureSeq"]:
    build_dir = os.path.join(AAS3Drf, "Data", dir_name)
    try:
        if not os.path.exists(build_dir):
            os.makedirs(build_dir)
    except OSError as e:
        pass



##***********************##
## Parse saved gene name ##
##***********************##
gene_file = os.path.join(AAS3Drf, "PredMod", "preData", "All_Human_INFO")
_GeneName = {}
with open(gene_file, "r") as ifile_gene:
    for line in ifile_gene:
        line = line.rstrip("\n")
        l_info = line.split("\t")
        if l_info[2] != "":
            _GeneName[l_info[0]] = l_info[2]
        else:
            _GeneName[l_info[0]] = "-"




##*********************************##
## Check and Complement input data ##
##*********************************##

## 20 standard residues
Lett_three = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
Lett_one = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
d_AA = dict(zip(Lett_three, Lett_one))

    
def parseTXT(txtfile, uniprot_ac):
    """
    Parse flat file
    """
    sequence = ''
    ref_gene = ''
    l_disulfide = []
    try:
        ifiletxt = open(txtfile)
        record = SwissProt.read(ifiletxt)
        sequence = record.sequence
        gene = record.gene_name
        features = record.features
        ifiletxt.close()
        
        if uniprot_ac in _GeneName:
            ref_gene = _GeneName[uniprot_ac]
        else:
            list_gene = []
            l_gene = gene.split(";")
            for item in l_gene:
                if item.find("Name=") != -1:
                    final_gene = item.split("=")[1].split(" ")[0]
                    list_gene.append(final_gene)
            ref_gene = "; ".join(list_gene)
            if ref_gene == "":
                ref_gene = "-"

        for feature in features:
            if feature[0] == "DISULFID":
                if feature[1] == "?" or feature[2] == "?" or str(feature[1]).startswith("?") or str(feature[2]).startswith("?"):
                    pass
                else:
                    if str(feature[1]).startswith("<") or str(feature[1]).startswith(">"):
                        start = int(str(feature[1])[1:])
                    else:
                        start = feature[1]
            
                    if str(feature[2]).startswith("<") or str(feature[2]).startswith(">"):
                        end = int(str(feature[2])[1:])
                    else:
                        end = feature[2]
                    
                    if start != end:
                        l_disulfide.append(start)
                        l_disulfide.append(end)
                    else:
                        l_disulfide.append(start)

    
    except:
        error = traceback.format_exc()
        print "***Error: Please check UniProt flat file***"
        print "***\n"+error+"***"
        print txtfile
    
    return sequence, ref_gene, l_disulfide


def parseXML(xmlfile, uniprot_ac):
    """
    Parse xml file    
    """ 
    sequence = ''
    ref_gene = ''
    l_disulfide = []
    try:
        ifileXML = open(xmlfile)
        record = SeqIO.read(ifileXML, "uniprot-xml")
        sequence = str(record.seq)

        if uniprot_ac in _GeneName:
            ref_gene = _GeneName[uniprot_ac]
        else:
            gene_list = []
            et_tree = ET.parse(xmlfile)
            uniprot_xml = et_tree.getroot()
            entry_xml = uniprot_xml.getchildren()[0]
            for element in entry_xml.getchildren():
                if element.tag == "{http://uniprot.org/uniprot}gene":
                    for children in element.getchildren():
                        if children.attrib['type'] == "primary":
                            gene_list.append(children.text)
            ref_gene = "; ".join(gene_list)
            if ref_gene == "":
                ref_gene = "-"

        features = record.features
        for feature in features:
            if feature.type == "disulfide bond":
                start_pos = feature.location.start.position
                end_pos = feature.location.end.position
                if start_pos == None or end_pos == None:
                    pass
                else:
                    if str(start_pos).startswith("<") or str(start_pos).startswith(">"):
                        start = int(str(start_pos)[1:]) +1
                    else:
                        start = start_pos +1
            
                    if str(end_pos).startswith("<") or str(end_pos).startswith(">"):
                        end = int(str(end_pos)[1:])
                    else:
                        end = end_pos
                    
                    if start != end:
                        l_disulfide.append(start)
                        l_disulfide.append(end)
                    else:
                        l_disulfide.append(start)
                    
    except:
        error = traceback.format_exc()
        print "***Error: Please check UniProt xml file***"
        print "***\n"+error+"***"
        print xmlfile
    
    return sequence, ref_gene, l_disulfide
                

def parsePDB(strufile, stru_id):
    """
    Parse pdb
    """
    try:
        ifile_pdb = PDBParser()
        structure = ifile_pdb.get_structure(stru_id, strufile)
        for model in structure:
            model_id = model.get_id()
            for chain in model:
                chain_id = chain.get_id()
                residues = chain.get_residues()
                l_pos = []
                l_seq = []
                for residue in residues:
                    residue_id = residue.get_id()
                    residue_pos = residue_id[1]
                    hetfield = residue_id[0]
                    if hetfield == " ":
                        l_pos.append(residue_pos)
                        resname = residue.get_resname()
                        l_seq.append(d_AA[resname])
        return l_pos, l_seq, chain_id
    
    except:
        error = traceback.format_exc()
        print "***Error: Please check pdb structure file***"
        print "***\n"+error+"***"
        print strufile
        return "XXXXX","XXXXX","XXXXX"



def getFASTAseq(file_dir, sequence, acc_id):
    filename = acc_id+".fasta"
    filepath = os.path.join(file_dir, filename)
    with open(filepath, "w") as ofile_seq:
        print >> ofile_seq, ">"+acc_id
        while sequence:
            cut = sequence[:75]
            sequence = sequence[75:]
            print >> ofile_seq,cut



def checkStruSeqPosANDaas(input_file, uniprot_dir, file3d_dir, work_path):
    """
    Check whether the sequence parsed from structure file exactly matches UniProt sequence region
    Check 20 standard residues
    Check whether the wild-type amino acid in AAS matches UniProt sequence
    """

    logging.info("Initializing...")
    logging.info("Input data checking...")

    aas_dict = {}
    with open(os.path.join(work_path, input_file), "r") as ifile_aas:
        for line in ifile_aas:
            line = line.rstrip()
            l = line.split("\t")
            key = tuple(l[:3])
            if key not in aas_dict:
                aas_dict[key] = [l[3]]
            else:
                aas_dict[key].append(l[3])

    d_out = {}
    d_disulfide = {}
    for key in aas_dict:
        uniprot_ac, stru_id, chain = key
#        uniprotfile = os.path.join(work_path, uniprot_dir, uniprot_ac+".txt")
#        uniprot_seq, gene_name, list_disulfide = parseTXT(uniprotfile, uniprot_ac)    
        uniprotfile = os.path.join(work_path, uniprot_dir, uniprot_ac+".xml")
        uniprot_seq, gene_name, list_disulfide = parseXML(uniprotfile, uniprot_ac)    
        strufile = os.path.join(work_path, file3d_dir, stru_id+".pdb")
        list_pos, list_aa, pdb_chain = parsePDB(strufile, stru_id)
        if (uniprot_seq == '') or (gene_name == '') or (list_pos == list_aa == pdb_chain == "XXXXX"):
            pass
        else:
            if len(list_pos) != 0:
                pos_start = list_pos[0]
                pos_end = list_pos[-1]
                if (len(list_aa) == pos_end-pos_start+1) and (uniprot_seq[(pos_start-1):pos_end] == "".join(list_aa)) and (chain == pdb_chain):
                    for aas in aas_dict[key]:
                        wt = aas[2]
                        aas_pos = aas[3:-1]
                        if uniprot_seq[int(aas_pos)-1] == wt and wt in Lett_one:
                            new_key = tuple([uniprot_ac, stru_id, gene_name, chain])
                            if new_key not in d_out:
                                d_out[new_key] = [aas]
                            else:
                                d_out[new_key].append(aas)
                        else:
                            print "***Error: Wild-type residue in AAS not match UniProt sequence OR Non-standard wild-type residue***"
                            print "\t".join(key)+"\t"+aas


                    getFASTAseq(os.path.join(AAS3Drf, "Data", "UniProtSeq"), uniprot_seq, uniprot_ac)
                    getFASTAseq(os.path.join(AAS3Drf, "Data", "StructureSeq"), "".join(list_aa), stru_id)
                    d_disulfide[uniprot_ac] = list_disulfide

                else:
                    print "***Error: Please check pdb structure file***"
                    print strufile
            else:
                print "***Error: Please check pdb structure file***"
                print strufile
                
    
    return d_out, d_disulfide
