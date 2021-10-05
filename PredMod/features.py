#-*- coding:utf-8 -*-

"""
Created on Sun Jan 5 10:40:25 2020 

@author: XiongYao
"""

import os, sys, logging, math
import numpy as np
from numpy import mean, var
from Bio.Blast import NCBIXML
from userConfig import AAS3Drf
from data import _AA_Freq, _rvis, _miyata, _F1, _Vol, _Hydro, _pI, _CDSS, _ERCE


## 20 standard residues
Lett_three = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
Lett_one = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
d_AA = dict(zip(Lett_three, Lett_one))


def parseNaccess(stru_id, chain):
    rsa_file = os.path.join(AAS3Drf, "Data", "Naccess", stru_id+".rsa")
    _rsa = {}
    if chain == " ":
        pos_index = 2
        value_index = 4
    else:
        pos_index = 3
        value_index = 5

    try:
        with open(rsa_file, "r") as ifile_rsa:
            for line in ifile_rsa:
                line = line.rstrip("\n")
                l_info = line.split()
                if l_info[0] == "RES":
                    key = (stru_id, l_info[pos_index], d_AA[l_info[1]])
                    value = l_info[value_index]
                    _rsa[key] = value
                else:
                    pass
    except IOError:
        pass
    return _rsa


def parseDSSP(stru_id, chain):
    dssp_file = os.path.join(AAS3Drf, "Data", "DSSP", stru_id+".dssp")
    dssp_ss = ["H","I","G","S","T","B","E"]
    _dssp = {}
    
    if chain == " ":
        value_index = 3
    else:
        value_index = 4

    try:
        with open(dssp_file, "r") as ifile_dssp:
            list_line = ifile_dssp.readlines()
            list_line = [line.rstrip("\n") for line in list_line]
            for line in list_line:
                if line.find("#") != -1:
                    start_index = list_line.index(line) +1
            for line_index in range(start_index, len(list_line)):
                l_info = list_line[line_index].split()
                key = (stru_id, l_info[1])
                value = l_info[value_index]
                if value not in dssp_ss:
                    value = "C"
                _dssp[key] = value
    except IOError:
        pass
                
    return _dssp        


def parseRingEntropy(stru_id):
    node_file = os.path.join(AAS3Drf, "Data", "RING", stru_id+".node")
    _entropy = {}

    try:
        with open(node_file, "r") as ifile_node:
            node_title = ifile_node.readline()
            for line in ifile_node:
                line = line.rstrip("\n")
                l_info = line.split("\t")
                for index in range(len(l_info)):
                    if l_info[index] == " ":
                        l_info[index] = "-"
                key = (stru_id, l_info[2], d_AA[l_info[3]])
                value = l_info[13]
                _entropy[key] = value
    except IOError:
        pass

    return _entropy
        
def parseRingContact(stru_id, chain):
    edge_file = os.path.join(AAS3Drf, "Data", "RING", stru_id+".edge")
    _contact = {}

    if chain == " ":
        start = "_:"
    else:
        start = chain+":"

    try:
        with open(edge_file, "r") as ifile_edge:
            for line in ifile_edge:
                if line.startswith(start):
                    line = line.rstrip("\n")
                    l_info = line.split()
                    node1 = l_info[0]
                    node2 = l_info[2]
                    node1_split = node1.split(":")
                    node2_split = node2.split(":")
                    if (stru_id, node1_split[1], d_AA[node1_split[3]]) not in _contact:
                        _contact[(stru_id, node1_split[1], d_AA[node1_split[3]])] = [(stru_id, node2_split[1], d_AA[node2_split[3]])]
                    else:
                        _contact[(stru_id, node1_split[1], d_AA[node1_split[3]])].append((stru_id, node2_split[1], d_AA[node2_split[3]]))
                    if (stru_id, node2_split[1], d_AA[node2_split[3]]) not in _contact:
                        _contact[(stru_id, node2_split[1], d_AA[node2_split[3]])] = [(stru_id, node1_split[1], d_AA[node1_split[3]])]
                    else:
                        _contact[(stru_id, node2_split[1], d_AA[node2_split[3]])].append((stru_id, node1_split[1], d_AA[node1_split[3]]))
    except IOError:
        pass

    return _contact



def parseDDG(stru_id):
    DDG_filename = "Dif_"+stru_id+".fxout"
    DDG_file = os.path.join(AAS3Drf, "Data", "FoldX", "BuildModel_OUT", DDG_filename)
    _out = {}

    try:
        with open(DDG_file, "r") as ifile_DDG:
            for line in ifile_DDG:
                line = line.rstrip("\n")
                if line.startswith(stru_id):
                    l_info = line.split()
                    aas_name = os.path.splitext(l_info[0])[0]
                    aas_num = aas_name[(aas_name.rfind("_")+1):]
                    _out[(stru_id, aas_num)] = l_info[1]
                else:
                    pass
    except IOError:
        pass
            
    mut_filename = "individual_list_"+stru_id+".mut"
    foldx_mut_file = os.path.join(AAS3Drf, "Data", "FoldX", "structure_mut", mut_filename)
    _DDG = {}
    try:
        with open(foldx_mut_file, "r") as ifile_mut:
            n = 1
            for line in ifile_mut:
                mut = line.rstrip("\n")[:-1]
                _DDG[stru_id, mut] = _out[(stru_id, str(n))]
                n += 1
    except:
        pass

    return _DDG



def addRVIS(gene_name):
    list_gene = gene_name.split("; ")
    list_rvis = []
    for gene in list_gene:
        if gene in _rvis:
            list_rvis.append(_rvis[gene])
        else:
            pass

    if len(list_rvis) != 0:
        rvis_score = str(sum(list_rvis)/len(list_rvis))
    else:
        rvis_score = "-"

    return rvis_score 


def addDDG(stru_id, wt, mt, pos, chain, d_DDG):
    if chain == " ":
        mut = wt+pos+mt
    else:
        mut = wt+chain+pos+mt

    if (stru_id, mut) in d_DDG:
        return d_DDG[(stru_id, mut)]
    else:
        return "-"


def addCDSS(stru_id, wt, mt, pos, rsa_score, d_dssp):
    dssp_ss = ["H","I","G","S","T","B","C","E"]
    full_ss = ["helix","helix","turn","turn","turn","coil","coil","beta"]
    d_state = dict(zip(dssp_ss, full_ss))

    if rsa_score == "-" or (stru_id, pos) not in d_dssp:
        cdss_score = "-"
    else:
        if float(rsa_score) > 18.0:
            flag_rsa = "exposed"
        else:
            flag_rsa = "buried"

        ss_state = d_dssp[(stru_id, pos)]
        cdss_score = _CDSS[(flag_rsa, d_state[ss_state], wt, mt)]
    return cdss_score


def addERCE(stru_id, wt, mt, pos, d_contact, d_dssp):
    dssp_ss = ["H","I","G","S","T","B","C","E"]
    ERCE_state = ["H","C","C","C","C","C","C","E"]
    d_state = dict(zip(dssp_ss, ERCE_state))

    if (stru_id, pos, wt) in d_contact and (stru_id, pos) in d_dssp:
        l_wt_ERCE = []
        l_partner = d_contact[(stru_id, pos, wt)]
        ss_state = d_dssp[(stru_id, pos)]
        wt_ss = d_state[ss_state]
        for (i, j, k) in l_partner:
            partner_ss = d_state[d_dssp[(i, j)]]
            l_wt_ERCE.append(float(_ERCE[(wt_ss, wt, partner_ss, k)]))
        erce_score = str(sum(l_wt_ERCE))
    else:
        erce_score = "-"

    return erce_score


def addNeighborPhy(stru_id, wt, mt, pos, d_contact, d_prop, feature):
    key_name = (stru_id, pos, wt)
    if key_name in d_contact:
        l_neighbor = []
        for (i, j, k) in d_contact[key_name]:
            l_neighbor.append(float(d_prop[k]))
        if feature == "d_Mean_Mutant":
            neighbor_mean = mean(l_neighbor)
            out_score = str(neighbor_mean - float(d_prop[mt]))
        if feature == "var_wild_NB":
            l_wild = l_neighbor[:]
            l_wild.append(float(d_prop[wt]))
            out_score = str(var(l_wild))
    else:
        out_score = "-"

    return out_score


def addSeqDis(pos, uniprot_ac, d_disulfide):
    cys_pos = d_disulfide[uniprot_ac]
    list_dis = []
    for position in cys_pos:
        dis = abs(pos - position)
        list_dis.append(dis)
    
    if len(list_dis) == 0:
        least_dis = "50"
    else:
        sort_dis = sorted(list_dis)
        least_dis = sort_dis[0]
        if least_dis > 50:
            least_dis = "50"
        else:
            least_dis = str(least_dis)

    return least_dis    


def getBLASTposFeatures(wt, mt, position, features, blast, flag):
    try:
        aas = blast['aas'][:, position-1]
        scores = blast['scores'][:, position-1]
        evalue = blast['evalues'][:, position-1]
    except:
        aas = np.array([], dtype=np.character)
        scores = np.array([], dtype=np.float)
        evalue = np.array([], dtype=np.float)
    
    aas_evalue = aas[evalue < 10**-45]
    weights_evalue = scores[evalue < 10**-45]
    aas_human = aas[blast['human']]
    weights_human = scores[blast['human']]
    aas_not_human = aas[np.logical_not(blast['human'])]
    weights_not_human = scores[np.logical_not(blast['human'])]

    blast_out = []
    for feature in features:
        feature_parts = feature.split("_")
        subset, aas_feature = feature_parts[:-1]
        
        if subset == "all":
            weights = scores if feature_parts[-1] == "w" else None
            fea_score = getAASfeature(wt, mt, aas_feature, aas, weights)

        elif subset == "eva":
            weights = weights_evalue if feature_parts[-1] == 'w' else None
            fea_score = getAASfeature(wt, mt, aas_feature, aas_evalue, weights)

        elif subset == "hum":
            weights = weights_human if feature_parts[-1] == 'w' else None
            fea_score = getAASfeature(wt, mt, aas_feature, aas_human, weights)

        elif subset == "nhu":
            weights = weights_not_human if feature_parts[-1] == 'w' else None
            fea_score = getAASfeature(wt, mt, aas_feature, aas_not_human, weights)

        blast_out.append(str(fea_score))
    
    if flag == 0:
        return ["-"]*8
    else:
        return blast_out


def getAASfeature(wt, mt, feature, aas, weights):
    if weights is None:
        weights = np.repeat(1, len(aas))
    
    def n_aln():
        return weights.sum()
    
    def wt_aas():
        return float(((aas == wt)*weights).sum())
    
    def mt_aas():
        return float(((aas == mt)*weights).sum())
    
    if feature == 'nwt':
        return wt_aas()
    elif feature == 'nmt':
        return mt_aas()
    elif feature == 'naa':
        return float(((aas != '-')*weights).sum())
    elif feature == 'nal':
        return n_aln()
    elif feature == 'rwt':
        return wt_aas() / len(aas) if len(aas) > 0 else 0.0
    elif feature == 'rmt':
        return mt_aas() / len(aas) if len(aas) > 0 else 0.0
    elif feature == 'pwm':
        return _pwm(wt, mt, wt_aas(), mt_aas())
       

def _pwm(wt, mt, wt_x, mt_x):
    try:
        return (math.log(mt_x/_AA_Freq[mt]) - 
                math.log(wt_x/_AA_Freq[wt]))
    except ValueError:
        return -10.


def getUniProtSEQ(seqfile):
    with open(seqfile, "r") as ifile_seq:
        list_seq = ifile_seq.readlines()
        list_seq = [line.rstrip("\n") for line in list_seq]
    seq = "".join(list_seq[1:])
    return seq

               
def parseBLASTxml(uniprot_ac):
    aas = []
    scores = []
    evalues = []
    human = []
    flag = 1

    blastXMLfile = os.path.join(AAS3Drf, "Data", "UniProtXML", uniprot_ac+".xml")
    uniprotSEQfile = os.path.join(AAS3Drf, "Data", "UniProtSeq", uniprot_ac+".fasta")
    sequence = getUniProtSEQ(uniprotSEQfile)

    try:
        blast_xml = open(blastXMLfile)
        records = NCBIXML.parse(blast_xml)
        for record in records:
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    new_aas = np.array(["-"]*len(sequence))
                    new_scores = np.full(len(sequence), 0.0)
                    new_evalue = np.full(len(sequence), 1.0)
                
                    start = hsp.query_start -1
                    end = hsp.query_end

                    if hsp.query == hsp.sbjct == sequence:
                        pass
                    else:
                        new_aas[start:end] = [s for q,s in zip(hsp.query, hsp.sbjct) if q != "-"]
                        new_scores[start:end] = hsp.bits
                        new_evalue[start:end] = hsp.expect
                        new_human = "Homo sapiens" in alignment.title
                        aas.append(new_aas)
                        scores.append(new_scores)
                        evalues.append(new_evalue)
                        human.append(new_human)

        blast_xml.close()
    except IOError:
        flag = 0
    except:
        pass

    return {'aas': np.array(aas),
            'scores': np.array(scores),
            'evalues': np.array(evalues),
            'human': np.array(human, dtype=np.bool)}, flag

    

def addFeatures(d_AAS, d_disulfide):
    logging.info("Feature extracting...")

    fea_list = ["all_nwt_no","nhu_nal_no","eva_nal_no","all_naa_w","all_pwm_w","hum_naa_w","eva_nwt_w","eva_rmt_w"] 
    outfile = open(os.path.join(AAS3Drf, "Data", "USER_DATA.tsv"),"w")
    out_title = ["UniProt_AC","Structure_ID","Chain","AAS","nwt","nal_no_human","nal_1e-45","naa_w","pwm_w","naa_human_w",
                 "nwt_1e-45_w","rmt_1e-45_w","Entropy_RING","RVIS","Miyata","Abs_dF1","Abs_dpI","SeqDis_DISULFID","dF1_Mean_Mutant",
                 "dVol_Mean_Mutant","varHP_Wild_NB","ERCE","CDSS","RSA","DDG"]
    print >> outfile,"\t".join(out_title)

    for key in d_AAS:
        uniprot_ac, stru_id, gene_name, chain = key
        d_MSA, flag = parseBLASTxml(uniprot_ac)
        d_rsa = parseNaccess(stru_id, chain)
        d_dssp = parseDSSP(stru_id, chain)
        d_entropy = parseRingEntropy(stru_id)
        d_contact = parseRingContact(stru_id, chain)
        d_DDG = parseDDG(stru_id)
        RVIS = addRVIS(gene_name)
        for variant in d_AAS[key]:
            aas, wt, mt, pos = variant[2:], variant[2], variant[-1], variant[3:-1]
            list_con = getBLASTposFeatures(wt, mt, int(pos), fea_list, d_MSA, flag)
            Entropy_RING = d_entropy[(stru_id, pos, wt)] if (stru_id, pos, wt) in d_entropy  else "-"
            miyata_score = _miyata[(wt, mt)]
            abs_dF1 = str(abs(float(_F1[mt]) - float(_F1[wt])))
            abs_dpI = str(abs(float(_pI[mt]) - float(_pI[wt])))
            seqdis_cys = addSeqDis(int(pos), uniprot_ac, d_disulfide)
            dF1_mean_mutant = addNeighborPhy(stru_id, wt, mt, pos, d_contact, _F1, "d_Mean_Mutant")
            dVol_mean_mutant = addNeighborPhy(stru_id, wt, mt, pos, d_contact, _Vol, "d_Mean_Mutant")
            varHP_wild_nb = addNeighborPhy(stru_id, wt, mt, pos, d_contact, _Hydro, "var_wild_NB")
            ERCE_score = addERCE(stru_id, wt, mt, pos, d_contact, d_dssp)
            rsa_score = d_rsa[(stru_id, pos, wt)] if (stru_id, pos, wt) in d_rsa  else "-"
            CDSS_score = addCDSS(stru_id, wt, mt, pos, rsa_score, d_dssp)
            DDG_score = addDDG(stru_id, wt, mt, pos, chain, d_DDG)
            print >> outfile,"\t".join([uniprot_ac, stru_id, chain, variant])+"\t"+"\t".join(list_con)+"\t"+Entropy_RING+\
                     "\t"+RVIS+"\t"+miyata_score+"\t"+abs_dF1+"\t"+abs_dpI+"\t"+seqdis_cys+"\t"+dF1_mean_mutant+"\t"+\
                     dVol_mean_mutant+"\t"+varHP_wild_nb+"\t"+ERCE_score+"\t"+CDSS_score+"\t"+rsa_score+"\t"+DDG_score

    outfile.close()
