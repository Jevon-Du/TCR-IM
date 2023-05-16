#!/usr/bin python
# _*_ coding:utf-8 _*_
#
# @Version : 1.0
# @Time    : 2023/2/10 11:20
# @Author  : Du Jiewen
# @File    : process_netMHC_MHC_seq.py
#
# 清洗netMHC pseudosequence
# IPD MHC sequence: 下载 -> 存储 -> 读取

import csv
from datetime import datetime
from pathlib import Path
import time
import os
import sys
import re
from typing import Union, Optional, List

import pandas as pd
import requests

def download_MHC_sequence(MHC_locus: str, 
                          saved_dir: Optional[Path], 
                          url: str = 'https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/{}_prot.fasta') -> Optional[str]:
    """
    Download standard MHC protein sequences from IPD-MHC Database[https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/]
    
    Args:
        MHC_locus: MHC locus name.
        url: eg. https://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/fasta/{A}_prot.fasta
        saved_dir: Path, optional. If saved_dir passed, save MHC sequence file to this path.
        
    Returns:
        MHC protein sequence.
    """
    resp = requests.get(url.format(MHC_locus), timeout=300)
    if resp.status_code == 200:
        if saved_dir:
            print(f"Save {MHC_locus}_prot.fasta to {saved_dir}")
            with open(saved_dir/f'{MHC_locus}_prot.fasta', 'w') as f:
                f.write(resp.text)
        else:
            return resp.text
    else:
        print(f"Download failed: status code {r.status_code}\n{r.text}")

        
def read_fasta(reader, sep='>'):
    item = []
    for line in reader:
        if sep in line:
            yield item
            item = []
        item.append(line.strip())
    yield item
    
def load_IPD_MHC_sequence(MHC_locus : str, 
                          data_dir : Path) -> List[list]:
    """
    Read MHC sequence from IPD-MHC Database
    
    Args:
        MHC_locus: which fasta files need to be read.
        data_dir: the directory where MHC protein fasta files are located.
    
    Returns:
        list: a list of MHC seuqence items.
    """
    IPD_MHCs = []
    
    file = data_dir/f'{MHC_locus}_prot.fasta'
    fasta_generator = read_fasta(open(file, 'r'))
    _ = next(fasta_generator)
    
    if MHC_locus != 'MHC':
        for item in fasta_generator:
            allele_name_4_digits = item[0].split(' ')[1]
            digits = allele_name_4_digits.split(':')
            if ('DR' or 'DQ' or 'DP') in allele_name_4_digits:
                allele_name_2_digits = digits[0] + ':' + digits[1]
                if len(digits)>2:
                    allele_name_3_digits = digits[0] + ':' + digits[1] + ':' + digits[2]
                else:
                    allele_name_3_digits = None
            else:
                allele_name_2_digits = 'HLA-' + digits[0] + ':' + digits[1]
                if len(digits)>2:
                    allele_name_3_digits = 'HLA-' + digits[0] + ':' + digits[1] + ':' + digits[2]
                else:
                    allele_name_3_digits = None
            seq = ''.join(item[1:])
            IPD_MHCs.append([allele_name_2_digits, allele_name_3_digits, allele_name_4_digits, seq])
    else:       
        species = ['BoLA-', 'DLA-', 'SLA-', 'Patr-', 'Gogo-', 'Mamu-']
        for item in fasta_generator:
            allele_name_4_digits = item[0].split(' ')[1]
            for x in species:
                if x in allele_name_4_digits:
                    digits = allele_name_4_digits.split(':')
                    try:
                        allele_name_2_digits = digits[0] + ':' + digits[1]
                        if len(digits)>2:
                                allele_name_3_digits = digits[0] + ':' + digits[1] + ':' + digits[2]
                        else:
                            allele_name_3_digits = None

                        seq = ''.join(item[1:])
                        IPD_MHCs.append([allele_name_2_digits, allele_name_3_digits, allele_name_4_digits, seq])
                    except:
                        # TODO(Du Jiewen): 非典型MHC allele处理： Gogo-DPB1*02, Gogo-DPB1*03, Gogo-DPB1*05, Gogo-H*01
                        # Patr-DOB*01, Patr-H*01, DLA-DRB1*09701
                        print(allele_name_4_digits)
                        break
    
    return IPD_MHCs


def load_netMHC_MHCI_data(MHC_pseudo_file: Path) -> list:
    """
    """
    mhc1_pseudos = []
    with open(MHC_pseudo_file, 'r') as f:
        for line in f:
            line = re.split('[\s]+', line.strip())
            mhc1_pseudos.append(line)
    
    def process_hla(x):
        if ':' not in x:
            return x[:5]+'*'+x[5:-2]+':'+x[-2:]
        elif x[5]==':':
            pass
        else:
            return x[:5]+'*'+x[5:]
    
    std_HLA_1 = [[x[0], process_hla(x[0]), x[1]] for x in mhc1_pseudos if 'HLA' in x[0]]
    
    def process_BoLA(x):
        if x.startswith('BoLA-N:'):
            return  ''
        else:
            return x.replace(':', '*')[:-2]+':'+x[-2:]

    std_BoLA_1 = [[x[0], process_BoLA(x[0]), x[1]] for x in mhc1_pseudos if ('BoLA' in x[0]) and (':' in x[0])]   
    
    def process_Patr(x):
        x = x.replace(':', '')
        return x[:6]+'*'+x[-4:-2]+':'+x[-2:]

    std_Patr_1 = [[x[0], process_Patr(x[0]), x[1]] for x in mhc1_pseudos if 'Patr' in x[0]]
    
    std_Gogo_1 = [[x[0], f'{x[0][:6]}*{x[0][6:8]}:{x[0][8:]}', x[1]] for x in mhc1_pseudos if 'Gogo' in x[0]]

    std_DLA_1 = [[x[0], f'{x[0][:6]}*{x[0][6:9]}:{x[0][9:]}', x[1]] for x in mhc1_pseudos if 'DLA' in x[0]]
    
    std_H2_1 = [[x[0], x[0].replace('H-2', 'H2'), x[1]] for x in mhc1_pseudos if ('H2' in x[0]) or ('H-2' in x[0])]
    
    def process_SLA(x):
        x = x.replace(':', '')
        if x[4:].isdigit():
            return x[:5]+'*'+x[-4:-2]+':'+x[-2:]
        else:
            return x        

    std_SLA_1 = [[x[0], process_SLA(x[0]), x[1]] for x in mhc1_pseudos if 'SLA' in x[0]]     

    def process_Mamu(x):
        if '*' in x:
            return  x[:-2]+':'+x[-2:]
        elif x[7]==':':
            if (not x.startswith('Mamu-A1')) and x[-1]=='0':
                x = x[:-1]
                # print(x)
            return  x[:7]+'*'+x[8:-2]+':'+x[-2:]
        elif x.startswith('Mamu-B'):
            if ':' not in x:
                if x[8:]:
                    a = x[:6]+'*'+x[6:8]+':'+x[8:]
                else:
                    a = x[:6]+'*'+x[-2:]
                return a
            elif x[6]==':':
                return x.replace(':', '*')[:-2]+':'+x[-2:]
            else:
                return 'Mamu-B*'+x.replace('Mamu-B', '')
        else:
            return x

    std_Mamu_1 = [[x[0], process_Mamu(x[0]), x[1]] for x in mhc1_pseudos if 'Mamu' in x[0]]

    # fix
    fixed_alleles = [['Mamu-A2*00102', 'Mamu-A2*01:02'], ['Mamu-A7*00103', 'Mamu-A7*01:03']]
    # Mamu-A2*01:02:01:01, Mamu-A7*01:03:01:01
    # 删除没有*和:的
    # IPD: Mamu-B*066:01:01:02
    # 除过A1外，A2~A7的格式为：Mamu-A2*12:01:01:02
    std_Mamu_1_fix = []

    for row in std_Mamu_1:
        if row[0]=='Mamu-A2*00102':
            row[1] = 'Mamu-A2*01:02'
        elif row[0]=='Mamu-A7*00103':
            row[1] = 'Mamu-A7*01:03'
        elif row[0]=='Mamu-AG:01':
            continue
        std_Mamu_1_fix.append(row)
    
    ## merge every species
    MHCI_pseudoseqs = []
    MHCI_pseudoseqs.extend(std_HLA_1)
    MHCI_pseudoseqs.extend(std_BoLA_1)
    MHCI_pseudoseqs.extend(std_Gogo_1)
    MHCI_pseudoseqs.extend(std_Patr_1)
    MHCI_pseudoseqs.extend(std_DLA_1)
    MHCI_pseudoseqs.extend(std_H2_1)
    MHCI_pseudoseqs.extend(std_SLA_1)
    MHCI_pseudoseqs.extend(std_Mamu_1_fix)

    all_mhc_1_names = set()
    MHCI_pseudoseqs_valid = []
    for row in MHCI_pseudoseqs:
        if ('*' in row[1]) and (':' in row[0]): 
            if row[1] not in all_mhc_1_names:
                all_mhc_1_names.add(row[1])
                MHCI_pseudoseqs_valid.append(row)
    
    return MHCI_pseudoseqs_valid

def merge_MHCI_seqs(netMHC_MHCI_seqs: list, 
                    IPD_MHCI_seqs: list) -> list:
    """
    Merge netMHC pseudo sequence with IPD sequence
    
    
    """
    ## 
    MHCI_seqs = []

    for netMHC_x, netMHC_standard_x, pseduoseq in netMHC_MHCI_seqs:
        for IPD_x in IPD_MHCI_seqs:
            if (netMHC_standard_x+':01:01' == IPD_x[2]):
                MHCI_seqs.append([netMHC_x, netMHC_standard_x, *IPD_x, pseduoseq])
                break
            elif (netMHC_standard_x+':01' == IPD_x[1]):
                MHCI_seqs.append([netMHC_x, netMHC_standard_x, *IPD_x, pseduoseq])
                break
            elif (netMHC_standard_x in IPD_x[0]):
                MHCI_seqs.append([netMHC_x, netMHC_standard_x, *IPD_x, pseduoseq])
                break
    
    return MHCI_seqs

def fix_HLA_II(x):
    if x[-1].isalpha():
        x = x[:-1]
    if 'DR' in x:
        # DRB5_0102
        # DRB5_0108N -> DRB5*01:08:01N 删除
        return ['DRA*01:01', f'{x[:4]}*{x[5:7]}:{x[7:]}']
    else:
        # HLA-DPA10103-DPB10101 -> 
        DPA = x.split('-')[1]
        DPB = x.split('-')[2]
        DPA_standard = DPA[:4]+'*'+DPA[4:-2]+':'+DPA[-2:]
        DPB_standard = DPB[:4]+'*'+DPB[4:-2]+':'+DPB[-2:]
        return [DPA_standard, DPB_standard]

def load_netMHC_MHCII_data(MHC_pseudo_file: Path) -> list:
    MHCII_pseudos = []
    with open(MHC_pseudo_file, 'r') as f:
        for line in f:
            line = re.split('[\s]+', line.strip())
            if 'BoLA' in line[0]:
                MHCII_pseudos.append([line[0], None, f'{line[0][:9]}*0{line[0][10:-2]}:{line[0][-2:]}', line[1]])
            elif 'H-2' in line[0]:
                MHCII_pseudos.append([line[0], line[0], None, line[1]])
            else:
                MHCII_pseudos.append([line[0], *fix_HLA_II(line[0]), line[1]])
    
    return MHCII_pseudos

def merge_MHCII_seqs(netMHC_MHCII_seqs: list, 
                     IPD_MHCII_seqs: list) -> list:
    """
    Merge netMHC pseudo sequence with IPD sequence
    
    
    """
    # HLA: DR只有DRB1,B3,B4,B5单链; DP:DPA1, DPB1双链; DQ:DQA1, DQB1双链

    # 扁平化 standard_netMHC_HLA_II_allele
    flatten_mhc2_allele = []
    for x in netMHC_MHCII_seqs:
        if x[1]:
            flatten_mhc2_allele.append(x[1])
        if x[2]:
            flatten_mhc2_allele.append(x[2])

    matched_mhc2_allele = {}
    for x in flatten_mhc2_allele:
        for IPD_x in IPD_MHCII_seqs:
            if (x +':01:01' == IPD_x[2]):
                matched_mhc2_allele[x] = [IPD_x[-2], IPD_x[-1]]
                break
            elif (x +':01' == IPD_x[1]):
                matched_mhc2_allele[x] = [IPD_x[-2], IPD_x[-1]]
                break
            elif (x in IPD_x[0]):
                matched_mhc2_allele[x] = [IPD_x[-2], IPD_x[-1]]
                break
    
    std_mhc2_with_IPD = []
    for x in netMHC_MHCII_seqs:
        std_mhc2_with_IPD.append([x[0], x[1], *matched_mhc2_allele.get(x[1], [None, None]), x[2], *matched_mhc2_allele.get(x[2], [None, None]), x[-1]])
    
    for x in netMHC_MHCII_seqs:
        if 'H-2' in x[0]:
            std_mhc2_with_IPD.append([x[0], None, None, None, None, None, None, x[-1]])
    
    return std_mhc2_with_IPD
        
if __name__ == '__main__':
    DATA_DIR = Path('/home/jiewen.du/data')
    IPD_dir = DATA_DIR / 'MHC_sequence'
    MHC_locus = ['A', 'B', 'C', 'E', 'F', 'G', 'DRA', 'DRB1', 'DRB345', 'DPA1', 'DPB1', 'DQA1', 'DQA2', 'DQB1', 'MHC']
    
    # load MHC sequence file exists.
    IPD_MHC_seqs = []
    for mhc in MHC_locus:
        if not (IPD_dir/f'{mhc}_prot.fasta').exists():
            download_MHC_sequence(mhc, IPD_dir)
        seqs = load_IPD_MHC_sequence(mhc, IPD_dir)
        IPD_MHC_seqs.extend(seqs)
    
    # load netMHC sequence
    netMHC_MHCI_seqs = load_netMHC_MHCI_data(DATA_DIR/'netMHC/NetMHCpan-4.1/NetMHCpan_train/MHC_pseudo.dat')
    netMHC_MHCII_seqs = load_netMHC_MHCII_data(DATA_DIR/'netMHC/NetMHCIIpan-4.0/NetMHCIIpan_train/pseudosequence.2016.all.X.dat')
    
    # merge sequence
    standard_MHCI_seqs = merge_MHCI_seqs(netMHC_MHCI_seqs, IPD_MHC_seqs)
    standard_MHCII_seqs = merge_MHCII_seqs(netMHC_MHCII_seqs, IPD_MHC_seqs)
    
    with open(DATA_DIR/'standard_MHCI_allele_sequence.csv', 'w') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(['MHCI_allele_netMHC', 'MHCI_allele_standard', 'MHCI_allele_IPD', 'MHCI_sequence', 'pseudosequence_netMHC'])
        for row in standard_MHCI_seqs:
            f_csv.writerow([row[0], row[1], row[-3], row[-2], row[-1]])
    
    MHCII_file = DATA_DIR / 'standard_MHCII_allele_sequence.csv'
    MHCII_file_header = ['MHCII_allele_netMHC', 'MHCII_A_allele_standard', 'MHCII_A_allele_IPD', 'MHCII_A_sequence', 
                          'MHCII_B_allele_standard', 'MHCII_B_allele_IPD', 'MHCII_B_sequence', 'pseudosequence_netMHC']
    with open(MHCII_file, 'w') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(MHCII_file_header)
        f_csv.writerows(standard_MHCII_seqs)
    
