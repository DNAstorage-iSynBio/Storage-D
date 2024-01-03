#! /usr/bin/python3
# -*- coding:utf-8 -*-

"""
@Author: Wang Yu

@Time: 2023/05/05 13:15:15
"""

import os
import argparse
import re
from collections import Counter
import numpy as np
import time

complementary_table = {"A":"T", "C":"G", "T":"A", "G":"C", "N":"N"}

def get_opts():
    group = argparse.ArgumentParser(description="Cluster and Filter sequences from PE sequencing reads")
    group.add_argument("-p1", "--pair1", help="PE2 fastq file", required=True)
    group.add_argument("-p2", "--pair2", help="PE2 fastq file", required=True)
    group.add_argument("-f", "--filter", help="value of proportion of low quality bases in the sequence", default=0.2, type=float)
    group.add_argument("-c", "--config", help="config path", required=True)
    group.add_argument("-t", "--type", help="type of sequence to assembly: gene or oligo (default=oligo)", default=False, required=False)
    return group

def getTime():
    """
    get current time
    """
    cur_time = time.asctime(time.localtime(time.time()))
    return cur_time

def readToDict(file):
    """
    read file to dict
    """
    _dict = dict()
    f = open(file)
    for i in f:
        i = i.strip().split("\t")
        key, value = i[0], i[1]
        _dict[key] = value
    
    return _dict

def SortDicValue(_dict):
    """
    sort dict key by value
    """
    dict_group = sorted(_dict.items(), key = lambda x:x[1], reverse = True)
    _dict = {i:j for i, j in dict_group}

    return _dict

def DictToFile(_dict, file):
    """
    save dict to file
    """
    f = open(file,'w')
    for key,value in _dict.items():
        f.write(key + "\t" + str(value) + "\n")

def GeneAssembly():
    """ 
    assembly genes sequences
    """
    result_dir = workdir + os.sep + "1_MergeReads"
    os.mkdir(result_dir)
    
    #transform phred33 to phred64
    command_1 = "Phred33_to_Phred64.pl {} > {}.fq\n".format(pe1)
    command_1 += "Phred33_to_Phred64.pl {} > {}.fq\n".format(pe2)
    os.system(command_1)

    #transform fastq file to frg file type
    command_2 = "fastqToCA -libraryname pmw1 -technology illumina-long -type illumina -insertsize 250 20 -mates {}.fq,{}.fq > out.frg".format(pe1, pe2)
    os.system(command_2)

    #run runCA program
    command_3 = "runCA -p genome -d result_gene out.frg"
    os.system(command_3)

    #run spades software
    command_4 = "spades.py -k 21,33,55,77,99,127 --careful --only-assembler -s genome.utg.fasta -o spades"
    os.system(command_4)

    #merge runCA and spades results
    commad_5 = "cat {}/result_gene/genome.asm {}/spades/scaffolds.fasta >{}/merge.fasta".format(result_dir,result_dir,result_dir)

    return result_dir+ os.sep + "merge.fasta"

def MergeReads(pe1, pe2):
    """ 
    merge pair-end reads with pear software
    """
    
    result_dir = workdir + os.sep + "1_MergeReads"
    os.system("mkdir -p {}".format(result_dir))
    result_dir_path = result_dir + os.sep + "reads"

    # run pear command
    command = "pear -f {} -r {} -o {} -y 4G -j 24\n".format(pe1, pe2, result_dir_path)

    log.write("step1: run MergeReads\t{}\n".format(getTime()))


    
    try:
        os.system(command)
    except:
        raise("pear fail run!\n")
    log.write("step1: done!\t{}\n".format(getTime()))

def FilterFaSeq(primer1, primer2):
    """
    select sequences between left primer and right primer
    """
    log.write("step2: filter sequence start\t{}\n".format(getTime()))

    # output file
    filter_dir = workdir + os.sep + file_name + os.sep + "2_FilterSeq"
    try:
        os.mkdir(filter_dir)
    except:
        raise("Please delete old 2_FilterSeq")
    filtered_seq_file = filter_dir + os.sep + "filtered_seq.fasta"
    filtered_seq_f = open(filtered_seq_file,"w")

    # mark duplication sequence and count duplication number
    _dir = workdir + os.sep + os.sep + "1_MergeReads"+ os.sep + "merge.fasta"

    seq_list = list()
    with open(_dir) as f:
        for i in f:
            i = i.strip()
            if i.startswith(">"):
                continue
            else:
                seq_list.append(i)

    seq_count = Counter(seq_list)
    seq_count_dict = dict(seq_count)                            # sequences and count number dict
    seq_count_file = filter_dir + os.sep + "seq.count"
    DictToFile(seq_count_dict, seq_count_file)

    # Match primer to sequences
    primer1_r = ''.join([complementary_table[i] for i in primer2[::-1]])
    primer2_r = ''.join([complementary_table[i] for i in primer1[::-1]])
    pattern_1 = re.compile("(\\S+{})(\\S+){}\\S+".format(primer1[-12:],primer1_r[:12]))
    pattern_2 = re.compile("(\\S+{})(\\S+){}\\S+".format(primer2[-12:],primer2_r[:12]))
    match_dict = dict()

    for seq in seq_count_dict.keys():
        match_seq_1 = re.finditer(pattern_1, seq)
        match_seq_2 = re.finditer(pattern_2, seq)

        seq1_list = [[i.group(2), i.span(2)] for i in match_seq_1]
        seq2_list = [[i.group(2), i.span(2)] for i in match_seq_2]

        try:
            match_dict[seq] = seq1_list[0][0]
        except:
            pass
        try:
            match_dict[seq] = seq2_list[0][0][::-1]
        except:
            pass

    count = 0
    seq_list = list(match_dict.values())
    seq_count = Counter(seq_list)
    com_length = len(seq_count.most_common(1)[0][0])
    filtered_seq_list = []
    for i in seq_list:
        if len(i) == com_length:
            filtered_seq_list.append(i)
            filtered_seq_f.write(">seq_" +str(count) + "\n" + i+ "\n")
            count += 1

    filtered_seq_count = Counter(seq_list)
    filtered_seq_count_dict = dict(filtered_seq_count)                            # sequences and count number dict
    filtered_seq_count_file = filter_dir + os.sep + "filtered_seq.count"
    DictToFile(filtered_seq_count_dict, filtered_seq_count_file)

    log.write("step2: done \t{}\n".format(getTime()))
    filtered_seq_f.close()
    return(filtered_seq_file)

def FilterSeq(primer1, primer2,file_name):
    """
    select sequences between left primer and right primer
    """
    log.write("step2: filter sequence start\t{}\n".format(getTime()))

    # output file
    _dir = workdir + os.sep + "2_FilterSeq"
    filter_dir = _dir + os.sep + file_name
    try:
        os.system("mkdir -p {}".format(filter_dir))
    except:
        raise("Please delete old 2_FilterSeq")

    filter_low_qual = filter_dir + os.sep + "low_qual_seq.tab"
    filter_low_qual_f = open(filter_low_qual,"w")
    filter_low_qual_f.write("sequence\tcount\tstart_pos\tend_pos\n")

    # mark duplication sequence and count duplication number
    _dir = workdir + os.sep + "1_MergeReads" + os.sep + "reads.assembled.fastq"
    print(_dir)
    command = "awk '{{if(NR%4 == 2 || NR%4 == 0) print}}' {} > {}/seq_qual.txt".format(_dir,filter_dir)
    os.system(command)

    # Match primer to sequences
    primer1_r = ''.join([complementary_table[i] for i in primer2[::-1]])
    primer2_r = ''.join([complementary_table[i] for i in primer1[::-1]])
    pattern_1 = re.compile("(\\S+{})(\\S+){}\\S+".format(primer1[int(len(primer1)/2):],primer1_r[:int(len(primer1_r)/2)]))
    pattern_2 = re.compile("(\\S+{})(\\S+){}\\S+".format(primer2[int(len(primer2)/2):],primer2_r[:int(len(primer2_r)/2)]))
    match_dict = dict()

    merged_list = []
    flag = 0
    seq_qual_dict = dict()
    merged_seq_qual_dict = dict()       #assembled read sum of qual within same read
    with open(filter_dir+os.sep+"seq_qual.txt") as f:
        for i in f:
            i = i.strip()
            if flag % 2 == 0:
                key = i
                merged_list.append(i)
                flag += 1
            else:
                seq_qual_dict[key] = i
                flag += 1
                if key not in merged_seq_qual_dict.keys():
                    merged_seq_qual_dict[key] = np.array([ord(q) for q in i])
                else:
                    merged_seq_qual_dict[key] += np.array([ord(q) for q in i])

    merged_seq_count = Counter(merged_list)
    merged_seq_count = dict(merged_seq_count)
    match_seq_count_dict = dict()
    for seq in set(merged_list):
        match_seq_1 = re.finditer(pattern_1, seq)
        match_seq_2 = re.finditer(pattern_2, seq)

        seq1_list = [[i.group(2), i.span(2)] for i in match_seq_1]
        seq2_list = [[i.group(2), i.span(2)] for i in match_seq_2]
        if seq1_list != []:
            try:
                _seq1, start, stop = seq1_list[0][0], seq1_list[0][1][0], seq1_list[0][1][1]
                match_dict[_seq1] += merged_seq_qual_dict[seq][start:stop]
                match_seq_count_dict[_seq1] += merged_seq_count[seq]
            except:
                _seq1, start, stop = seq1_list[0][0], seq1_list[0][1][0], seq1_list[0][1][1]
                match_dict[_seq1] = merged_seq_qual_dict[seq][start:stop]
                match_seq_count_dict[_seq1] = merged_seq_count[seq]
        elif seq2_list != []:
            try:
                _seq2, start, stop = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]]), seq2_list[0][1][0], seq2_list[0][1][1]
                match_dict[_seq2] += merged_seq_qual_dict[seq][start:stop][::-1]
                match_seq_count_dict[_seq2] += merged_seq_count[seq]
            except:
                _seq2, start, stop = "".join([complementary_table[i] for i in seq2_list[0][0][::-1]]), seq2_list[0][1][0], seq2_list[0][1][1]
                match_dict[_seq2] = merged_seq_qual_dict[seq][start:stop][::-1]
                match_seq_count_dict[_seq2] = merged_seq_count[seq]
        else:
            continue
    
    match_seq_qual_dict = dict()
    for key in match_dict.keys():
        match_seq_qual_dict[key] = "".join([chr(int(q/match_seq_count_dict[key])) for q in match_dict[key]])

    seq_count_file = filter_dir + os.sep + "matched_seq.count"
    DictToFile(match_seq_count_dict, seq_count_file)

    seq_list = []
    # filter low quality sequence
    for key, value in match_seq_qual_dict.items():
        if len([i for i in value if i< "5"]) > len(value) * proportion:
            filter_low_qual_f.write(key + "\t" + str(value) +"\n")
            match_seq_count_dict.pop(key)
        else:
            seq_list.append(key)

    match_seq_count_dict = SortDicValue(match_seq_count_dict)
    com_length = len(list(match_seq_count_dict.keys())[0])
    print("--------------length: ", com_length)
    for i in seq_list:
        if len(i) != com_length:
            match_seq_count_dict.pop(i)

    match_seq_count_dict = SortDicValue(match_seq_count_dict)
    filtered_seq_count_file = filter_dir + os.sep + "filtered_seq.count"
    DictToFile(match_seq_count_dict, filtered_seq_count_file)

    log.write("step2: done \t{}\n".format(getTime()))
    filter_low_qual_f.close()

def Blast(filtered_seq_file):
    """
    blast filtered sequence
    """
    log.write("step3: BLAST start \t{}\n".format(getTime()))


    # blasti
    _dir = workdir + os.sep + "3_Blast"
    blast_dir = _dir + os.sep + file_name

    try:
        os.system("mkdir -p {}".format(blast_dir))
    except:
        raise("can't mkdir blast_dir")

    blast_out = blast_dir + os.sep + "blast.out"

    pre_blast_seq = blast_dir + os.sep + "ref.fa"
    ref_file = open(pre_blast_seq,'w')
    seq_dict = dict()
    flag = 0
    with open(filtered_seq_file) as f:
        for i in f:
            i = i.strip()
            seq = i.split("\t")[0]
            id = ">seq_" + str(flag)
            ref_file.write(id + "\n" + seq + "\n")
            seq_dict[id] = seq
            flag += 1
    f.close()
    ref_file.close()

    command = "makeblastdb -in {} -dbtype nucl\n".format(pre_blast_seq)
    command += "blastn -query {} -db {} -out {} -outfmt 6 -num_threads 24".format(pre_blast_seq, pre_blast_seq, blast_out)
    print(command)
    os.system(command) 
    log.write("step3: blast done \t{}\n".format(getTime()))

    # select
    log.write("step4: select start \t{}\n".format(getTime()))
    
    seq_out = list()
    flag = ""
    with open(blast_out) as f:
        for i in f:
            i = i.strip().split("\t", 4)
            seq1, seq2 = i[0], i[1]
            seq2_n = int(seq2.replace("seq_", ""))
            if seq1 != flag:
                try:
                    seq_out.append(min(group_list))
                except:
                    pass
                group_list = [seq2_n]
                flag = seq1
            else:
                group_list.append(seq2_n)
    
    f.close()
    seq_out.append(min(group_list))
    seq_out = set(seq_out)
    seq_out = sorted(list(seq_out))
    select_file = blast_dir + os.sep + "select_seq.fasta"
    with open(select_file, 'w') as f:
        for i in seq_out:
            id = ">seq_" + str(i)
            f.write(">seq_" + str(i) + "\n" + seq_dict[id] + "\n")
    f.close()
    log.write("step4: select done \t{}\n".format(getTime()))        

if __name__ == "__main__":
    group = get_opts()
    opts = group.parse_args()
    pe1 = opts.pair1
    pe2 = opts.pair2
    proportion = opts.filter
    config = opts.config
    type = opts.type

    p = os.path.split(pe1)[1]
    workdir = os.getcwd() + os.sep + "result_" + p.split("_")[0]
    log = open("result.log","w")
    
    if not type:
        # 1. merge pair-end reads with pear software
        MergeReads(pe1, pe2)
    elif type == "gene":
        GeneAssembly(pe1,pe2)
    
    config_files = os.listdir(config)
    for file in config_files:
        file_name = file.replace(".config", "")
        with open(config + os.sep + file) as f:
            for i in f:
                log.write(file_name+"\tstart run: \t{}\n".format(getTime()))
                i = i.strip()
                if i.startswith("leftPrimer"):
                    primer1 = i.split("\t")[1]
                elif i.startswith("rightPrimer"):
                    primer2 = i.split("\t")[1]

            # 2. select sequences between left primer and right primer
            if not type:
                filtered_seq_file = FilterSeq(primer1, primer2, file_name)
            elif type == "gene":
                filtered_seq_file = FilterFaSeq(primer1, primer2, file_name)
    
            # 3. blast filtered sequence
            filtered_seq_file = workdir + os.sep + "2_FilterSeq" + os.sep + file_name + os.sep + "filtered_seq.count"
            Blast(filtered_seq_file)
            log.write(file_name+"\tdone ! \n")

    log.close()
