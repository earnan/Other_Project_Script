#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   3951_codon_change.py
#         Author:   yujie
#    Description:   3951_codon_change.py
#        Version:   1.0
#           Time:   2022/03/21 11:58:48
#  Last Modified:   2022/03/21 11:58:48
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
from xml.dom.minidom import Element
from Bio import SeqIO
import os
import re
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   3951_codon_change')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='snp-file', type=str, default="F:/3951答疑/global_align/archive/Hibiscus_sabdariffa_var_altissima_vs_Hibiscus_cannabinus_MK404537.1.snp.xls", required=False)
optional.add_argument('-s', '--sample',
                      metavar='[file]', help='sample-fasta', type=str, default="F:\\3951答疑\\global_align\\archive\\Hibiscus_sabdariffa_cds.fasta", required=False)
optional.add_argument('-r', '--ref',
                      metavar='[file]', help='ref-fasta', type=str, default="F:\\3951答疑\\global_align\\archive\\Hibiscus_cannabinus_cds.fasta", required=False)
optional.add_argument('-snote',
                      metavar='[file]', help='sample-note', type=str, default="F:\\3951答疑\\global_align\\archive\\Hibiscus_sabdariffa_var_altissima.ann", required=False)
optional.add_argument('-rnote',
                      metavar='[file]', help='ref-note', type=str, default="F:\\3951答疑\\global_align\\archive\\Hibiscus_cannabinus.ann", required=False)
optional.add_argument(
    '-f', '--flag',  metavar='[flag]', help='flag', type=str, default='Y', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

#################################################################
# 格式化成2016-03-20 11:45:39形式
begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))
#################################################################
# 前置脚本 from_gbk_get_cds


def read_fasta_to_dic(infasta):  # 增加位置的字典
    with open(infasta, 'r') as f:
        seq_id = ''  # 基因名
        dict_seq = {}  # 基因名-序列
        dict_len = {}  # 基因名-长度
        dict_pos = {}  # 基因名-位置
        d_pos = {}
        for line in f:
            seq_pos = []  # 某基因对应的位置组成的列表
            s_pos = []
            l_n = [0]  # 计算同名基因是第几个
            if line.startswith('>'):
                seq_id = line.strip('\n').split()[2].split('=')[
                    1].strip(']')  # 基因名
                if seq_id in dict_seq.keys():
                    l_n.append(0)
                    seq_id = seq_id+'-'+str(len(l_n))  # 基因名+1 ycf1-2形式
                l_tmp = re.findall(
                    r'\d+', line.strip('\n').split()[1].lstrip(
                        '[').rstrip(']'))  # .sort()  # 位置打散成一个个起点或终点
                l_tmp2 = []
                [l_tmp2.append(int(i)) for i in l_tmp]  # 全部转换成数字,放进l_tmp2,未排序
                [s_pos.append(i) for i in l_tmp2]
                d_pos[seq_id] = s_pos
                l_tmp2.sort()
                # print(l_tmp2)
                [seq_pos.append(i) for i in l_tmp2]
                dict_pos[seq_id] = seq_pos
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{0} Item Quantity: {1} {2} {3}'.format(os.path.basename(infasta),
                                                  len(dict_seq), len(dict_len), len(dict_pos)))
    return dict_seq, dict_len, dict_pos, d_pos


def judgment_section(number, list):  # 判断是否位于该基因上
    s = ''
    if len(list) == 6:
        if (number >= list[0] and number <= list[1]) or (number >= list[2] and number <= list[3]) or (number >= list[4] and number <= list[5]):
            s = 'Y'
    elif len(list) == 4:
        if (number >= list[0] and number <= list[1]) or (number >= list[2] and number <= list[3]):
            s = 'Y'
    elif len(list) == 2:
        if (number >= list[0] and number <= list[1]):
            s = 'Y'
    else:
        s = False
    return s


def read_file_to_dic(infile, s_dict_pos):  # 把snp结果读成字典然后与已有字典比较,判断位于哪个基因上
    with open(infile, 'r') as f:
        seq_id = ''
        d_point = {}
        d_base = {}
        tmp_list = []
        f.readline()  # 跳过第一行
        for line in f:
            if len(line.split()) == 8:
                # print(line.split())
                # print(line.split()[2])
                s_point_pos = int(line.split()[2])
                s_base = line.split()[3]
                r_point_pos = int(line.split()[5])
                r_base = line.split()[6]
                if judgment_section(s_point_pos, s_dict_pos[line.split()[7]]):
                    seq_id = line.split()[7]
                elif (line.split()[7]+'-2') in s_dict_pos.keys() and (judgment_section(s_point_pos, s_dict_pos[line.split()[7]+'-2'])):
                    seq_id = line.split()[7]+'-2'
                """
                else:
                    if (line.split()[7]+'-2') in s_dict_pos.keys():
                        print(s_point_pos, line.split()[
                              7], s_dict_pos[line.split()[7]], line.split()[
                              7]+'-2', s_dict_pos[line.split()[7]+'-2'])
                    else:
                        print(s_point_pos, line.split()[
                              7], s_dict_pos[line.split()[7]])
                    print(
                        'seq_id will be set to {} y/n'.format(line.split()[7]))
                    if args.flag == 'Y':
                        seq_id = line.split()[7]
                    else:
                        flag = input()
                        if flag == 'y' or flag == 'Y':
                            seq_id = line.split()[7]
                        else:
                            seq_id = flag
                    print(seq_id)
                """
                # print(seq_id)
                if seq_id in d_point.keys():
                    d_point[seq_id].append([s_point_pos, r_point_pos])
                    d_base[seq_id].append([s_base, r_base])
                else:
                    d_point[seq_id] = []
                    d_base[seq_id] = []
                    d_point[seq_id].append([s_point_pos, r_point_pos])
                    d_base[seq_id].append([s_base, r_base])
    print(len(d_point), len(d_base))
    """
                if tmp_id != seq_id:
                    tmp_list.append([tmp_s_point_pos, tmp_r_point_pos])
                    d_point[tmp_id] = tmp_list
                    tmp_list = []
                else:
                    tmp_list.append([s_point_pos, r_point_pos])

                tmp_s = s_point_pos
                tmp_r = r_point_pos
                print(len(d_point))
            tmp_id = seq_id
            tmp_s_point_pos = tmp_s
            tmp_r_point_pos = tmp_r
    """
    return d_point, d_base


(s_dict_seq, s_dict_len, s_dict_pos, s_d_pos) = read_fasta_to_dic(args.sample)
(r_dict_seq, r_dict_len, r_dict_pos, r_d_pos) = read_fasta_to_dic(args.ref)
# print(s_dict_pos)
# print('\n')
# print(s_d_pos)  # 原始顺序
# print(r_dict_pos)

d_point, d_base = read_file_to_dic(args.infile, s_dict_pos)
# print(d_point)
# print('\n')
# print(d_base)
# s_dict_seq, s_dict_len, s_dict_pos名字里有-2形式
# 根据草图修改
# def find_codon(d_point, s_dict_seq, s_dict_pos):
# for ele in d_point.keys()


def judg_n(n, l):  # judg_n(655,[333,1000])#判断是否位于该区间
    if (n >= l[0] and n <= l[1]) or (n <= l[0] and n >= l[1]):
        s = True
    else:
        s = False
    return s


def judgment_segmentation(number, list):  # 原始顺序,判断属于基因几个区间中的哪一段
    if len(list) == 6:
        if judg_n(number, [list[0], list[1]]):
            s = 1
        elif judg_n(number, [list[2], list[3]]):
            s = 2
        elif judg_n(number, [list[4], list[5]]):
            s = 3
    elif len(list) == 4:
        if judg_n(number, [list[0], list[1]]):
            s = 1
        elif judg_n(number, [list[2], list[3]]):
            s = 2
    elif len(list) == 2:
        if judg_n(number, [list[0], list[1]]):
            s = 1
    else:
        s = False
    return s  # s=1-[0]   2-[2] 3-[4]


def trans2acid(codon):  # 翻译成氨基酸
    """
    The Bacterial and Plant Plastid Code (11):
    Stnd    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    This    AAs = FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
    Starts      = ---M---------------M------------MMMM---------------M------------
    Base1       = TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
    Base2       = TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
    Base3       = TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG
    """
    genetic_code_number = 11
    acid = ''
    code_table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
                  'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                  'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                  'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', }
    acid = code_table[codon]
    return acid


def find_codon(point_pos, d_pos, key, seq):  # 具体查找密码子
    s = judgment_segmentation(point_pos, d_pos[key])  # 处于第几段
    # print(s)
    codon = ''
    if s == 1:
        gene_start = d_pos[key][0]
        n = 1+abs(gene_start-point_pos)  # n代表cds中位置
        # print(n)
        if n % 3 == 2:
            codon = seq[n-2:n+1]
        elif n % 3 == 0:
            codon = seq[n-3:n]
        elif n % 3 == 1:
            codon = seq[n-1:n+2]
    elif s == 2:
        gene_start = d_pos[key][2]
        # print(gene_start)
        n = (1+abs(d_pos[key][0]-d_pos[key][1]))+1+abs(gene_start-point_pos)
        # print(n)
        if n % 3 == 2:
            codon = seq[n-2:n+1]
        elif n % 3 == 0:
            codon = seq[n-3:n]
        elif n % 3 == 1:
            codon = seq[n-1:n+2]
    elif s == 3:
        gene_start = d_pos[key][4]
        # print(gene_start)
        n = (1+abs(d_pos[key][0]-d_pos[key][1])) + \
            (1+abs(d_pos[key][2]-d_pos[key][3]))+1+abs(gene_start-point_pos)
        # print(n)
        if n % 3 == 2:
            codon = seq[n-2:n+1]
        elif n % 3 == 0:
            codon = seq[n-3:n]
        elif n % 3 == 1:
            codon = seq[n-1:n+2]
    # print(codon)
    return codon, n


def read(file):
    dic = {}
    with open(file, 'r') as f:
        for line in f:
            l_n = [0]
            if line.startswith('CDS'):
                seq_id = line.strip('\n').split()[2]  # 基因名
                if seq_id in dic.keys():
                    l_n.append(0)
                    seq_id = seq_id+'-'+str(len(l_n))  # 基因名+1 ycf1-2形式
                dic[seq_id] = line.strip('\n').split()[1]
    print(len(dic))
    return dic


r_dic = read(args.rnote)
s_dic = read(args.snote)

p = 0
q = 0
print('Gene:\tPosition:\tRef_base--Sample_base:\tCodon_mutate:\taa_mutate:\tMutate_type:\tRef_location:\tSample_location:\t')
for key in d_point.keys():  # 循环控制查找
    p += 1
    # if key == 'ndhB-2':
    # print(key)
    s_seq = s_dict_seq[key]
    r_seq = r_dict_seq[key]
    for i in range(len(d_point[key])):
        q += 1
        s_point_pos = d_point[key][i][0]  # 第i组snp中s的位点
        r_point_pos = d_point[key][i][1]
        s_base = d_base[key][i][0]
        r_base = d_base[key][i][1]

        s_codon, s_n = find_codon(s_point_pos, s_d_pos, key, s_seq)
        s_acid = trans2acid(s_codon)
        r_codon, r_n = find_codon(r_point_pos, r_d_pos, key, r_seq)
        r_acid = trans2acid(r_codon)
        if s_acid == r_acid:
            flag = 'Synonymous'
        else:
            flag = 'Nonsynonymous'
        #print('Gene:{0} Position:{1}↔{2} Ref_base↔Sample_base:{3}↔{4} Codon_mutate:{5}↔{6} aa_mutate:{7}↔{8} Mutate_type:{9} Position_Start:{10} Position_End:{11}'.format(key, r_point_pos, s_point_pos, r_base, s_base, r_codon, s_codon, r_acid, s_acid, flag, 0, 0))
        print(key+'\t'+str(r_point_pos)+'--'+str(s_point_pos)+'\t'+r_base+'--' + s_base + '\t' +
              r_codon+'--'+s_codon+'\t'+r_acid+'--'+s_acid+'\t'+flag+'\t['+str(r_dic[key])+']\t['+str(s_dic[key])+']')
        #print('样本 cds位置{0}: {1} {2}  参考 cds位置{3}: {4} {5}'.format(            s_n, s_codon, s_acid, r_n, r_codon, r_acid))


print('{0}个基因共计{1}个snp位点'.format(p, q))
###############################################################
print('\n')
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
###############################################################
