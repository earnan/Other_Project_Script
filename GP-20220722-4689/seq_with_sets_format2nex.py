#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   seq_with_sets_format2nex.py
#         Author:   yujie
#    Description:   seq_with_sets_format2nex.py
#        Version:   1.0
#           Time:   2022/08/29 15:59:53
#  Last Modified:   2022/08/29 15:59:53
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
#import pretty_errors
import re
import sys
import time
#import copy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(
    add_help=False, usage='\n\
\n\
##########################################################\n\
#\n\
#       Filename:   seq_with_sets_format2nex.py\n\
#         Author:   yujie\n\
#    Description:   seq_with_sets_format2nex.py\n\
#        Version:   1.0\n\
#           Time:   2022/08/29 16:00:13\n\
#  Last Modified:   2022/08/29 16:00:13\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   seq_with_sets_format2nex.py\n\
Function:\n\
1.Geneflow -gf -i [aln.fa] -g [group_list] -o [nex] \n\
2.Mismatch_distribution -mm -i [all.aln.fa] -g [group_list] -o [nex] \n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\Other_Project_Script\GP-20220722-4689\seq_with_sets_format2nex.py\n\
Path: /share/nas1/yuj/script/Other_Project_Script/GP-20220722-4689/seq_with_sets_format2nex.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--inaln', metavar='[aln.fa]', help='aln.fa', type=str, default='F:/4689/21_sets_seq.aln.fas', required=False)
optional.add_argument(
    '-g', '--grouplist', metavar='[group.list]', help='分组文件,不通用', type=str, default='F:/4689/group_mapping.txt', required=False)
optional.add_argument(
    '-o', '--outnex', metavar='[out nex]', help='outfile', type=str, default='F:/4689/21_sets_seq.nex', required=False)
optional.add_argument('-gf', '--geneflow', help='默认否,运行则-gf',
                      action='store_true', required=False)
optional.add_argument('-mm', '--mismatch', help='默认否,运行则-mm',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[help_information]')
args = parser.parse_args()

# 文件路径
aln_fa_file_path = args.inaln  # i
group_list_file_path = args.grouplist  # g
nex_out_path = args.outnex  # o


def format_fasta(seq, num):  # 格式化字符串
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        if (index + 1) % num == 0:  # 可以用来换行
            format_seq += "\n"
    return format_seq


# ##################################################################################nex的CHARACTERS部分,读取比对文件,返回formatted_accession_seq_dict
def get_formatted_seq_dict(aln_fa_file_path):
    with open(aln_fa_file_path, 'r') as fa_handle:
        # CHARACTERS部分,nex存放具体序列的部分,在中间
        all_seq_dict = {}  # 直接读取文件 登录号(有括号信息)  :  序列,无换行
        formatted_accession_seq_dict = {}  # 登录号(有括号信息)  :  序列,换行
        previous_current_note = ''  # 上一个id
        previous_final_note = ''  # 上一个id
        all_seq_dict[previous_current_note] = ''  # 赋空值,第一次要写的序列为空
        format_seq = ''
        int_ntax = 0  # 序列总个数

        for line in fa_handle:
            if line.startswith('>'):
                '''写入序列'''
                int_ntax += 1
                format_seq = format_fasta(
                    all_seq_dict[previous_current_note], 70)
                if format_seq != '':  # 由于是先写序列,后写id,直接写入,那么第一行会空出一行
                    # 第一次要写的序列为空,正好可以判断为空就不写入
                    formatted_accession_seq_dict[previous_final_note] = format_seq+'\n'+'\n'
                '''读取成字典'''
                current_note = line.strip().lstrip('>')  # 有括号信息的登录号
                previous_current_note = current_note
                all_seq_dict[current_note] = ''
                previous_final_note = current_note
                formatted_accession_seq_dict[previous_final_note] = ''
                # note = "'{}'{}".format(current_note, '   ')  # 后面有3空格的id
            else:
                all_seq_dict[current_note] += line.strip().upper()

        str_seq_len = len(all_seq_dict[previous_current_note])
        # 第509个序列  最后一个序列,要单独拿出来存入字典
        format_seq = format_fasta(all_seq_dict[previous_current_note], 70)
        formatted_accession_seq_dict[previous_final_note] = format_seq+'\n'+'\n'
    return formatted_accession_seq_dict, int_ntax, str_seq_len


# ####################################################1.geneflow子函数,对每组字符串进行解析,获得该组所有序列号
def geneflow_get_all_accession(tmp_flg, group_dict, tmp_count):
    tmp_str = group_dict[tmp_flg][0]  # NC1-NC5,MK1,MW3 原始信息
    tmp_number = group_dict[tmp_flg][1]  # 22 25  24  个数
    list1 = tmp_str.split(',')  # 没有修改过的登录号列表,后续要进行一些修改
    list_all = []
    for i in list1:
        if i.find('~') >= 0:
            prefix = re.findall(r'[A-Z]+', i)[0]  # 匹配大写字母  登录号前缀
            tmp_list = re.findall(r'\d+', i)  # 匹配数字
            if len(tmp_list) == 4:
                # 说明有版本号
                for j in range(int(tmp_list[0]), int(tmp_list[2])+1):
                    list_all.append(prefix+str(j)+'.1')
            elif len(tmp_list) == 2:
                # 说明无版本号
                for j in range(int(tmp_list[0]), int(tmp_list[1])+1):
                    list_all.append(prefix+str(j))
            else:
                print('-------------------XXXXXXXXXXXXXXXXXXXXXXX')
        else:
            list_all.append(i.strip())  # 记得去除空格
    if len(list_all) == int(tmp_number):
        tmp_count += int(tmp_number)
    return tmp_count, list_all  # 完全没涉及到括号问题


def geneflow_get_group_sets_dict(group_list_file_path):
    with open(group_list_file_path, 'r', encoding='utf-8') as group_handle:
        group_dict = {}  # 初始
        group_accession_dict = {}  # A : MK   NC1  NC2
        group_taxset_dict = {}  # TaxSet {t5} = 9-10;
        tmp_count = 0  # 总计数,每次都有累加
        start = 0
        for line in group_handle:
            if not line.startswith('Variety'):
                # 以下切片
                line_content = line.strip().split('\t')
                group_name = line_content[-1]
                Numbers = line_content[2]
                Accession_number = line_content[3]
                # 赋值
                group_taxset_dict[group_name] = ''  # TaxSet {t5} = 9-10;
                group_dict[group_name] = [
                    Accession_number, Numbers]  # 登录号   个数,初始信息

                tmp_count, list_all = geneflow_get_all_accession(
                    group_name, group_dict, tmp_count)  # 解析
                group_accession_dict[group_name] = list_all  # 登录号列表

                group_taxset_dict[group_name] = '{}-{}'.format(
                    start+1, start+int(Numbers))  # TaxSet {t5} = 9-10;
                start = start+int(Numbers)
    return group_accession_dict, group_taxset_dict


# #################################################2.错配分布Mismatch distribution子函数
def mismatch_get_group_sets_dict(group_list_file_path):
    with open(group_list_file_path, 'r', encoding='utf-8') as group_handle:
        group_accession_dict, group_taxset_dict = {}, {}
        all_accession_list = []
        start = 0
        int_ntax_group = 0
        for line in group_handle:
            if not line.startswith('>'):
                group_taxset_dict[line.split()[0]] = ''
                group_accession_dict[line.split()[0]] = ''
                group_taxset_dict[line.split()[0]] = '{}-{}'.format(
                    start+1, start+int(line.split()[1]))  # TaxSet {t5} = 9-10;
                start = start+int(line.split()[1])
            elif line.startswith('>'):
                int_ntax_group += 1
                all_accession_list.append(line.strip().lstrip('>'))
    for i in group_accession_dict.keys():
        start_index = int(group_taxset_dict[i].split('-')[0])-1
        end_index = int(group_taxset_dict[i].split('-')[1])
        group_accession_dict[i] = all_accession_list[start_index:end_index]
    return group_accession_dict, group_taxset_dict, int_ntax_group


# ################################################################################输出的nex里需要用的变量,返回这几个str前缀变量,然后写入
if args.geneflow and args.inaln and args.grouplist and args.outnex:  # 基因流
    formatted_accession_seq_dict, int_ntax, str_seq_len = get_formatted_seq_dict(
        aln_fa_file_path)
    group_accession_dict, group_taxset_dict = geneflow_get_group_sets_dict(
        group_list_file_path)
    str_ntax = str(int_ntax)
    str_matrix = ''
    str_taxlabels = ''
    for i in group_accession_dict.keys():  # A B C D E
        # a组对应的登录号,无括号信息,因为分组文\文件里是范围,压根没写括号
        acc_list = group_accession_dict[i]
        for j in formatted_accession_seq_dict.keys():  # 键全为登录号,有括号
            if j.split('(')[0] in acc_list:  # 表明属于这一组
                str_taxlabels += ("'"+j+"'"+'\n')
                str_matrix += ("'"+j+"'"+'   '+formatted_accession_seq_dict[j])
    str_taxlabels = str_taxlabels.strip()
    str_sets = ''
    for k, v in group_taxset_dict.items():
        str_sets += '   TaxSet {} = {};\n'.format(k, v)

elif args.mismatch and args.inaln and args.grouplist and args.outnex:  # 错配分布
    formatted_accession_seq_dict, int_ntax, str_seq_len = get_formatted_seq_dict(
        aln_fa_file_path)
    group_accession_dict, group_taxset_dict, int_ntax_group = mismatch_get_group_sets_dict(
        group_list_file_path)
    if int_ntax_group != int_ntax:
        print('数量不对应')
    str_ntax = str(int_ntax)
    str_matrix = ''
    str_taxlabels = ''
    for i in group_accession_dict.keys():  # A B C D E
        print(i)
        acc_list = group_accession_dict[i]  # a组对应的登录号,有括号信息
        for j in formatted_accession_seq_dict.keys():  # 键全为登录号,有括号
            if j in acc_list:  # 表明属于这一组
                str_taxlabels += ("'"+j+"'"+'\n')
                str_matrix += ("'"+j+"'"+'   '+formatted_accession_seq_dict[j])
    str_taxlabels = str_taxlabels.strip()
    str_sets = ''
    for k, v in group_taxset_dict.items():
        str_sets += '   TaxSet {} = {};\n'.format(k, v)


s = "#NEXUS\n\
[File generated by DnaSP Ver. 6.12.03, from file: {0}    Aug 1, 2022]\n\
\n\
BEGIN TAXA;\n\
DIMENSIONS NTAX={1};\n\
TAXLABELS\n\
{2};\n\
END;\n\
\n\
BEGIN CHARACTERS;\n\
DIMENSIONS NCHAR={3};\n\
FORMAT DATATYPE=DNA  MISSING=? GAP=- ;\n\
MATRIX\n\
{4}\n\
;\n\
END;\n\
\n\
BEGIN SETS;\n\
{5}\
END;\n\
\n\
BEGIN CODONS;\n\
   CODESET * UNTITLED = Universal: all;\n\
END;\n\
\n\
BEGIN CODONUSAGE;\n\
END;\n\
\n\
BEGIN DnaSP;\n\
   Genome= Haploid;\n\
   ChromosomalLocation= Mitochondrial;\n\
   VariationType= DNA_Seq_Pol;\n\
   Species= ---;\n\
   ChromosomeName= ---;\n\
   GenomicPosition= 1;\n\
   GenomicAssembly= ---;\n\
   DnaSPversion= Ver. 6.12.03;\n\
END;\n\
".format(aln_fa_file_path, str_ntax, str_taxlabels, str_seq_len, str_matrix.rstrip(), str_sets)

with open(nex_out_path, 'w') as nex_handle:
    nex_handle.write(s)

# ############################################################范例
'''
NEXUS
[File generated by DnaSP Ver. 6.12.03, from file: {}    Aug 1, 2022]

BEGIN TAXA
DIMENSIONS NTAX = {str_ntax}

TAXLABELS
{str_taxlabels};
END;

BEGIN CHARACTERS;
DIMENSIONS NCHAR={str_seq_len};
FORMAT DATATYPE=DNA  MISSING=? GAP=- ;
MATRIX
{str_matrix}
;
END;

BEGIN SETS;
{str_sets}
END;

BEGIN CODONS;
   CODESET * UNTITLED = Universal: all;
END;

BEGIN CODONUSAGE;
END;

BEGIN DnaSP;
   Genome= Haploid;
   ChromosomalLocation= Autosome;
   VariationType= DNA_Seq_Pol;
   Species= ---;
   ChromosomeName= ---;
   GenomicPosition= 1;
   GenomicAssembly= ---;
   DnaSPversion= Ver. 6.12.03;
END;
'''
