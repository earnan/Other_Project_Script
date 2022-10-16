#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   3structure_2str.py
#         Author:   yujie
#    Description:   3structure_2str.py
#        Version:   1.0
#           Time:   2022/09/02 15:04:54
#  Last Modified:   2022/09/02 15:04:54
#        Contact:   hi@arcsona.cn
#        License:   GNU General Public License v3.0
#
##########################################################
from Bio import SeqIO
from Bio.Seq import Seq
#from icecream import ic
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
#       Filename:   3structure_2str.py\n\
#         Author:   yujie\n\
#    Description:   3structure_2str.py\n\
#        Version:   1.0\n\
#           Time:   2022/09/02 15:05:04\n\
#  Last Modified:   2022/09/02 15:05:04\n\
#        Contact:   hi@arcsona.cn\n\
#        License:   Copyright (C) 2022\n\
#\n\
##########################################################\n\
\n\
\npython3   3structure_2str.py\n\
Function:\n\
1. -i1 [aln fa ]  -i2 [group list ]  -o1 [new aln fa ]   -o2 [new group list ]\n\
\n\
##########################################################\n\
Path: E:\OneDrive\jshy信息部\Script\Other_Project_Script\GP-20220722-4689\3structure_2str.py\n\
Path: /share/nas1/yuj/script/chloroplast/assembly/3structure_2str.py\n\
Version: 1.0\n\
##########################################################\n\
'
)
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i1', '--infile1', metavar='[aln fa]', help='老师提供的比对好的文件', type=str, default='F:\\结果\\1_21个牦牛群体_基因交流值\\21_sets_seq.aln.fas', required=False)
optional.add_argument(
    '-i2', '--infile2', metavar='[group list]', help='老师提供的分组', type=str, default='F:\\结果\\1_21个牦牛群体_基因交流值\\group.list', required=False)
optional.add_argument(
    '-o1', '--outfile1', metavar='[new_aln_fa]', help='样品按名字排序后的比对文件', type=str, default='F:\\4689\\result\\3_Structure\\sortesd_21_sets_seq.aln.fas', required=False)
optional.add_argument(
    '-o2', '--outfile2', metavar='[new group list]', help='生成文件,第一列样品名  第二列分组名', type=str, default='F:\\Project\\4689\\result\\1_gene_flow\\group_list.txt', required=False)
optional.add_argument('-h', '--help', action='help', help='[help_information]')
args = parser.parse_args()

aln_fa_file_path = args.infile1
group_list_file_path = args.infile2
out_new_aln_fa = args.outfile1
new_group_list_path = args.outfile2


def format_fasta(seq, num):  # 格式化字符串,
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
                    # list_all.append(prefix+str(j))
            elif len(tmp_list) == 2:
                # 说明无版本号
                for j in range(int(tmp_list[0]), int(tmp_list[1])+1):
                    list_all.append(prefix+str(j))
            else:
                print('-------------------XXXXXXXXXXXXXXXXXXXXXXX')
        else:
            list_all.append(i.strip())  # 记得去除空格
            # list_all.append(i.strip().rstrip('.1'))  # 记得去除空格和.1
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
                group_accession_dict[group_name] = list_all

                group_taxset_dict[group_name] = '{}-{}'.format(
                    start+1, start+int(Numbers))  # TaxSet {t5} = 9-10;
                start = start+int(Numbers)
    return group_accession_dict, group_taxset_dict


# ###############################################################################################
formatted_accession_seq_dict, int_ntax, str_seq_len = get_formatted_seq_dict(
    aln_fa_file_path)
group_accession_dict, group_taxset_dict = geneflow_get_group_sets_dict(
    group_list_file_path)
# #############################################################################################
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

# ######################################################生成id被排序的比对文件
with open(out_new_aln_fa, 'w') as out_handle:
    for i in group_accession_dict.keys():
        acc_list = group_accession_dict[i]  # 一组登录号
        for j in formatted_accession_seq_dict.keys():  # 键有括号
            if j.split('(')[0] in acc_list:  # 表明属于这一组
                out_handle.write('>'+j+'\n')
                out_handle.write(formatted_accession_seq_dict[j].strip()+'\n')


# #########################################################生成 第一列样品名  第二列分组名 文件

with open(new_group_list_path, 'w') as out_handle:
    for i in group_accession_dict.keys():
        acc_list = group_accession_dict[i]  # 一组登录号
        for j in formatted_accession_seq_dict.keys():  # 键有括号
            if j.split('(')[0] in acc_list:  # 表明属于这一组
                out_handle.write(j+'\t'+i+'\n')
