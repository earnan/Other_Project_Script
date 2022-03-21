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
from Bio import SeqIO
import os
import re
import time

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   3951_codon_change')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--infile',
                      metavar='[file]', help='snp-file', type=str, default="F:/3951答疑/global_align/fasta/Hibiscus_sabdariffa_var_altissima_vs_Hibiscus_cannabinus_MK404537.1.snp.xls", required=False)
optional.add_argument('-s', '--sample',
                      metavar='[file]', help='sample-fasta', type=str, default="F:\\3951答疑\\global_align\\fasta\\Hibiscus_sabdariffa_cds.fasta", required=False)
optional.add_argument('-r', '--ref',
                      metavar='[file]', help='ref-fasta', type=str, default="F:\\3951答疑\\global_align\\fasta\\Hibiscus_cannabinus_cds.fasta", required=False)
# optional.add_argument('--input',
# metavar='[input]', help='input', type=str, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


begin_time = time.time()
start_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('Start Time : {}'.format(start_time))
# 前置脚本 cp_from_gbk_get_cds


def read_fasta_to_dic(infasta):  # 增加位置的字典
    with open(infasta, 'r') as f:
        seq_id = ''
        dict_seq = {}
        dict_len = {}
        dict_pos = {}
        for line in f:
            seq_pos = []
            l_n = [0]
            if line.startswith('>'):
                seq_id = line.strip('\n').split()[2].split('=')[1].strip(']')
                if seq_id in dict_seq.keys():
                    l_n.append(0)
                    seq_id = seq_id+'-'+str(len(l_n))
                l_tmp = line.strip('\n').split()[1].lstrip(
                    '[').rstrip(']').split('..')
                [seq_pos.append(int(i)) for i in l_tmp]
                dict_pos[seq_id] = seq_pos
                dict_seq[seq_id] = ''
                dict_len[seq_id] = ''
            else:
                dict_seq[seq_id] += line.strip('\n')
                dict_len[seq_id] += str(len(line.strip('\n')))
    print('{0} Item Quantity: {1} {2} {3}'.format(os.path.basename(infasta),
                                                  len(dict_seq), len(dict_len), len(dict_pos)))
    return dict_seq, dict_len, dict_pos


def read_file_to_dic(infile, s_dict_pos):
    with open(infile, 'r') as f:
        seq_id = ''
        d_point = {}
        f.readline()
        tmp_list = []
        for line in f:
            if len(line.split()) == 8:
                # print(line.split())
                # print(line.split()[2])
                s_point_pos = int(line.split()[2])
                r_point_pos = int(line.split()[5])

                if s_point_pos >= s_dict_pos[line.split()[7]][0] and s_point_pos <= s_dict_pos[line.split()[7]][1]:
                    seq_id = line.split()[7]
                elif s_point_pos <= s_dict_pos[line.split()[7]][0] and s_point_pos >= s_dict_pos[line.split()[7]][1]:
                    seq_id = line.split()[7]
                elif (line.split()[7]+'-2') in s_dict_pos.keys() and (s_point_pos >= s_dict_pos[line.split()[7]+'-2'][0] and s_point_pos <= s_dict_pos[line.split()[7]+'-2'][1]):
                    seq_id = line.split()[7]+'-2'
                elif (line.split()[7]+'-2') in s_dict_pos.keys() and (s_point_pos <= s_dict_pos[line.split()[7]+'-2'][0] and s_point_pos >= s_dict_pos[line.split()[7]+'-2'][1]):
                    seq_id = line.split()[7]+'-2'
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
                    flag = input()
                    if flag == 'y' or flag == 'Y':
                        seq_id = line.split()[7]
                    else:
                        seq_id = flag
                    print(seq_id)
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

    return d_point


(s_dict_seq, s_dict_len, s_dict_pos) = read_fasta_to_dic(args.sample)
(r_dict_seq, r_dict_len, r_dict_pos) = read_fasta_to_dic(args.ref)
print(s_dict_pos)
print(r_dict_pos)
d_point = read_file_to_dic(args.infile, s_dict_pos)


print(d_point)

# print('\n')
end_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
print('End Time : {}'.format(end_time))
print('Already Run {}s'.format(time.time()-begin_time))
###########
