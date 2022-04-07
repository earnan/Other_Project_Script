#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   377-1.py
#         Author:   yujie
#    Description:   377-1.py
#        Version:   1.0
#           Time:   2022/04/06 17:44:42
#  Last Modified:   2022/04/06 17:44:42
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################
import argparse
import os
from icecream import ic
import linecache

parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   筛选共有的ssr')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--indir', metavar='[dir]', help='结果目录', type=str, default='F:/3777/stat', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def read_file(file):  # fa将单文件读取为字典及列表
    dict_ssr_type = {}  # 每个ssr类型为键,值为重复次数及其个数组成的字典
    list_repeats = []  # 重复次数
    with open(file, 'r') as f:
        n = 0
        for line in f:
            n += 1
            if line.startswith('Frequency of identified SSR motifs'):
                n1 = n
            elif line.startswith('Frequency of classified repeat types (considering sequence complementary)'):
                n2 = n

    for i in range(n1+3, n2-1):
        line = linecache.getline(file, i).strip()
        if line.startswith('Repeats'):
            [list_repeats.append(i) for i in line.split('\t')]
        else:
            dict_ssr_type[line.split('\t')[0]] = []  # 'a':[]
            for index, ele in enumerate(line.split('\t')):
                if index < len(line.split('\t'))-1 and ele.isdigit():
                    dict_ssr_type[line.split('\t')[0]].append(
                        int(list_repeats[index]))  # 'TTT': {4: 13, 5: 2}
    return dict_ssr_type


def find(list1, list2):  # 寻找俩列表共有元素
    tmp = []
    for x in list1:
        if x in list2:
            tmp.append(x)
    return tmp


# 批量赋值给字典
createVar = locals()
list_name_dict = []
for file_name in os.listdir(args.indir):
    file = os.path.join(args.indir, file_name)
    name = 'dict_'+os.path.basename(file).split('.')[0]
    list_name_dict.append(name)
    createVar[name] = read_file(file)

# 查找共有的ssr类型
list0 = dict_Camellia_sinensis_L_O_Kuntze_cv_Xillian_1.keys()
n = 0
for file_name in os.listdir(args.indir):
    n += 1
    file = os.path.join(args.indir, file_name)
    name = 'dict_'+os.path.basename(file).split('.')[0]
    list0 = find(list0, createVar[name].keys())
dict_total = {}
for i in list0:
    dict_total[i] = {}

# 统计共有ssr类型在不同物种中情况
for file_name in os.listdir(args.indir):
    n += 1
    file = os.path.join(args.indir, file_name)
    name = 'dict_'+os.path.basename(file).split('.')[0]
    for i in dict_total.keys():
        dict_total[i][os.path.basename(file).split('.')[
            0]] = createVar[name][i]
ic(dict_total)

# 找共有的键
# dict_Camellia_sinensis_L_O_Kuntze_cv_Xillian_1
# dict_KF156839
# dict_KF562708
# dict_KJ806280
# dict_MH019307
# dict_MH394406
# dict_MH394410
# dict_MN086819
# dict_MT773375
# dict_MW148820
# dict_MZ043860
# 再找里面不同的重复次数
