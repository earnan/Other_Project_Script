#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   377-1ssr筛选.py
#         Author:   yujie
#    Description:   377-1ssr筛选.py
#        Version:   1.0
#           Time:   2022/04/07 13:30:27
#  Last Modified:   2022/04/07 13:30:27
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
    '-i', '--indir', metavar='[dir]', help='输入目录', type=str, default='F:/3777/stat', required=False)
optional.add_argument(
    '-o', '--outdir', metavar='[dir]', help='输出目录', type=str, default='F:/3777/out', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

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


def read_file(file):
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
            dict_ssr_type[line.split('\t')[0]] = []  # 形如'a':[]
            for index, ele in enumerate(line.split('\t')):
                if index < len(line.split('\t'))-1 and ele.isdigit():
                    dict_ssr_type[line.split('\t')[0]].append(
                        int(list_repeats[index]))  # 形如'TTT': [4,5]   也可改成'TTT': {4: 13, 5: 2}
    return dict_ssr_type


def find(list1, list2):  # 寻找俩列表共有元素
    tmp = []
    for x in list1:
        if x in list2:
            tmp.append(x)
    return tmp


def merge(list1, list2):  # 列表并集
    tmp = []
    tmp = list(set(list1) | set(list2))
    return tmp


# 总#######################################################################
# 批量赋值给字典
createVar = locals()
list_species = []
for file_item in os.listdir(args.indir):
    file = os.path.join(args.indir, file_item)
    access_id = os.path.basename(file).split('.')[0]
    list_species.append(access_id)
    createVar[access_id] = read_file(file)


# 分########################################################################
# 查找共有ssr的最小单元
dict_total_unique = {}  # 各物种独有ssr
list0 = Camellia_sinensis_L_O_Kuntze_cv_Xillian_1.keys()
for access_id in list_species:
    list0 = find(list0, createVar[access_id].keys())
    dict_total_unique[access_id] = {}  # 各物种独有ssr

for access_id in list_species:
    tmp_list = set(createVar[access_id].keys())-set(list0)
    for i in tmp_list:
        dict_total_unique[access_id][i] = createVar[access_id][i]

with open(os.path.join(args.outdir, '1.txt'), 'w') as f:
    list1 = []
    for i in list_species:
        tmp = dict_total_unique[i].keys()
        list1 = merge(list1, tmp)
    list1.sort()
    [f.write('\t'+j) for j in list1]
    f.write('\n')
    for i in list_species:
        f.write(i)
        f.write('\n')


# ic(dict_total_unique)
# TODOdict_total_unique准备写入文件

dict_total_share = {}  # 共有ssr
for i in list0:
    dict_total_share[i] = {}


#总##################################################################################
# 统计共有ssr最小单元在不同物种中的重复次数
for access_id in list_species:
    for i in list0:
        dict_total_share[i][access_id] = createVar[access_id][i]
# ic(dict_total_share)
# TODO
#分################################################################################

# 进一步统计
dict_total_share_same = {}  # 没有多态性的ssr类型及其重复次数
dict_total_share_unique = {}  # 在共有ssr中重复次数不同
for i in list0:  # ssr类型
    dict_total_share_unique[i] = {}
    list1 = dict_total_share[i]['Camellia_sinensis_L_O_Kuntze_cv_Xillian_1']
    for access_id in list_species:  # 物种
        list1 = find(list1, dict_total_share[i][access_id])
    dict_total_share_same[i] = list1

    for access_id in list_species:  # 物种
        tmp = set(dict_total_share[i][access_id])-set(list1)
        dict_total_share_unique[i][access_id] = list(tmp)

# ic(dict_total_share_same)
# ic(dict_total_share_unique)
# TODOdict_total_share_same写入文件
# TODOdict_total_share_unique写入文件
