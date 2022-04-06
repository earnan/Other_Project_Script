#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename: 编程处理模板.py
#         Author: yuj@genepioneer.cn
#    Description: sample
#  Last Modified: 2022-01-12 16:29:29
#
# Copyright (C) 2021xxxx genepioneer Corporation
##########################################################
import argparse
import os
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   筛选交叉物种')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--input', metavar='[dir]', help='运行时在输入目录下', type=str, required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


# readfasta的输入文件为fa格式,物种很长包含基因的名字
# 单文件处理
def readfasta(input_file):  # fa将单文件读取为字典及列表,最终子函数,后续均采用这个
    seq_id = ''
    seq_dict = {}
    species_id_list = []
    seq = []
    for line in input_file:
        if line.startswith('>'):
            seq_id = line.strip('\n')
            species_id_list.append(line.replace("\n", "").replace(">", ""))
            seq_dict[seq_id] = ''
        else:
            seq_dict[seq_id] += line.strip('\n')
    for value in seq_dict.values():
        seq.append(value)
    return seq, species_id_list, seq_dict


# 获取属+种+变种名,取前4个组合为名字
# 单文件处理
def get_new_species_id_list(species_id_list):  # 获取新的物种名变种名
    for species_id in species_id_list:
        content = species_id.split()
        id = content[1]+'_'+content[2]+'_'+content[3]+'_'+content[4]
        # print(id)
        new_species_id_list.append(id)
    return new_species_id_list


def find(list1, list2):  # 寻找列表共有元素
    tmp = []
    for x in list1:
        if x in list2:
            tmp.append(x)
    return tmp


i = 0
createVar = locals()  # 动态生成变量
for item in os.listdir(args.input):
    print(item)
    if os.path.isfile(os.path.abspath(item)):
        i += 1
        # print(item)
        input_file = open(item, 'r')
        (seq, species_id_list, seq_dict) = readfasta(input_file)
        # print('原个数{}'.format(len(species_id_list)))
        new_species_id_list = []  # 清零
        new_species_id_list = get_new_species_id_list(species_id_list)
        # print('新个数{}'.format(len(new_species_id_list)))
        #createVar['species_id_list_' +item.replace('.txt', '').replace('-', '_')] = new_species_id_list
        createVar['species_id_list_' + str(i)] = new_species_id_list
        print(i)
    if i == 4:
        break

alt = createVar['species_id_list_1']  # alt
print(len(alt))
atp = createVar['species_id_list_2']  # atp
print(len(atp))
his = createVar['species_id_list_3']  # his
print(len(his))
its = createVar['species_id_list_4']  # its
print(len(its))


id1 = find(alt, atp)
print('以下物种同时具有alt,atp两条数据{}'.format(id1))

id2 = find(alt, his)
print('以下物种同时具有alt,his两条数据{}'.format(id2))

id3 = find(alt, its)
print('以下物种同时具有alt,its两条数据{}'.format(id3))

id4 = find(atp, his)
print('以下物种同时具有atp,his两条数据{}'.format(id4))

id5 = find(alt, its)
print('以下物种同时具有alt,its两条数据{}'.format(id5))

id6 = find(his, its)
print('以下物种同时具有his,its两条数据{}'.format(id6))


"""
id2 = find(id1, his)
print('以下物种同时具有alt,atp,his三条数据{}'.format(id2))

id3 = find(id2, its)
print('同时具有4项基因数据的菌株为{}'.format(id3))
"""
