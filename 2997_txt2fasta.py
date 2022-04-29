#!/usr/bin/python3
# -*- coding : utf-8 -*-
##########################################################
#
#       Filename:   2997_txt2fasta.py
#         Author:   yujie
#    Description:   2997_txt2fasta.py
#        Version:   1.0
#           Time:   2022/04/29 15:15:36
#  Last Modified:   2022/04/29 15:15:36
#        Contact:   hi@arcsona.cn
#        License:   Copyright (C) 2022
#
##########################################################

import argparse
import os
parser = argparse.ArgumentParser(
    add_help=False, usage='\npython3   格式化为fa标准文件,追加修改内容(替换)')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument('-i', '--input',
                      metavar='[file/dir]', help='目录的话末尾不要加\\,linux不加/', type=str, required=False)
#optional.add_argument('-c', '--check', metavar='[适用于linux]', type=bool, help="linux下使用'-c 1'", default='', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()


def formatting(file):  # 绝对路径
    with open(file, 'r+') as f:  # 可读可写模式
        seq_id = ''
        seq_str = ''
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip('\n')
            else:
                seq_str += line.strip('\n')
        f.write(seq_id+'\n'+seq_str+'\n')
    return 0


def get_judge_from_str(abs_path):  # 判断是目录还是单文件,对每个文件进行同样操作
    item_list = []

    if os.path.isfile(abs_path):
        file = abs_path
        formatting(file)

    elif os.path.isdir(abs_path):
        print('目录')
        for item in os.listdir(abs_path):
            item = os.path.join(abs_path, item)  # 确保是绝对路径
            item_list.append(item)
            get_judge_from_str(item)

    return len(item_list)  # 文件个数


if __name__ == '__main__':  # 该块程序用于安全测试,当整个文件作为模块导入使用时,这块不会运行
    path = args.input
    abs_path = os.path.abspath(path)  # 确保是绝对路径
    number = get_judge_from_str(abs_path)
    print(number)
    print('Done')
