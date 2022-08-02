from Bio import SeqIO
from Bio.Seq import Seq
from icecream import ic
import argparse
import linecache
import os
import re
import time
import sys

from numpy import str_
parser = argparse.ArgumentParser(
    add_help=False, usage='\
\npython3   format2nex.py\n\
step1\n\
step2\n\
V1.0')
optional = parser.add_argument_group('可选项')
required = parser.add_argument_group('必选项')
optional.add_argument(
    '-i', '--infile', metavar='[infile]', help='infile', type=str, default='E:/', required=False)
optional.add_argument(
    '-o', '--outfile', metavar='[outfile]', help='outfile', type=str, default='F:/', required=False)
optional.add_argument('-c1', '--flag1', help='run step 1?默认是,不运行则-c1',
                      action='store_false', required=False)
optional.add_argument('-c2', '--flag2', help='run step 2?默认否,运行则-c2 ',
                      action='store_true', required=False)
optional.add_argument('-h', '--help', action='help', help='[帮助信息]')
args = parser.parse_args()

# ###################################################格式化字符串


def format_fasta(seq, num):
    format_seq = ""
    for index, char in enumerate(seq):
        format_seq += char
        if (index + 1) % num == 0:  # 可以用来换行
            format_seq += "\n"
    return format_seq


# ############################################################生成tmp
# ####################  70一换行
# open('F:/4689/tmp.nex', 'w') as tmp_handle:
with open('F:/4689/21_sets_seq.aln.fas', 'r') as fa_handle:
    # CHARACTERS部分
    tmp_all_seq_dict = {}  # 直接读取文件 用到的字典
    accession_seq_final_dict = {}  # 登录号  :  序列,有括号信息

    previous_tmp_note = ''  # 上一个id
    previous_final_note = ''  # 上一个id
    tmp_all_seq_dict[previous_tmp_note] = ''  # 赋空值
    # accession_seq_final_dict[previous_tmp_note] = ''  # 赋空值
    format_seq = ''  # 定义
    # 其他部分
    int_NTAX = 0  # 序列总个数

    for line in fa_handle:
        if line.startswith('>'):
            '''写入序列'''
            int_NTAX += 1
            format_seq = format_fasta(tmp_all_seq_dict[previous_tmp_note], 70)
            if format_seq != '':  # 由于是先写序列,后写id,直接写入,那么第一行会空出一行
                # 第一次要写的序列为空,正好可以判断为空就不写入
                # tmp_handle.write(format_seq+'\n'+'\n')
                accession_seq_final_dict[previous_final_note] = format_seq+'\n'+'\n'
            '''读取成字典'''
            tmp_note = line.strip().lstrip('>')  # 有括号信息
            previous_tmp_note = tmp_note
            tmp_all_seq_dict[tmp_note] = ''
            previous_final_note = tmp_note  # .split('(')[0]
            accession_seq_final_dict[previous_final_note] = ''

            note = "'{}'{}".format(tmp_note, '   ')
            # tmp_handle.write(note)  # 后面有3空格的id

        else:
            tmp_all_seq_dict[tmp_note] += line.strip().upper()

    format_seq = format_fasta(
        tmp_all_seq_dict[previous_tmp_note], 70)  # 第509个序列  最后一个序列
    # tmp_handle.write(format_seq+'\n')
    # 最后一个序列
    accession_seq_final_dict[previous_final_note] = format_seq+'\n'+'\n'


# ###################################################对字符串进行解析,获得所有序列号
# #################################打开分组文件  解析用函数


def get_all_accession(tmp_flg, group_dict, tmp_count):
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


# ###############################################################################解析分组
with open('F:/4689/group_mapping.txt', 'r', encoding='utf-8') as group_handle:
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

            tmp_count, list_all = get_all_accession(
                group_name, group_dict, tmp_count)  # 解析
            group_accession_dict[group_name] = list_all  # 登录号列表

            group_taxset_dict[group_name] = '{}-{}'.format(
                start+1, start+int(Numbers))  # TaxSet {t5} = 9-10;
            start = start+int(Numbers)


# #######################################################最终nex里需要用的变量

str_seq_len = len(tmp_all_seq_dict[previous_tmp_note])
str_NTAX = str(int_NTAX)
str_TAXLABELS = ''

str_MATRIX = ''
for i in group_accession_dict.keys():  # A B C D E
    acc_list = group_accession_dict[i]  # a组对应的登录号,无括号信息
    for j in accession_seq_final_dict.keys():  # 键全为登录号,有括号
        if j.split('(')[0] in acc_list:  # 表明属于这一组
            str_TAXLABELS += ("'"+j+"'"+'\n')
            str_MATRIX += ("'"+j+"'"+'   '+accession_seq_final_dict[j])
str_TAXLABELS = str_TAXLABELS.strip()

str_sets = ''
for k, v in group_taxset_dict.items():
    str_sets += '   TaxSet {} = {};\n'.format(k, v)

# #############################################################################写入 nex文件
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
".format('21_sets_seq.aln.fas', str_NTAX, str_TAXLABELS, str_seq_len, str_MATRIX.rstrip(), str_sets)

with open('F:/4689/out/21_sets_seq.nex', 'w') as nex_handle:
    nex_handle.write(s)

# ############################################################范例
'''
NEXUS
[File generated by DnaSP Ver. 6.12.03, from file: {}    Aug 1, 2022]

BEGIN TAXA
DIMENSIONS NTAX = {str_NTAX}

TAXLABELS
{str_TAXLABELS};
END;

BEGIN CHARACTERS;
DIMENSIONS NCHAR={str_seq_len};
FORMAT DATATYPE=DNA  MISSING=? GAP=- ;
MATRIX
{str_MATRIX}
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
