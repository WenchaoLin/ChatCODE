'''
计算所给序列的单一变异位点和简约信息位点个数。并计算ATCG的平均含量，AT含量，GC含量
Singleton variable sites（单一变异位点）
Parsimony informative sites（简约信息位点）
'''
import sys
from collections import Counter

def count_variable_sites(sequences):
    positions = []
    for i in range(len(sequences[0])):
        bases = set(sequence[i] for sequence in sequences)
        if len(bases) > 1:
            positions.append(i)

    singletons = 0
    parsimony_informative = 0
    for position in positions:
        bases = [sequence[position] for sequence in sequences]
        base_counts = Counter(bases)
        if base_counts.most_common(1)[0][1] == 1:
            singletons += 1
        if len(base_counts) > 1:
            parsimony_informative += 1

    return singletons, parsimony_informative

def calculate_base_content(sequences):
    base_counts = Counter(''.join(sequences))
    total_bases = sum(base_counts.values())

    A_content = base_counts['A'] / total_bases * 100
    T_content = base_counts['T'] / total_bases * 100
    C_content = base_counts['C'] / total_bases * 100
    G_content = base_counts['G'] / total_bases * 100
    AT_content = (A_content + T_content) / 2
    GC_content = (G_content + C_content) / 2

    return A_content, T_content, C_content, G_content, AT_content, GC_content

def process_fasta_file(fasta_file):
    sequences = []
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        current_seq = ''
        for line in lines:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                current_seq = ''
            else:
                current_seq += line.strip().upper()  # 将序列转为大写
        if current_seq:
            sequences.append(current_seq)

    singletons, parsimony_informative = count_variable_sites(sequences)
    A_content, T_content, C_content, G_content, AT_content, GC_content = calculate_base_content(sequences)

    print(f"Singleton variable sites: {singletons}")
    print(f"Parsimony informative sites: {parsimony_informative}")
    print(f"A content: {A_content:.1f}%")
    print(f"T content: {T_content:.1f}%")
    print(f"C content: {C_content:.1f}%")
    print(f"G content: {G_content:.1f}%")
    print(f"A+T content: {AT_content:.1f}%")
    print(f"G+C content: {GC_content:.1f}%")

# 指定要处理的 FASTA 文件路径
fasta_file = sys.argv[1]
# 处理 FASTA 文件
process_fasta_file(fasta_file)
