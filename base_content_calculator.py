import os
import sys
import argparse
import json

def calculate_base_content(fasta_file):
    sequences = {}
    with open(fasta_file, 'r') as file:
        lines = file.readlines()
        current_seq = ''
        for line in lines:
            if line.startswith('>'):
                current_seq = line[1:].strip()
                sequences[current_seq] = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
            else:
                sequence = line.strip().upper()  # 将序列转为大写
                for base in sequence:
                    if base in sequences[current_seq]:
                        sequences[current_seq][base] += 1

    return sequences

def print_base_content_table(sequences, output_file, output_type):
    if output_type == 'json':
        with open(output_file, 'w') as file:
            json.dump(sequences, file, indent=4)
    elif output_type == 'tsv':
        with open(output_file, 'w') as file:
            file.write('Sample\tThymine\tCytosine\tAdenine\tGuanine\n')
            for sample, bases in sequences.items():
                thymine = round(bases['T'] / sum(bases.values()) * 100, 1)
                cytosine = round(bases['C'] / sum(bases.values()) * 100, 1)
                adenine = round(bases['A'] / sum(bases.values()) * 100, 1)
                guanine = round(bases['G'] / sum(bases.values()) * 100, 1)
                file.write(f'{sample}\t{thymine}\t{cytosine}\t{adenine}\t{guanine}\n')
    else:
        print("Invalid output type. Please choose 'json' or 'tsv'.")

def main():
    parser = argparse.ArgumentParser(description='Calculate ATCG base content in a FASTA file.')
    parser.add_argument('-i', '--input', help='Input FASTA file path', required=True)
    parser.add_argument('-t', '--type', help='Output type: json or tsv', default='tsv')
    parser.add_argument('-o', '--output', help='Output file path')
    args = parser.parse_args()

    fasta_file = args.input
    output_type = args.type
    output_file = args.output

    sequences = calculate_base_content(fasta_file)

    if output_file:
        print_base_content_table(sequences, output_file, output_type)
        print(f'Results saved to {output_file}.')
    else:
        print_base_content_table(sequences, '', output_type)

if __name__ == '__main__':
    
    main()
