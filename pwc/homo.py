import re
from Bio import SeqIO

def find_pattern_positions(fasta_file, output_file):
    """
    查找FASTA文件中的mononucleotide和dinucleotide模式的位置
    
    mononucleotide规则：
    - G或C连续出现7个及以上
    - A或T连续出现10个及以上
    
    dinucleotide规则：
    - 两个不同碱基交替出现10次及以上（如AGAGAGAGAGAGAGAGAGAG）
    
    参数:
        fasta_file: 输入的FASTA文件路径
        output_file: 输出结果文件路径
    """
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq = str(record.seq).upper()
            seq_id = record.id
            
            # 查找连续的G或C (7个及以上)
            for match in re.finditer(r'([GC])\1{6,}', seq):
                start = match.start()  # 转为0-based坐标
                end = match.end()
                pattern = match.group()
                out.write(f"{seq_id}\t{start}\t{end}\t{pattern}\tmononucleotide\t{len(pattern)}\n")
            
            # 查找连续的A或T (10个及以上)
            for match in re.finditer(r'([AT])\1{9,}', seq):
                start = match.start()  # 转为0-based坐标
                end = match.end()
                pattern = match.group()
                out.write(f"{seq_id}\t{start}\t{end}\t{pattern}\tmononucleotide\t{len(pattern)}\n")
            
            # 查找交替碱基模式 (两个不同碱基交替出现10次及以上)
            for match in re.finditer(r'(([ACGT]{2})\2{9,})', seq):
                full_pattern = match.group(1)
                unit = match.group(2)
                start = match.start()  # 转为0-based坐标
                end = match.end()
                # 检查是否为两个不同的碱基
                if unit[0] != unit[1]:
                    out.write(f"{seq_id}\t{start}\t{end}\t{full_pattern}\tdinucleotide\t{len(full_pattern)}\n")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(
            "Usage: python3 homo.py <reference.fasta> <homo_regions.tsv>\n"
            "\n"
            "Scan the reference for homopolymer / dinucleotide-repeat regions:\n"
            "  - G/C runs >= 7 bp, A/T runs >= 10 bp\n"
            "  - Alternating dinucleotide repeats >= 10 units\n"
            "\n"
            "Output columns: chrom  start  end  pattern  type  length (0-based BED coordinates)\n"
            "Use the output with split_variants.py to partition variants into homo.bed / nonhomo.bed.",
            file=sys.stderr,
        )
        sys.exit(1)

    find_pattern_positions(sys.argv[1], sys.argv[2])

