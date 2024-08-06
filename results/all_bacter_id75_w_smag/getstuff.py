#!/usr/bin/env python

import os

superclusters = "GCF_000012725_000000000001_gene_clusters.csv"
genomes = ["GCF_000013885", "GCA_001897285", "GCF_002028545", "GCF_000152905", "GCF_000012725", "GCA_001896955", "GCF_006539545"]
dir_path = os.getcwd()
outputs = []

cur_acc = []
cur_starts = []
cur_names = []
contigs=[]
with open(superclusters, 'r') as file:
    file.readline()
    for line in file:
        line_parts = line.split(',')
        acc = int(line_parts[-1])
        start = int(line_parts[1])
        if len(line) > 6:
            name=""
            for i in range(4,len(line_parts)-1):
                name= name + ',' +line_parts[i]
            name=name[1:]
        else:
            name=line_parts[4]
        
        if not cur_acc:
            cur_acc.append(acc)
            cur_starts.append(start)
            cur_names.append(name)
            continue

        if acc - 10 <= cur_acc[-1]:
            cur_acc.append(acc)
            cur_starts.append(start)
            cur_names.append(name)
            continue

        if len(cur_acc) > 4:
            start1 = cur_starts[0] - 2000
            stop = cur_starts[-1] + 3000
            if start1 < 0:
                start1 = 0
                stop=cur_starts[-1] +300
            output = [("GCF_000012725_000000000001", start1, stop)]

            for genome in genomes:
                for file_name in os.listdir(dir_path):
                    if file_name.startswith(genome) and file_name.endswith('.csv'):
                        c_starts = []
                        contig = os.path.splitext(file_name)[0]
                        contig=contig.replace('_gene_clusters', '')
                        if contig not in contigs:
                            contigs.append(contig)
                        if contig == 'GCF_000012725_000000000001':
                            continue
                        with open(os.path.join(dir_path, file_name), 'r') as f:
                            for name1 in cur_names:
                                f.seek(0)
                                for f_line in f:
                                    f_line_parts = f_line.split(',')
                                    if name1 in f_line:
                                        c_starts.append(int(f_line_parts[1]))

                        if c_starts:
                            print(c_starts, contig, cur_names)
                            start2 = min(c_starts) - 2000
                            if start2 < 0:
                                output.append([contig])
                                continue
                            end1 = max(c_starts) + 3000
                            output.append((contig, start2, end1))
            
            outputs.append(output)

        cur_acc=[acc]
        cur_names=[name]
        cur_starts=[start]

with open("get_clinked.txt", 'w') as output_file:
    count=0
    for op in outputs:
        count+=1
        output_string = "clinker"
        output_string2= " -r"

        for item in op:
            if len(item) > 1:
                genome_id = item[0]
                start_pos = item[1]
                end_pos = item[2]
                output_string += f" {genome_id}.gff3"
                output_string2 += f" {genome_id}:{start_pos}-{end_pos}"
            else:
                output_string += f" {item[0]}.gff3"

        print(f'{output_string+output_string2} -p clusters/cluster{count}.html', file=output_file)

with open('contigs.txt','w') as f:
    for c in contigs:
        print(c, file=f)
