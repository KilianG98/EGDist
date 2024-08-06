import re
import subprocess
import pandas as pd

#greps first line matching a pattern from gff file
def grep(pattern, filename):

    command = ["grep", "-m 1",pattern, filename]
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        # Print the output of grep command
        return(result.stdout)
    except subprocess.CalledProcessError:
        print("Pattern not found or error occurred")

#defines colums oi
def get_columns(header, genomes):
    columns=[pos for pos, gene in enumerate(header) if gene in genomes]
    return columns

#gets accesion number and returns info from the gff file
def get_contig(accession, gff_dir):
    gff=gff_dir+'_'.join(accession.split('_')[:2]) + ".gff"
    line=grep(accession,gff)
    
    match = re.search(r'product=(.*)', line)
    gene_info = match.group(1)

    line=line.split()
    #print(line)

    #greps contig information from the header of gff file
    line2=grep(line[0],gff)
    c_len=line2.split()[-1]

    #return contig, gene lentg, gene start, strand, gene product, length of contig
    return line[0], int(line[4]) -int(line[3]), int(line[3]),line[6],gene_info, c_len
        
#function that returns a df containing information about contig sizes, genenumber, orientation...
import pandas as pd

def get_eggnog_inf(access,eggnog):

    grep_command = f"grep '^{access}' {eggnog}"
    awk_command = "awk -F'\t' '{for(i=8; i<=NF; i++) printf $i \"\t\"; printf \"\\n\"}'"

    # Run the grep command and pipe its output to the awk command
    grep_process = subprocess.Popen(grep_command, shell=True, stdout=subprocess.PIPE)
    awk_process = subprocess.Popen(awk_command, shell=True, stdin=grep_process.stdout, stdout=subprocess.PIPE)

    # Capture the output of the awk command
    output, _ = awk_process.communicate()

    # Decode the output from bytes to string
    output_str = output.decode('utf-8')
    
    
    egg_info=output_str.split('\t')

    if egg_info==['']:
        return['-' for i in range(0,14)]
    
    return egg_info
import pandas as pd


def get_df(gff_files, csv_reader, genomes, eggnog,supercont):
    data = {}
    header = next(csv_reader)
    columns_oi = get_columns(header, genomes)
    supercol=get_columns(header,(supercont))[0]
    
    # Define lists to hold sorted data before inserting into DataFrame
    contigs = []
    gene_sizes_list = []
    gene_starts_list = []
    gene_accessions_list=[]
    strands_list = []
    gene_products_list = []
    contig_lengths = []

    #egg: Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
    egg_desc=[]
    egg_name=[]
    egg_GOs=[]
    egg_EC=[]
    kegg_ko=[]
    kegg_pw=[]
    kegg_module=[]
    kegg_reaction=[]
    kegg_rclass=[]
    egg_brite=[]
    kegg_TC=[]
    cazy=[]
    bigg_reaction=[]
    pfams=[]



    count=1
    
    for line in csv_reader:
        gene_info=line[2]
        if gene_info in gene_products_list:
            count+=1
            gene_info+= '_'+str(count)

        for col in columns_oi:
            if line[col] != '':

                contig, gene_size, gene_start, strand, gene_inf, contig_len = get_contig(line[col], gff_files)
                eggnog_inf=get_eggnog_inf(line[supercol],eggnog)
        
                contigs.append(contig)
                gene_sizes_list.append(gene_size)
                gene_starts_list.append(gene_start)
                strands_list.append(strand)
                gene_products_list.append(gene_info)
                contig_lengths.append(contig_len)
                gene_accessions_list.append(int(line[col].split('_')[-1]))

                egg_desc.append(eggnog_inf[0])
                egg_name.append(eggnog_inf[1])
                egg_GOs.append(eggnog_inf[2])
                egg_EC.append(eggnog_inf[3])
                kegg_ko.append(eggnog_inf[4])
                kegg_pw.append(eggnog_inf[5])
                kegg_module.append(eggnog_inf[6])
                kegg_reaction.append(eggnog_inf[7])
                kegg_rclass.append(eggnog_inf[8])
                egg_brite.append(eggnog_inf[9])
                kegg_TC.append(eggnog_inf[10])
                cazy.append(eggnog_inf[11])
                bigg_reaction.append(eggnog_inf[12])
                pfams.append(eggnog_inf[13])

    # Sort the data based on gene start positions
    sorted_data = sorted(zip(contigs, gene_sizes_list, gene_starts_list, strands_list, gene_products_list, contig_lengths,gene_accessions_list, egg_desc, egg_name, egg_GOs,egg_EC,kegg_ko,kegg_pw, kegg_module, kegg_reaction, kegg_rclass, egg_brite, kegg_TC, cazy, bigg_reaction, pfams)
                         , key=lambda x: x[2])

    # Extract sorted data into separate lists
    contigs, gene_sizes_list, gene_starts_list, strands_list, gene_products_list, contig_lengths,gene_accesions, egg_desc, egg_name, egg_GOs,egg_EC,kegg_ko,kegg_pw, kegg_module, kegg_reaction, kegg_rclass, egg_brite, kegg_TC, cazy, bigg_reaction, pfams= zip(*sorted_data)
    
    # Populate the data dictionary
    for contig, gene_size, gene_start, strand, gene_info, contig_len,gene_accession, egg_d, egg_n, egg_GO, egg_E, kegg_k, kegg_p, kegg_m, kegg_r, kegg_rc, egg_b, kegg_T, cazy, bigg_r, pfam in zip(contigs, gene_sizes_list, gene_starts_list, strands_list, gene_products_list, contig_lengths,gene_accesions, egg_desc, egg_name, egg_GOs, egg_EC, kegg_ko, kegg_pw, kegg_module, kegg_reaction, kegg_rclass, egg_brite, kegg_TC, cazy, bigg_reaction, pfams):
        if contig not in data:
            data[contig] = {'Contig Length': contig_len, 'Gene Count': 0, 'Gene Size': [], 'Gene Start': [], 'Strand': [], 'Gene Product': [], 'Accession':[], 'Egg Description': [], 'Egg Name': [], 'Egg GOs': [], 'Egg EC': [], 'KEGG KO': [], 'KEGG Pathway': [], 'KEGG Module': [], 'KEGG Reaction': [], 'KEGG Reaction Class': [], 'Egg BRITE': [], 'KEGG TC': [], 'CAZy': [], 'BiGG Reaction': [], 'PFAMs': []}
        data[contig]['Gene Size'].append(gene_size)
        data[contig]['Gene Start'].append(gene_start)
        data[contig]['Strand'].append(strand)
        data[contig]['Gene Product'].append(gene_info)
        data[contig]['Accession'].append(gene_accession)
        data[contig]['Egg Description'].append(egg_d)
        data[contig]['Egg Name'].append(egg_n)
        data[contig]['Egg GOs'].append(egg_GO)
        data[contig]['Egg EC'].append(egg_E)
        data[contig]['KEGG KO'].append(kegg_k)
        data[contig]['KEGG Pathway'].append(kegg_p)
        data[contig]['KEGG Module'].append(kegg_m)
        data[contig]['KEGG Reaction'].append(kegg_r)
        data[contig]['KEGG Reaction Class'].append(kegg_rc)
        data[contig]['Egg BRITE'].append(egg_b)
        data[contig]['KEGG TC'].append(kegg_T)
        data[contig]['CAZy'].append(cazy)
        data[contig]['BiGG Reaction'].append(bigg_r)
        data[contig]['PFAMs'].append(pfam)
        data[contig]['Gene Count'] += 1


    df = pd.DataFrame(data).transpose().sort_index()
    return df

def create_gene_presence_df(csv_reader, genomes):
    gene_presence_data = {}
    header = next(csv_reader)
    positions = get_columns(header, genomes)

    # Iterate through CSV reader to populate gene_presence_data dictionary
    count = 1
    for line in csv_reader:
        gene_info = line[2]
        if gene_info in gene_presence_data:
            count += 1
            gene_info += '_' + str(count)

        for col in positions:
            if line[col] != '':
                if gene_info not in gene_presence_data:
                    gene_presence_data[gene_info] = [0] * len(genomes)
                gene_presence_data[gene_info][genomes.index(header[col])] = 1

    # Create DataFrame from gene_presence_data dictionary
    gene_presence_df = pd.DataFrame(gene_presence_data, index=genomes).transpose()
    gene_presence_df.columns = list(genomes)


    return gene_presence_df
