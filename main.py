from cfg import *
import gene_density
import graphs
import csv
import pandas as pd
import os

#prints exclusive genes in the representative genome to csv
def get_supercontig(df):
    row = df.loc[super_cont]
    gps=row['Gene Product']
    data={}
    for  i in range (len(gps)):
        data[gps[i]] = {'Accession': row['Accession'][i],'Egg Description': row['Egg Description'][i], 'Egg Name': row['Egg Name'][i],'Gene Size': row['Gene Size'][i], 'Gene Start': row['Gene Start'][i], 'Strand': row['Strand'][i],  'Egg GOs': row['Egg GOs'][i], 'Egg EC': row['Egg EC'][i], 'KEGG KO': row['KEGG KO'][i], 'KEGG Pathway': row['KEGG Pathway'][i], 'KEGG Module': row['KEGG Module'][i], 'KEGG Reaction': row['KEGG Reaction'][i], 'KEGG Reaction Class': row['KEGG Reaction Class'][i], 'Egg BRITE': row['Egg BRITE'][i], 'KEGG TC': row['KEGG TC'][i], 'CAZy': row['CAZy'][i], 'BiGG Reaction': row['BiGG Reaction'][i], 'PFAMs': row['PFAMs'][i]}
    
    supercontig_df = pd.DataFrame.from_dict(data, orient='index')
    supercontig_df.to_csv(f'{outdir}/supercontig_genes.csv')
    

def main():
    #make outdir if not exists
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        
    outname1 = f'{outdir}/ex_genes_per_contig.csv'
    outname2 = f'{outdir}/gpa.csv'
    
    with open(csv_file, 'r') as file:
        csv_reader = csv.reader(file, quoting=csv.QUOTE_ALL)
        
        #read and print all genes, contigs
        df = gene_density.get_df(gff_files, csv_reader, genomes, eggnog, super_cont)
        df.to_csv(outname1, index=True)
        #filter for contigs with 4 or more esclusive genes
        filtered_df = df[df['Gene Count'] >= 4]
        get_supercontig(df)
    
        
        file.seek(0)
        #gene presence/absence
        gpa = gene_density.create_gene_presence_df(csv_reader, genomes)
        gpa.to_csv(outname2, index=True)
        #plot distribution
        graphs.plot_gene_distribution(filtered_df, gpa, outdir)
        graphs.plot_gene_clusters(filtered_df, gpa, outdir)
        #graphs.get_pa_matrix(gpa)

if __name__ == "__main__":
    main()
