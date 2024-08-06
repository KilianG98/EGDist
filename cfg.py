'''config file, I defined the two inputs I used'''
bacter=True

if bacter:
    gff_files="/home/kilian/uni/masters_thesis/data/gff/bacter/"
    csv_file="/home/kilian/uni/masters_thesis/data/roary/bacter/xan_id75_cd80_v2_trim/nitrobacter_exclusive_genes_w_smag.csv"
    genomes=("GCF_000012725","GCF_006539545","GCF_000152905","GCF_002028545","GCA_001896955","GCA_001897285","GCF_000013885")
    outdir="/home/kilian/uni/masters_thesis/scripts/github/results/nbacter_id75"
    eggnog="/home/kilian/uni/masters_thesis/data/exclusive_eggnog/bacter/eggnog_exclusive_annot.emapper.annotations"
    super_cont="GCF_000012725_000000000001"
    
else:
    gff_files="/home/kilian/uni/masters_thesis/data/gff/nitrococcus/"
    csv_file="/home/kilian/uni/masters_thesis/data/roary/nitrococcus/nc_id60_cd_v2/nc_exclusive_genes.csv"
    genomes=("GCF_000153205", "GCA_024102725", "SMAGOTU_00516")
    outdir="/home/kilian/uni/masters_thesis/scripts/github/results/ncoccus_60"
    eggnog="/home/kilian/uni/masters_thesis/data/exclusive_eggnog/nc/eggnog_exclusive_annot.emapper.annotations"
    super_cont="GCF_000153205_000000000001"
