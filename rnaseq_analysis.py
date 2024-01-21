import pandas as pd
import numpy as np
import sys
import subprocess
from Bio.Seq import Seq

print(sys.argv)
file_prefix, mode, post = sys.argv[1], sys.argv[2], sys.argv[3]

#the post option refers to performing analysis separately from alignment
if post == "post":
    file_prefix = "oldalignments/"+file_prefix
else:
    file_prefix = "alignments/"+file_prefix

#Remove unnecessary flags for pandas to process SAM files correctly
with open(file_prefix+"Chimericofftarget.sam","r") as source:
    with open(file_prefix+"CleanedChimericofftarget.sam","w") as dest:
        for line in source:
            cleaned = str(line).split("NH:i")[0]
            dest.write(cleaned+'\n')

gtsm_chim_reads_df = pd.read_csv(file_prefix+"CleanedChimericofftarget.sam",sep="\t", header=None)

# Dictionaries and variables containing neccessary information about each target

reporter_cargo = "gttaatcgaattgaactgaa".upper()
endogenous_cargo = "tccggcagcgaaacccctgg".upper()

gene_dict = {
    "Reporter" : {
        "cargo_seq" : reporter_cargo,
        "off_targ_nt" : 1943,
        "cis_chrom" : "ReporterFullCis",
        "cis_sd" : 1398,
        "cis_sa" : 1942,
        "trans_chrom" : "ReporterFullTrans",
        "trans_sd" : 1398,
        "trans_sa" : 1942
    },
    "ITGB1" : {
        "cargo_seq" : endogenous_cargo,
        "off_targ_nt" : 1006,
        "cis_chrom" : "chr10",
        "cis_sd" : 32922736,
        "cis_sa" : 32923584,
        "trans_chrom" : "ITGB1_trans",
        "trans_sd" : 157,
        "trans_sa" : 1005
    },
    "SMARCA4" : {
        "cargo_seq" : endogenous_cargo,
        "off_targ_nt" : 6394,
        "cis_chrom" : "chr19",
        "cis_sd" : 11035133,
        "cis_sa" : 11041306,
        "trans_chrom" : "SMARCA4_trans",
        "trans_sd" : 220,
        "trans_sa" : 6393
    },
    "TFRC" : {
        "cargo_seq" : endogenous_cargo,
        "off_targ_nt" : 1931,
        "cis_chrom" : "chr3",
        "cis_sd" : 196069569,
        "cis_sa" : 196071395,
        "trans_chrom" : "TFRC_trans",
        "trans_sd" : 104,
        "trans_sa" : 1930
    }
}

gene_attributes = gene_dict[mode]

cargo_seq = gene_attributes["cargo_seq"]

on_target_seq = "aatttgaaggggacactctggttaatcgaattgaactgaa".upper()
unprocessed_seq = "ACTCTTCTTTTTTTTCTCAG".upper()

# Important columns of the SAM file
# column 0 is name
# column 3 is position
# column 9 is sequence

off_targ_nt = gene_attributes["off_targ_nt"]

ts_reads = gtsm_chim_reads_df.loc[gtsm_chim_reads_df[3] == off_targ_nt]
r_names = ts_reads.iloc[:,0].values
r_seqs = ts_reads.iloc[:,9].values
r_mate_map = ts_reads.iloc[:,6].values
# print(len(r_names))

def find_cargoi(seq,query,mismatch_tolerance=4):
    cargoi = Seq(seq).find(query)
    if cargoi !=-1:
        return cargoi
    else:
        matches = []
        for i in range(len(seq) - len(query) + 1):
            # Calculate the number of mismatches between the query sequence and the target sequence, starting at index i.
            mismatches = 0
            for j in range(len(query)):
                if query[j] != seq[i + j]:
                    mismatches += 1
            if mismatches <= mismatch_tolerance:
                matches.append((i,mismatches))
        matches.sort(key = lambda x: x[1])
        try:
            return matches[0][1]
        except:
            print(matches)
            return 1

#run samtools to extract length of BAM
total_nonchim_reads = int(subprocess.check_output(["samtools","view","-c",file_prefix+"Aligned.sortedByCoord.out.bam"]))//2 
#total_nonchim_reads = 1000
#Number of (unfiltered) putative off-target trans-splicing reads
total_chim  = len(r_names)
filtered_chim = 0
full_ot_seqs = []
read_dict = {}

#filter putative off-target trans-splicing reads by length of the upstream exon. only save the non-cargo portion of the read
with open(file_prefix+"otseqs.fasta", "w") as f:
    for i in range(len(r_names)):
        name = r_names[i]
        read = r_seqs[i]
        cargo_i = find_cargoi(read,cargo_seq,4)
        ot_seq = read[:cargo_i]
        if len(ot_seq) >= 20 and not name in read_dict:
            f.write(">"+name+"\n")
            f.write(ot_seq+"\n")
            full_ot_seqs.append((name,read))
            filtered_chim+=1
            read_dict[name] = ot_seq

#also save all putative off-target trans-splicing reads
with open(file_prefix+"fullotseqs.fasta","w") as f:
    for (name,read) in full_ot_seqs:
        f.write(">"+name+"\n")
        f.write(read+"\n")

#Open the STAR splice junction database to capture the number of cis reads
sj_df = pd.read_csv(file_prefix+"SJ.out.tab",sep="\t",header=None)

cis_chrom, cis_sd, cis_sa = gene_attributes["cis_chrom"], gene_attributes["cis_sd"], gene_attributes["cis_sa"]
trans_chrom, trans_sd, trans_sa = gene_attributes["trans_chrom"], gene_attributes["trans_sd"], gene_attributes["trans_sa"]

#Extract the number of reads corresponding to the target cis splice junctions
cis_count = sj_df.loc[(sj_df[0] == cis_chrom) & (sj_df[1] == cis_sd) & (sj_df[2] == cis_sa)].loc[:,6].values[0]
trans_count = sj_df.loc[(sj_df[0] == trans_chrom) & (sj_df[1] == trans_sd) & (sj_df[2] == trans_sa)].loc[:,6].values

if len(trans_count) > 0:
	trans_count = trans_count[0]
else:
	trans_count = 0

#Use BLAST to map the upstream exons. Any upstream exons that are not mappable are filtered out
if filtered_chim>0:
    print("Starting BLAST")
    subprocess.run(["blastn", "-db", "hg38_db", "-query",file_prefix+"otseqs.fasta", "-outfmt", "6", "-max_target_seqs", "5","-out", file_prefix+"mappedot.tab"])
    print("Finishing BLAST")
    blasted_df = pd.read_csv(file_prefix+"mappedot.tab", sep="\t",header = None)
    grouped = blasted_df.loc[blasted_df.groupby([0])[10].idxmin()]
    #print(grouped.head())
    chroms = grouped.loc[:,1].values
    starts = grouped.loc[:,8].values
    ends = grouped.loc[:,9].values
    names = grouped.loc[:,0].values
    total_ot = len(names)
    by_cis = {}
    def parse_key(key):
        s_i = key.find("strand")
        d_i = key.find("splice")
        return (key[:s_i],key[s_i+len("strand"):d_i],key[d_i+len("splice"):])
    for i in range(total_ot):
        chrom = chroms[i]
        start = starts[i]
        end = ends[i]

        if start>end:
            strand = 2
            sd = end-1
        else:
            strand = 1
            sd = end+1

        key = str(chrom)+"strand"+str(strand)+"splice"+str(sd)
        if key in by_cis:
            by_cis[key]+=1
        else:
            by_cis[key] = 1
    spec_data = []
    #add up the number of reads corresponding to the same off-target splice junction
    for key, num_off in by_cis.items():
        (chrom, strand, sd) = parse_key(key)
        chrom = str(chrom)
        strand = int(strand)
        sd = int(sd)
        junction_df = sj_df.loc[(sj_df[0] == chrom) & (sj_df[strand] == sd) & (sj_df[3] == strand)]
        #print(junction_df)
        num_cis = np.sum(junction_df.loc[:,6].values)
        eff_ot = num_off/(num_cis+num_off)
        spec_data.append([chrom,strand,sd,num_cis,num_off,eff_ot])
    normalized_spec = pd.DataFrame(data=spec_data,columns=["Chromosome","Strand","Splice Donor","Number of Cis Junctions","Number of Off-targets","Off-target Efficiency"])
    normalized_spec.to_csv(file_prefix+"normalizedspecificity.csv")
else:
    total_ot = 0

efficiency = trans_count/(cis_count+trans_count)
specificity = trans_count/(trans_count + total_ot)

with open(file_prefix+"finalcounts.csv", "w") as f:
    f.write("Total Non-Chimeric Reads:,"+str(total_nonchim_reads)+"\n")
    f.write("Efficiency:,"+str(efficiency)+"\n")
    f.write("Specificity:,"+str(specificity)+"\n")
    f.write("Cis Count:,"+str(cis_count)+"\n")
    f.write("Trans Count:,"+str(trans_count)+"\n")
    f.write("Off-target Count:,"+str(total_ot)+"\n")
    f.write("Chimeric Count:,"+str(total_chim)+"\n")