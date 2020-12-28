"""
Adjustments:
1. Change "ind" to be a list of your samples, in the order of the bam file you used for angsd.
2. Change "mafs_infile" and "counts_infile" to your .mafs.gz and .counts.gz files
3. Anywhere cchr and cpos are defined, change the list positions to whatever the last two column numbers are in your file.
4. In the first while loop, change the chromosome and position variables to whatever the chrom and pos is of the last line in your counts file
"""
import gzip
import random

filename_start = 1
filename_end = 2
currchrom = "null"

ind=["Sub1_unadmix","Tmq3_admixed","Har03a_indicus","Ga5_gaur_outgroup"]

maf_infile = "ABtest_4Cow_baseCts.mafs.gz"
counts_infile = "ABtest_4Cow_baseCts_pos.counts.gz"

for i in range(filename_start, filename_end):
    #print("file num:",str(i))
    outfile = "ABtest_alladmixCow_v3.vcf"
    f = gzip.open(counts_infile, 'rt')
    g = gzip.open(maf_infile, 'rt')
    with open(outfile, 'w') as h:
        h.write("##fileformat=VCFv4.1" + "\n")
        h.write('##INFO=<ID=.,Number=1,Type=Character,Description="description">' + "\n")
        h.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' + "\n")
        vcf_header = ["#CHROM",	"POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"] + ind
        vcf_header_print = "\t".join(vcf_header)
        h.write(vcf_header_print + "\n")
        f.readline()
        countline = f.readline()
        countspline = countline.split()
        cchr,cpos = countspline[16:18]
        g.readline()
        mafline = g.readline()
        mafspline = mafline.split()
        mchr,mpos = mafspline[0:2]
        while cchr != "chrM" or cpos != "15412":
            if cchr == mchr:
                currchrom = cchr
                while int(cpos) < int(mpos):
                    countline = f.readline()
                    countspline = countline.split()
                    cchr,cpos = countspline[16:18]
                while int(mpos) < int(cpos):
                    mafline = g.readline()
                    mafspline = mafline.split()
                    mchr,mpos = mafspline[0:2]
                if cpos == mpos:
                    snp_set = set()
                    ind_snps = []
                    for i in range(0,len(countspline)-2,4):
                        snp = 0
                        A_count = int(countspline[i])
                        C_count = int(countspline[i+1])
                        G_count = int(countspline[i+2])
                        T_count = int(countspline[i+3])
                        SNPs_for_ind = "A" * A_count + "C" * C_count + "G" * G_count + "T" * T_count
                        if len(SNPs_for_ind) == 0:
                            snp = "."
                        else:
                            snp = random.choice(SNPs_for_ind)
                            snp_set.add(snp)
                        ind_snps.append(snp)
                    ref_snp = mafspline[4]
                    snp_list = list(snp_set)
                    #print(snp_list)
                    if ref_snp in snp_list and len(snp_set)==2: 
                        #print("true")
                        if snp_list[0] == ref_snp:
                            alt_snp = snp_list[1]
                        else:
                            alt_snp = snp_list[0]
                        vcf_line = [mchr, mpos, ".", ref_snp, alt_snp, ".", "PASS", ".", "GT"]
                        for snp in ind_snps:
                            if snp == ref_snp:
                                vcf_line.append("0")
                            elif snp == alt_snp:
                                vcf_line.append("1")
                            elif snp == ".":
                                vcf_line.append(".")
                        vcf_line_print = "\t".join(vcf_line)
                        h.write(vcf_line_print + "\n")
                    countline = f.readline()
                    countspline = countline.split()
                    cchr,cpos = countspline[16:18]
                    mafline = g.readline()
                    mafspline = mafline.split()
                    mchr,mpos = mafspline[0:2]
            elif mchr != currchrom:
                while cchr != mchr:
                    countline = f.readline()
                    countspline = countline.split()
                    cchr,cpos = countspline[16:18]
            elif cchr != currchrom:
                while mchr != cchr:
                    countline = f.readline()
                    countspline = countline.split()
                    cchr,cpos = countspline[16:18]
        print(cchr,cpos)
        if cchr == mchr and cpos == mpos:
            snp_set = set()
            ind_snps = []
            for i in range(0,len(countspline)-2,4):
                snp = 0
                A_count = int(countspline[i])
                C_count = int(countspline[i+1])
                G_count = int(countspline[i+2])
                T_count = int(countspline[i+3])
                SNPs_for_ind = "A" * A_count + "C" * C_count + "G" * G_count + "T" * T_count
                if len(SNPs_for_ind) == 0:
                    snp = "."
                else:
                    snp = random.choice(SNPs_for_ind)
                    snp_set.add(snp)
                ind_snps.append(snp)
            ref_snp = mafspline[4]
            snp_list = list(snp_set)
            if ref_snp in snp_list and len(snp_set)==2: 
                #print("true")
                if snp_list[0] == ref_snp:
                    alt_snp = snp_list[1]
                else:
                    alt_snp = snp_list[0]
                vcf_line = [mchr, mpos, ".", ref_snp, alt_snp, ".", "PASS", ".", "GT"]
                for snp in ind_snps:
                    if snp == ref_snp:
                        vcf_line.append("0")
                    elif snp == alt_snp:
                        vcf_line.append("1")
                    elif snp == ".":
                        vcf_line.append(".")
                vcf_line_print = "\t".join(vcf_line)
                h.write(vcf_line_print + "\n")
                
f.close()
g.close()
h.close()
