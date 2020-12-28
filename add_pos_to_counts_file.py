import gzip
sites_infile = "ABtest_4Cow_baseCts.counts.gz"
pos_infile = "ABtest_4Cow_baseCts.pos.gz"

#sites_infile = "abtest_head.counts.gz"
#pos_infile = "abtest_snps_head.txt"

outfile = "ABtest_4Cow_baseCts_pos.counts"

pf = gzip.open(pos_infile, 'rt')
f = open(outfile, 'w')
with gzip.open(sites_infile, 'rt') as sf:
    next(sf)
    pf.readline()
    header = ["ind0_A","ind0_C","ind0_G","ind0_T","ind1_A","ind1_C","ind1_G","ind1_T","ind2_A","ind2_C","ind2_G","ind2_T","ind3_A","ind3_C","ind3_G","ind3_T","chr","pos"]
    headerline = "\t".join(header)
    f.write(headerline)
    for line in sf:
        posline = pf.readline()
        pspline = posline.split()
        chrom, pos = pspline[0:2]
        dataline = "\n"+line[:-1]+chrom+"\t"+pos
        f.write(dataline)
