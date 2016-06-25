import vcf

def print_alleles(data, sep = "\t"):
    num_alleles = len(data.GT.split("/"))
    if  num_alleles == 1:
        print sep.join([str(x) for x in [0., 0, 0., 0, 0.]]),
    elif num_alleles == 2:
        print sep.join([str(x) for x in [data.DP, rec.ALT[0], data.HF, 0, 0.]]),
    elif num_alleles == 3:
        print sep.join([str(x) for x in [data.DP, rec.ALT[0], data.HF[0], rec.ALT[1], data.HF[1]]]),

reader = vcf.Reader(open('VCF_file.vcf.gz', 'r'))
sep = "\t"
print sep.join(["CHROM", "POS", "REF"]), sep,
print sep.join([reader.samples[0] + "-" + i  for i in ["DEPTH", "ALT1", "VAF1", "ALT2", "VAF2"]]), sep,
print sep.join([reader.samples[1] + "-" + i  for i in ["DEPTH", "ALT1", "VAF1", "ALT2", "VAF2"]]),
print
for rec in reader:
    print rec.CHROM, sep, rec.POS, sep, rec.REF, sep,
    print_alleles(rec.samples[0].data, sep)
    print sep,    
    print_alleles(rec.samples[1].data, sep)
    print