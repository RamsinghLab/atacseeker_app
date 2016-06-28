import os
import pandas as pd
import sys
import vcf

def print_alleles(data, sep = "\t"):
    num_alleles = len(data.GT.split("/"))
##    print data.GT
    if  num_alleles == 1 and data.GT is "0":
        print sep.join([str(x) for x in ["-", 0., "-", 0.]]),
    elif num_alleles == 2 or data.GT is  "1":
##        print data.GT, rec.ALT
        print sep.join([str(x) for x in [rec.ALT[0], data.HF, "-", 0.]]),
    elif num_alleles >= 3:
        print sep.join([str(x) for x in [rec.ALT[0], data.HF[0], rec.ALT[1], data.HF[1]]]),

if len(sys.argv) < 4:
    print "Usage: python", sys.argv[0], "vcf_file sample_1 sample_2"
    sys.exit()

reader = vcf.Reader(open(sys.argv[1], 'r'))
s1 = int(sys.argv[2])
s2 = int(sys.argv[3])
sep = "\t"

print >> sys.stderr, "*" * 80
print >> sys.stderr, "samples used:"
print >> sys.stderr, "s1", str(s1), ":", reader.samples[s1]
print >> sys.stderr, "s2", str(s2), ":", reader.samples[s2]
print >> sys.stderr, "*" * 80

df1 = pd.read_table(os.path.join(os.path.dirname(sys.argv[1]), reader.samples[s1] + "-table.txt"), index_col="Position")
df2 = pd.read_table(os.path.join(os.path.dirname(sys.argv[1]), reader.samples[s2] + "-table.txt"), index_col="Position")

print "REF", sep,
print sep.join(["S-" + i  for i in ["ALT1", "VAF1", "ALT2", "VAF2"]]), sep,
print sep.join(["N-" + i  for i in ["ALT1", "VAF1", "ALT2", "VAF2"]]), sep,
print sep.join(["CHROM", "POS", "S-(A,C,G,T)", "N-(A,C,G,T)"])
for rec in reader:
    position = int(rec.POS)
    depth1 = df1.loc[position, "BaseCount(A,C,G,T)"]
    depth2 = df2.loc[position, "BaseCount(A,C,G,T)"]
    print rec.REF, sep,
    print_alleles(rec.samples[s1].data, sep)
    print sep,    
    print_alleles(rec.samples[s2].data, sep)
    print sep, rec.CHROM, sep, position, sep, depth1, sep, depth2
