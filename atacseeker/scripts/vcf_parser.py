import vcf
import sys

def print_alleles(data, sep = "\t"):
    num_alleles = len(data.GT.split("/"))
    if  num_alleles == 1:
        print sep.join([str(x) for x in ["-", 0., "-", 0.]]),
    elif num_alleles == 2:
        print sep.join([str(x) for x in [rec.ALT[0], data.HF, "-", 0.]]),
    elif num_alleles == 3:
        print sep.join([str(x) for x in [rec.ALT[0], data.HF[0], rec.ALT[1], data.HF[1]]]),

reader = vcf.Reader(open('VCF_file.vcf', 'r'))
sep = "\t"

if len(sys.argv) < 3:
    print "Usage: python", sys.argv[0], "sample_1 sample_2"
    sys.exit()

s1 = int(sys.argv[1])
s2 = int(sys.argv[2])

print >> sys.stderr, "*" * 80
print >> sys.stderr, "samples used:"
print >> sys.stderr, "s1", str(s1), ":", reader.samples[s1]
print >> sys.stderr, "s2", str(s2), ":", reader.samples[s2]
print >> sys.stderr, "*" * 80

print "REF", sep,
print sep.join(["S-" + i  for i in ["ALT1", "VAF1", "ALT2", "VAF2"]]), sep,
print sep.join(["N-" + i  for i in ["ALT1", "VAF1", "ALT2", "VAF2"]]), sep,
print sep.join(["CHROM", "POS"])
for rec in reader:
    print rec.REF, sep,
    print_alleles(rec.samples[s1].data, sep)
    print sep,    
    print_alleles(rec.samples[s2].data, sep)
    print sep, rec.CHROM, sep, rec.POS
