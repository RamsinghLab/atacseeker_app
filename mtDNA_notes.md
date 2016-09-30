# Mitochondrial DNA Analysis #

## mtDNA Reference ##

This [page](http://haplogrep.uibk.ac.at/blog/rcrs-vs-rsrs-vs-hg19) has some information on the different references available for mtDNA. 

## Literature ##

### mtDNA copy number ###
Katrina's talk
- mentioned this [paper](http://onlinelibrary.wiley.com/doi/10.1002/bies.201500082/abstract) related to mitochondria.
- `mtDNA copy number / cell = 2 * (mtDNA average covearge) / (autosomal DNA average coverage ) `

Other papers that might be of interest:
- [2012 - Quantitative assessment of mitochondrial DNA copies from whole genome sequencing](http://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-S7-S5)
- [2014 - Accurate mitochondrial DNA sequencing using off-target reads provides a single test to identify pathogenic point mutations](http://www.nature.com/gim/journal/v16/n12/full/gim201466a.html)
- [2015 - Assessing Mitochondrial DNA Variation and Copy Number in Lymphocytes of ~2,000 Sardinians Using Tailored Sequencing Analysis Tools](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005306)
- [2016 - Mitochondrial DNA copy number variation across human cancers](https://elifesciences.org/content/5/e10769)

### Heteroplasmy Detection ###

- [2010 - Detecting heteroplasmy from NGS data - Li M, Stoneking M](http://www.cell.com/ajhg/abstract/S0002-9297(10)00370-8)
	- pathogenic heteroplasmy levels are around `0.7-0.9`
	- many non-disease related heteroplasmies have been discovered
	- [MIA](https://github.com/mpieva/mapping-iterative-assembler) assembly used, rCRS reference, 76bps reads, `QS < 10` in at most 2 bases, mean coverage of `76x`
	- criteria for calling heteroplasmic positions
		- `QS >= 20`
		- `QS >= 15` at the 5bp flanking either side
		- MAF of `0.1`
		- minor allele observed in atleast __two reads__ in each direction
- [2012 - NGS of mtDNA uncovers high heteroplasmy](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1002737)
	- rCRS is 16,569 bp
		- samll (6.8%) non-coding displacement loop (D-loop) or control region
		- large (93.2%) coding region housing 37 genes 
			- 22 tRNAs, 13 proteins and 2 rRNA that encode proteins critical to electron transport chain
	- higher variability than nuclear genome
	- hypervariability within the D-loop as compared to rest of mtGenome
	- mtDNA haplogroups are markers of individual ancestry
	- 100-1000 fold higher mutation rate in mtDNA as compared to nuclear genome is owing to the lack of DNA repair system
	- 32% of intial reads were eliminated
	- 120-fold level of sequence coverage
	- 71 sites were found to be heteroplasmic across 40 samples of European and African ancestry
	- 17/71 (29%) sites were validated by experiment but FP is not believed by authors
	- Picardi and Pesole show that ~1% of all reads map to the mtgenome and not to known _numts_
- [2012 - mtDNA gleaned from WES]()
	- Higher level of heteroplasmy proposed from 1K genome data than that considered here
	- This discrepancy is due to the lack of [NumtS](https://en.wikipedia.org/wiki/Numt) sequence filtering in the 1K Genomes study
- [2013 - MITObim](http://nar.oxfordjournals.org/content/41/13/e129)
	- other tools [MIA](https://github.com/mpieva/mapping-iterative-assembler) and [IMAGE](https://sourceforge.net/projects/image2)
	- pre-processing NGS data
		- Error-correcting NGS read using [SOAPdenovo2](https://github.com/aquaskyline/SOAPdenovo2)
		- quality trimming using [MIRA](http://sourceforge.net/projects/mira-assembler)
		- [MITOS](http://mitos.bioinf.uni-leipzig.de/index.py) webserver was used for automated annotation of obtained mtDNA genomes
- [2014 - MToolBox](http://www.ncbi.nlm.nih.gov/pubmed/25028726)
	- [RSRS](http://www.cell.com/ajhg/abstract/S0002-9297(12)00146-2): Reconstructed Sapiens Reference Sequence
	- [rCRS](http://www.nature.com/ng/journal/v23/n2/full/ng1099_147.html): revised Cambridge Reference Sequence
	- Outputs `contigs.fa` and finds out haplogroup to investigate private variants that might be of clinical interest
- [2016 - mtDNA heteroplasmy using ChIP-seq data](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0996-y)
	- several identical mtDNA copies (2-10) are present in each individual mitochondrion, which means a single cell can contain hundreds to thousands of copies of mtDNA
	- het detection filters
		- `QS > 23`
		- minimum coverage threshold of 20
		- minimum heteroplasmy level of 0.15
		- average coverage is 60 (SD 25)
	- 79 individuals from 16 speicies
		- 107 positions in 45 individuals across 14 species
		- 57% of individuals express heteroplasmy
		- 44 positions in intergenic regions, 39 in non-protein-coding genes, 24 in protein coding genes
		- in protein coding gene changes, 13 are synonymous and 11 non synonymous varaints
		- in humans, 5% of positions are associated with disease
	- [MITOMAP](http://www.mitomap.org/MITOMAP) to check heteroplasmic positions for disease
	- liver tisue has one of the highest relative number of heteroplasmies compared to other human tissue

I think the pipeline I want to implement is:
- Use BWA Aligner to map reads to hg19RCRS reference
- Pull out chrRSRS reads and convert `bam` to `fastq`
- Use [MITObim](https://github.com/chrishah/MITObim) or [MIA](https://github.com/mpieva/mapping-iterative-assembler) to make a consensus 
- align reads using bwa against the consensus
- use the heteroplasmy python script to detect heteroplasmies

This, I think, should account for indels, use rCRS as reference and employ required filtering criteria for heteroplasmy detection. The thing I am worried about is how do we map the heteroplasmies back to RCRS. Also, we could use the consensus for normal cells to call heteroplasmies in senescent. Also, I might want to see how I can take care of NUMTs.

## ATACseq mtDNA ##

Because in ATACseq a lot of the reads are biased for the mitochondrial DNA, it seems like an obvious thing to use ATACseq reads to call variants on mt-DNA. There are a few caveats though. The mt-DNA is haploid, inherited completely from the mother. However, most callers like `gatk` and `samtools` assume a diploid genome. 

This [post](http://gatkforums.broadinstitute.org/gatk/discussion/1214/can-i-use-gatk-on-non-diploid-organisms) says that post-version 3.3 `gatk` should be able to handle non-diploid genome. But a later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3345/question-about-ploidy-when-mtdna-variants-calling) categorically denied any work by gatk relating to mt-DNA. Also, this later [post](http://gatkforums.broadinstitute.org/gatk/discussion/3651/gatk-for-mitochondrial-dna-mutations) claims that the user might have to play around with the ploidy settings of `gatk` but does not offer any conclusive answers. 

`samtools mpileup` explicitly says that it is for diploid genomes - as seen here:

```shell
Asifs-iMac:~ asifzubair$ samtools mpileup

Usage: samtools mpileup [options] in1.bam [in2.bam [...]]

  ...
  ...
  ...

Notes: Assuming diploid individuals.

```

However, a question that begs to be asked is that - if indeed we go out looking for heteroplasmy, shouldn't a heterozygous position by `samtools` give us evidence of the same. I fear my reasoning might be naive. 

In any case, I tried using samtools to call variants following the [CUREFFI](http://www.cureffi.org/2012/09/07/an-alternative-exome-sequencing-pipeline-using-bowtie2-and-samtools) folks. However, as I am interested in only reads that align with the mt-DNA, I filtered for those first. 

```shell
## Get chrM reads
for file in *bam; do  
    samtools view -b $file chrM > $(basename $file .bam)".chrM.bam"; 
done

## call variants
samtools mpileup -d 200 -D -B -f ../../tmp/reference/chrM.fa -b bamlist.txt -u | \
bcftools view -  > variants.vcf
```

However, when I looked at the VCF all calls were reported as missing. Also, in some cases `samtools` complained that the mate was missing. I then got counsel from this biostar [page](https://www.biostars.org/p/74386/) which says to use the `-r` option in place of filtering, which made my command even simpler.

```shell
samtools mpileup -d 200 -D -B -f ../../tmp/reference/hg19.fa -r chrM -b bamlist.txt -u | \
bcftools view -  > variants.vcf
```

Still, no cigar ! 

Based on a recommendation by [GGD](http://www.gettinggeneticsdone.com/2014/03/visualize-coverage-exome-targeted-ngs-bedtools.html), I checked that `bedtools` can be used for checking [coverage](https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md) and I tried to do the same for one ATACseq sample. I can't remember the exact command, but it was probably something like this.

```shell
# Calculate a histogram of coverage for each chromosome
    # as well as genome-wide.
    # Estimate: 30 minutes
    bedtools genomecov \
             -ibam NA19146.bam \
    > NA19146.coverage.hist.txt
```

I remember that when I ran the above command, per nucleotide coverage in chrM was quite high - above 100 bases. 

__4/13__: I finally got [MToolbox](https://github.com/mitoNGS/MToolBox) to work. The trick was getting the right indexes for `gsnap`. I used `gsnapv7` which is the same as the `MToolbox`'s page. A thing I had difficulty with is that my compiled version of `gsnap` will not work on all nodes on the cluster. 

To get all the mtDNA reads from aligned bams, I did the following:
```shell
samtools view -h atac.bam "chrM" > out.sam
samtools view -Sb out.sam > out.bam
```
[Biostar](https://www.biostars.org/p/56246) has some comments on using the `-f 4` flag with `samtools`. __must__ check with Tim.
