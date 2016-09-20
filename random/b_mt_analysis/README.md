## mtDNA analysis of Buenrostro data ##

### Sample information ###

| Sample Id|Donor Id| mt-DNA Mutation | Cell Type| No. of Variants Calls | No. of Het Calls |
|:--------:|:------:|:---------------:|:--------:|:---------------------:|:----------------:|
|	SRR2920576 |SU484|	MT-ND5	|	AML, pHSC| 31 | 4 |
|	SRR2920583 |SU575|	MT-ND1  |	AML, LSC | 20 | 1 |
|	SRR2920584 |SU575|	MT-ND1	|	AML, pHSC| 20 | 1 |
|	SRR2920586 |SU583|	MT-CO1	|	AML, LSC | 75 | 3 |
|	SRR2920587 |SU583|	MT-CO1	|	AML, pHSC| 75 | 4 |
|	SRR2920588 |SU583|	MT-CO1  |	AML, pHSC| 75 | 3 |
|	SRR2920592 |SU623|	MTCYB	  |	AML, pHSC| 34 | 2 |
|	SRR2920594 |SU654|	MT-ATP6	|	AML, LSC | 43 | 3 |
|	SRR2920595 |SU654|	MT-ATP6	|	AML, pHSC| 43 | 2 | 


### Files ###
- `var.vcf`
  - raw variant files
- `results.vcf` / `results.tsv`
  - variants filtered with `QUAL > 20`
- `[Sample ID].tsv`
  - sample wise variant calls
-  `[Sample ID]_func_impact.tsv`
  - functional annotation of variants using [MITO MAP](http://www.mitomap.org/MITOMAP)
- `[Sample ID].het.csv`
  - Het Calls for positions where minor allele frequency is above 0.1
