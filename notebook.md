## Data

Main data to analyze, from Jon Moller email:

Attached are two .zip files containing the phage resistance phenotype for each phage and the presence/absence of significant SNPs, COGs, and k-mers (determined by the GWAS for each phage). The files either contain a quantitative phage resistance phenotype (0 is sensitive while 1 or more is resistant) or a qualitative phage resistance phenotype (0 is sensitive, 1 is semi-sensitive, and 2 is resistant). The row names are strains while the column names are the features or the phenotype (pheno). Files with "STCC" in the name include sequence type (ST) and clonal complex (CC) as additional predictors (columns).


From Tim Read:

[...] a collection of about 250 Staph genomes and experimentally determined growth parameters for 7 bacteriophages.  You can think of these are continuous measurements with low numbers indicating phage resistance.  Like Michelle, he has performed GWAS and identified interesting features and incorporated them into random forest/ gradient boosted tree models.

Files:
- `phage_qual.zip`
- `phage_quant.zip`


Each file (pXXX) corresponds to a different phage. Phage is a virus that infects bacteria. After the experiment, they see how much bacteria is still there:
- close to 0 -> bacteria sensitive to phage
- close to 1 -> bacteria resistant to phage

We have 8 types of phages. Each type of phage has two files: one with "STCC" which contain two extra columns: one for ST=sequence type, and one for CC=clonal complex

The response is the column `pheno`.

In `quant` folder, the response `pheno` is binary:
- 0-0.5 -> 0
- 0.5-1 -> 1

In `qual` folder, the response `pheno` is ordinal:
- 0.1-0.4 -> 0
- 0.5-0.7 -> 1
- 0.7-1 -> 2

Each column corresponds to a 
- kmer (0/1=presence/abscence of this kmer sequence), e.g. `AACGCGATGGTATTTACAAC`
- SNP (0/1=presence/absence) e.g. `X1_6822_T_A`
- group of genes (0/1= presence/absence of that group of genes) e.g. `group_7666`

Each phage files has different columns. Each phage file has the kmers/SNPs/groups of genes that are thought to be associated with that phage.

## Jon random forest results
File: `phageHostRange_predictiveAccuracy.xlsx`

Note on data from Jon's slides: 
"OD is optical density at 600 nm. These are the turbidities of "S", "SS", or "R" strains determined with the manual assay."

"S = Sensitive (OD600 0.1-0.4), SS = Semi-sensitive (OD600 0.4-0.7), R = Resistant (OD600 0.7+). Essentially a partition of the quantitative OD phenotype"

# Preliminary analyses

We will start with `p002ypresabsSTCC_quant.csv`.

```{r}
> dat=read.table("data/phage_quant/p002ypresabsSTCC_quant.csv", header=TRUE, sep=",")
> str(dat)
'data.frame':	255 obs. of  1761 variables
```

