# Final shortlist (phenotype-driven)

- Master: `results/19_master/master_scored_all.tsv`
- Phenotype: `results/20_phenotype/phenotype_summary.md`
- Gene panels: `results/20_phenotype/gene_panels.tsv`

## Candidates

### 1. GBE1

- model/kind: ar_homo, variant
- key: `chr3-81648954-C-G`
- impact: `MODERATE` | ClinVar: `Uncertain_significance` | gnomAD AF: `not reported / extremely rare`
- scores: tech=4, panel=8, inh=1, noise=0, final=13
- phenotype note: GSD IV (Andersen); hepatic glycogen accumulation; AR; liver involvement
- IGV: [ ] proband genotype  [ ] parents segregation  [ ] coverage/artefacts
- ACMG draft: 
- Phenotype fit (why/why not): 

#### Disease association
GBE1 is associated with Glycogen Storage Disease type IV (Andersen disease), autosomal recessive, characterized by hepatic glycogen accumulation and progressive liver dysfunction.

#### Mechanism
GBE1 encodes glycogen branching enzyme. Loss-of-function leads to abnormal glycogen structure and storage in hepatocytes.

#### Phenotype fit
+ Increased hepatic glycogen directly matches disease mechanism.
+ Liver dysfunction consistent with GSD IV.
+ AR inheritance matches homozygous model.
- No clear neuromuscular involvement reported (if applicable).

#### IGV validation ![Evidence picture](image.png)

Proband: homozygous ALT (47/47 reads, 100% G, DP=47
Father: heterozygous (14 C / 13 G, balanced)
Mother: heterozygous (19 C / 22 G, balanced)
No strand bias or alignment artefacts observed.
Segregation consistent with autosomal recessive inheritance.

#### ACMG draft

PM2 (variant absent from gnomAD)
PP1 (segregation consistent with autosomal recessive inheritance)
PP3 (multiple in silico tools predict deleterious effect)  
PP4 (phenotype highly specific for GSD IV)
![alt text](<Varsome_results.png>)

Likely pathogenic candidate variant in GBE1.

Conclusion:
Homozygous missense variant in GBE1 (p.Arg198Thr) identified.
Variant absent in population databases.
Segregation consistent with AR inheritance.
Phenotype highly specific for Glycogen Storage Disease type IV.
Classified as likely pathogenic candidate.

### 2. DENND4B

- model/kind: comphet, pair
- key: `DENND4B|chr1-153934829-G-GC|chr1-153934248-C-G`
- impact: `HIGH` | ClinVar: `VUS_OR_OTHER` | AF: `.`
- scores: tech=10, panel=0, inh=1, noise=0, final=11
- phenotype note: .
- No established Mendelian liver/glycogen disorder in OMIM.
- Gene function (vesicle trafficking in keratinocytes) does not explain hepatic glycogen accumulation.
- ClinVar evidence weak/conflicting.
Conclusion: Unlikely to explain phenotype.

### 3. FHAD1

- model/kind: comphet, pair
- key: `FHAD1|chr1-15374592-C-T|chr1-15369450-C-T`
- impact: `HIGH` | ClinVar: `NONE` | AF: `.`
- scores: tech=9, panel=0, inh=1, noise=0, final=10
- phenotype note: .
- IGV: [ ] proband genotype  [ ] parents segregation  [ ] coverage/artefacts
- ACMG draft: 
- Phenotype fit (why/why not): 

### 4. KDM4A

- model/kind: comphet, pair
- key: `KDM4A|chr1-43671586-C-A|chr1-43660291-CA-C`
- impact: `HIGH` | ClinVar: `NONE` | AF: `.`
- scores: tech=9, panel=0, inh=1, noise=0, final=10
- phenotype note: .
- IGV: [ ] proband genotype  [ ] parents segregation  [ ] coverage/artefacts
- ACMG draft: 
- Phenotype fit (why/why not): 
