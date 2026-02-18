awk -F'\t' 'BEGIN{OFS="\t"}
NR==1{next}
{
  chrom=$1; pos=$2;
  if(chrom !~ /^chr/) chrom="chr"chrom;
  print chrom, pos, pos
}' results/13_gnomad/master_variants.tsv \
| sort -k1,1 -k2,2n > results/13_gnomad/master_variants.regions.tsv

wc -l results/13_gnomad/master_variants.regions.tsv