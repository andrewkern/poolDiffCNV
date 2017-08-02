# Overview
This repository contains a few Python scripts that can be used to detect copy number variants (CNVs)—large duplication (tandem duplications only, sorry!) or deletion polymorphisms—that differ substantially in allele frequency between two different pooled population genomic samples. These scripts were used in order to find CNV candidates for spatially varying selection in the following papers:

Schrider, DR, Hahn, MW, and Begun, DJ. (2016) Parallel Evolution of Copy-Number Variation Across Continents in Drosophila melanogaster. Molecular Biology and Evolution. 33: 1308-1316.

Schrider, DR, Begun, DJ, and Hahn MW (2013) Detecting highly differentiated copy number variants from pooled population sequencing. Pacific Symposium on Biocomputing. 18:344-355.

Briefly, this uses short-read pool-seq data to detect CNVs that appear to be differentiated in allele frequency between separate pool-seq samples (e.g. from different locales, environments). This is done using the following five steps:

0. Map data from each pooled sample to a reference genome, generating paired-end SAM-formatted alignments using your preferred short read aligner (probably bwa-mem).
1. In each population sample, find discordant read pairs suggestive of a CNV (i.e. read pairs mapped too far apart in the case of a deletion, or "everted" read pairs in the case of a tandem duplication).
2. Cluster discordant read pairs into candidate CNVs. This step is necessary because many discordant read pairs may support the same CNV.
3. Combine clusters across pooled samples. This step is necessary determine, for each candidate CNV, how much support do we have for the CNV in each population sample.
4. Count the number of properly mapped read pairs within each candidate CNV. This gives us information about the read depth of each CNV, which in the case of deletions will be lower in samples where the deletion is at high frequency, and for duplications will be higher in the sample with higher allele frequency.


Once these steps are completed, you will have the number of paired-ends supporting a CNV from each population sample, and also the read depth in the CNV region in these two samples. Using this information one can then identify candidate CNVs that appear to be significantly differentiated according to both of these measures. This step could be performed in several ways. You may devise your own strategy or see Schrider (2016) for one possible strategy. The resulting set of outliers should contain CNVs differentiated in allele frequency across population samples, which therefore may represent candidates for local adaptation. Below I walk through each step with example command lines (skipping step zero).

Note: this code was all written with for Python version 2.7 and would require some tweak to work if you are using version 3 or above.

# Step 1: In each population sample, find discordant read pairs suggestive of a CNV.

For tandem duplications, we are looking for read pairs mapped in the "-+" orientation, as opposed to the normal "+-" orientation where the reads are "facing" each other. If for you are using some long-mate end reads where the expected orientation is "-+" then this script will need tweaking, but for most users this should not be an issue. Anyway, we find these "everted" read pairs thusly:

cat samFileNameForPoolA | findEvertedInserts.py > everted_inserts_poolA.tsv
cat samFileNameForPoolB | findEvertedInserts.py > everted_inserts_poolB.tsv

where samFileNameForPoolA is the pull path to our SAM file with paired-end mapping information for our first pooled sample, and samFileNameForPoolB is the same for our second sample. These are intermediate files used by the next step so you should never have to look at them if things are working properly. But if you are curious (or things seem to be going wrong), the output is tab-separated, with each line corresponding to an everted read pair. The fields are: read id, chromosome, left read's starting position, left read's end, left read's strand, left read's mapping quality, left read's sequence, left read's quality scores, right read's starting position, right read's end, right read's strand, right read's mapping quality, right read's sequence, right read's quality scores.

For deletions, we need to find read pairs in the proper "+-" orientation but that are too far apart—according to some insert size cutoff. This cutoff must be supplied as an argument. In other words, this program outputs all read pairs mapped in proper orientation and whose rightmost read's end position minus the leftmost read's starting position is greater than the specified cutoff. Usage:

cat sameFileNameForPoolA | findDistantInserts.py insertSizeCutoffA > distant_inserts_poolA.tsv
cat sameFileNameForPoolB | findDistantInserts.py insertSizeCutoffB > distant_inserts_poolB.tsv

where the input files are the same as above, and insertSizeCutoffA and insertSizeCutoffB are our upper insert size thresholds for each pool. Again, this script generates intermediate files, but the output is tab-separated, with each line corresponding to a pair of reads that map too far apart. The fields are: read id, chromosome, left read start, left read end, left read strand, right read start, right read end, right read strand, the word "deletion" (very helpful!), one unused field, another unused field, and the inferred insert size based on their mapping locations.

With respect to selecting the insert size cutoffs, I usually go with the 99th percentile of insert size from the distribution of all properly mapped read pairs, but you could adjust this depending on how stringent you want to be—the tradeoff is between specificity (and number of candidate to sift through in later steps) and sensitivity to smaller deletions. You could parse this information out of the SAM file yourself, or use something like PICARD (https://broadinstitute.github.io/picard/). Your mapping program may also give you some information that could be useful for selecting this threshold (e.g. mean insert size and standard deviation).

# Step 2: Cluster discordant read pairs into candidate CNVs.

Now that we have read pairs indicative of CNVs, we have to combine this evidence together to come up with sets of candidate CNVs that are worth a closer look. Though one could devise a more principled (and hopefully practical) manner of accomplishing this, the scripts in this step perform this using clustering heuristics that seem to work well enough. For clustering duplications:

cat everted_inserts_poolA.tsv | python clusterEvertedInserts.py insertSizeCutoffA > clustered_everted_inserts_poolA.tsv
cat everted_inserts_poolB.tsv | python clusterEvertedInserts.py insertSizeCutoffB > clustered_everted_inserts_poolB.tsv

where everted_inserts_pool*.tsv and insertSizeCutoff* are as described above, and clustered_everted_inserts_pool*.tsv contain our tab-delimited intermediate output form this step. Each line in the output corresponds to one cluster of everted read pairs. The first field is a comma-separated list of the format c,s,e where c is the chromosome arm/scaffold, s is the estimated starting point of the duplication, and e is the estimated endpoint. The second field gives the number of read pairs supporting the putative duplication. The remaining fields are comma-separated lists giving information about the individual read pairs supporting the event.

When clustering discordant read pairs, you may also find the occasional suspicious cluster that either looks like a super-massive duplication and/or has far too many inserts supporting it. Such cases are probably assembly artifacts that should be removed. This is an optional step that you would have to perform yourself, but it should be easy enough to awk, modify the clustering scripts, or write your own filtering script.

Now, for deletions:

cat distant_inserts_poolA.tsv | python clusterDistantInserts.py insertSizeDiffCutoffA minNumSupportingInserts minLenCutoff > clustered_distant_inserts_poolA.tsv
cat distant_inserts_poolB.tsv | python clusterDistantInserts.py insertSizeDiffCutoffB minNumSupportingInserts minLenCutoff > clustered_distant_inserts_poolB.tsv

where distant_inserts_pool*.tsv and insertSizeDiffCutoff* are as described above. minNumSupportingInserts specifies the minimum number of distant read pairs required to retain a cluster of deletion-supporting read pairs, and minLenCutoff specifies the minimum length of deletions. All clusters not meeting these criteria will be filtered out. These filtering steps are helpful because there will be probably be a lot of pairs in the distantPairFile (i.e. about 1% of all properly oriented read pairs if you used the upper 1% tail as your cutoff for findDistantInserts.py). Thus, very short deletion candidates and/or those supported by very few inserts may often be false positives. I used cutoffs of 2 inserts and a length of at least 50 for Schrider et al. (2016), but you may find that different cutoffs are more appropriate for your dataset. Everted read pairs on the other hand are usually far less common so these thresholds may not be necessary and I have omitted them from findEvertedInserts.py.

Again, the output is tab-delimited with the following fields: the coordinates of the candidate deletion (comma-separated), the number of inserts supporting the event, and all remaining fields containing information about the individual read pears supporting the event (the formatting of these read-pair fields differs slightly from that produced by findEvertedInserts.py, but you probably don't need to pay too much attention to this information anyway).

In both of these scripts, I am assuming that duplicate reads have already been removed. I am also not filtering any chromosome arms out here (i.e. if there are any tiny scaffolds or scaffolds of randomly placed sequence, you may wish to remove them).

# Step 3: Combine clusters across pooled samples.

In order to find CNVs differentiated in allele frequency across samples, we have to have a way of comparing our clusters across these two samples, which we do as follows for duplications:

python combinedEvertedClustersAcrossPools.py clustered_everted_inserts_poolA.tsv clustered_everted_inserts_poolB.tsv normConstA normConstB distanceCutoff > clustered_everted_inserts_merged.tsv

where the .tsv files are our output from the previous step, normConst1 and normConst2 specify corrections for coverage for the respective samples (discussed below), and distanceCutoff specifies how much the coordinates between events from the two samples are allowed to differ for them to be considered the same event. The usage for is the same. For deletions, we have a different script with the same usage:

python combineDistantClustersAcrossPools.py clustered_distant_inserts_poolA.tsv clustered_distant_inserts_poolB.tsv normConstA normConstB distanceCutoff > clustered_distant_inserts_merged.tsv

If coverage differs across pools, then normConst1 and normConst2 can be use to correct for this. After counting the number of read pairs supporting a mutation (which will be proportional to its allele frequency), the scripts multiply this value by the appropriate normConst to make for a more fair comparison. If dA and dB are the total number of read pairs mapped in proper orientation in your two pooled samples, respectively, normConstA should be ((dA+dB)/2)/dA, and normConstB should be ((dA+dB)/2)/dB. If for whatever reason you do not wish to correct for differences in coverage, you can set these both to 1.

Output: Both scripts output a tab-separated line for each cluster with the following fields:

1) The coordinates of the cluster, comma-separated.
2) The coordinates of the cluster found in pool 1 ("NA" if the cluster was only found in pool 2)
3) The (corrected) number of read pairs supporting the event in pool 1 (0 if only found in pool 2)
4) The coordinates found in pool 2 ("NA" if the cluster was only found in pool 1)
5) The (corrected) number of read pairs supporting the event in pool 2 (0 if only found in pool 1)

These scripts are our way of making a guess as to whether any two putative duplications/deletions found in two different samples are in fact the same mutation, and this guess and one based on some arbitrary rules. So it would be prudent to spend some time going through the output to make sure the results appear to make sense.

# Step 4: Count the number of properly mapped read pairs within each candidate CNV.

The last bit of information we need to find differentiated CNVs is the read depth within each candidate event in each pooled sample. True discordant CNVs will have disparities across pools in both the number of discordantly mapped read pairs supporting the breakpoints (always higher in the population with higher allele frequency), and the number of number of properly mapped read pairs within the region (will differ significantly across pools within differentiated CNVs). Moreover, these disparities will be in the same "direction" (e.g. a duplication with higher frequency in pool A will have more everted read pairs and higher read depth; a deletion with higher frequency in pool A will have more distant read pairs and lower read depth).

We count the number of reads within each region using the following commands for duplications:

python countReadPairsInCNV.py clustered_everted_inserts_merged.tsv samFileNameForPoolA maskedRegionsFileName > clustered_everted_inserts_merged_counts_Aonly.tsv
python countReadPairsInCNV.py clustered_everted_inserts_merged_counts_Aonly.tsv samFileNameForPoolB maskedRegionsFileName > clustered_everted_inserts_merged_counts_both.tsv

and for deletions:

python countReadPairsInCNV.py clustered_distant_inserts_merged.tsv samFileNameForPoolA maskedRegionsFileName > clustered_distant_inserts_merged_counts_Aonly.tsv
python countReadPairsInCNV.py clustered_distant_inserts_merged_counts_Aonly.tsv samFileNameForPoolB maskedRegionsFileName > clustered_distant_inserts_merged_counts_both.tsv

where clustered_everted_inserts_merged.tsv, clustered_distant_inserts_merged.tsv, samFileNameForPoolA, and samFileNameForPoolB are as described above. maskedRegionsFileName is a tab-delimited file listing the coordinates of masked regions (chrom+"\t"+start+"\t"+end) you want to ignore when counting read depth (e.g. repeats). If you do not have this information or for whatever reason do not wish to use it, supply a blank file. The output format is the same as the input file but with two additional fields appended to the end: the number of reads within the CNV in the .sam file (ignoring reads whose starting position is masked), and the total number of unmasked positions within the CNV. So for both dups and dels you will have to run this script twice: once for each pool's corresponding .sam file. You can run it for one pool, and then again for the second pool but using output from the first. You will end up with your original input plus four fields (depth in sample A, unmasked length, depth in sample B, unmasked length again).

# Beyond

At this point you want to identify candidates where both the number of discordant read pairs and the read depth within the event support differentiation in allele frequency between the two samples. One possible strategy would be to take the 5% tails from the empirical distributions of each of these, and filter out all but the candidates that exceed both of these thresholds. The empirical distribution of differences in (or ratios of) the numbers of discordant read pairs could be parsed from the output of step 3. For read depth, I favor a threshold that depends on the length of the event, for example, a threshold for many length bins (i.e. for each length bin, take the read depth ratio of random genomic regions of this length, and identify the tails of these distributions). Schrider et al. (2016) did something very similar to this.
