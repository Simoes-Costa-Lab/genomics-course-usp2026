---
title: Practical Session
parent: CUT&RUN
nav_order: 2
---

# Practical Session

In this section we will walk through the initial analysis of CUT&RUN data and perform an integration of the CUT&RUN and ATAC-seq datasets. This will allow you to compare the two types of data together.

As mentioned in the differential accessibility script, the narrowPeak files contain a lot of information about the peaks, in addition to just their genomic coordinates. For this analysis we are specifically interested in the genomic coordinates of unique CUT&RUN peaks. The CUT&RUNs we are analyzing were performed to assess the occupancy of the transcription factor c-Jun, which is a key member of the AP-1 transcription factor complex. c-Jun binds DNA as a homodimer with other members of the AP-1 complex such as those in the Fos family of proteins. We have performed the CUT&RUN for c-Jun in experimental duplicates in Day 0 stem cells. We have also performed a CUT&RUN using an IgG antibody in the same cell type.

**QUESTIONS FOR TEH CLASS**

- What is the purpose of performing a CUT&RUN for IgG?
- What does this control help us determine?
- Are there any other controls that you can think of for CUT&RUN experiments?

The first thing we will do is obtain the narrowPeak files containing information on the peaks that were called by macs2. Next, we will extract from the narrowPeak files the genomic coordinates of unique peaks for both c-Jun experimental replicates as well as for the IgG control replicate. Most of the commands below are being run directly on the terminal and not within R studio. Take note of the packages and functions that we will use to perform this analysis.

```bash
mkdir ~/cutrun #create a dir for the outputs
cd ~/cutrun #this commands changes the working directory

ln -s /home/course/cutrun/pg_o91_con_rep1_cjun_cst.mm39_nodups_peaks.narrowPeak . 
ln -s /home/course/cutrun/pg_o91_con_rep2_cjun_cst.mm39_nodups_peaks.narrowPeak .
ln -s /home/course/cutrun/pg_o91_con_rep1_igg.mm39_nodups_peaks.narrowPeak .
ln -s /home/course/cutrun/HOMER_pg_o91_con_cjun_CONSENSUS .
ln -s /home/course/cutrun/pg_o91_con_rep1_cjun_cst.mm39_nodups.bam.bw .
ln -s /home/course/cutrun/pg_o91_con_rep2_cjun_cst.mm39_nodups.bam.bw .

#the command below extracts the genomic coordinates of the unique peaks from the narrowPeak file (chromosome, start, end) and writes them as a bed file
#a bed file is a file that coantins at least three columns (chromsome, start, end) for any genomic position
#in our case, we use bed files as files that contain information on the genomic coordinates of our peaks
cut -f 1-3 ./pg_o91_con_rep1_cjun_cst.mm39_nodups_peaks.narrowPeak | uniq > ./pg_o91_con_rep1_cjun.bed

#the step below sorts the peaks by their positions
#it improves memory usage by the subsequent bedtools commands that we will use
#it is only really necessary for very large datasets but we'll do it anyway so you can make sure to perform it if you have datasets that contain more than, say, 10000 regions
bedtools sort -i ./pg_o91_con_rep1_cjun.bed > ./pg_o91_con_rep1_cjun.sorted.bed

#the command below is doing the same thing for the second c-Jun experimental replicate
cut -f 1-3 ./pg_o91_con_rep2_cjun_cst.mm39_nodups_peaks.narrowPeak | uniq > ./pg_o91_con_rep2_cjun.bed
bedtools sort -i ./pg_o91_con_rep2_cjun.bed > ./pg_o91_con_rep2_cjun.sorted.bed

#the command below is extracting the peaks for the IgG control CUT&RUN
cut -f 1-3 ./pg_o91_con_rep1_igg.mm39_nodups_peaks.narrowPeak | uniq > ./pg_o91_con_rep1_igg.bed
bedtools sort -i ./pg_o91_con_rep1_igg.bed > ./pg_o91_con_rep1_igg.sorted.bed
```

Lets see how many peaks we have in each condition/sample. To do this we simply have to count the rows (or lines) of our bed files using the wc -l (word count -line) function. We can do this because our bed files simply contain a list of the genomic coordinates of our CUT&RUN peaks.
```bash
wc -l ./pg_o91_con_rep1_cjun.bed #6618
wc -l ./pg_o91_con_rep2_cjun.bed #6518
wc -l ./pg_o91_con_rep1_igg.bed #182
```

Next, we want to subtract the IgG peaks from our c-Jun peaks for each replicate. To do this we will use the bedtools subtract function (https://bedtools.readthedocs.io/en/latest/content/tools/subtract.html).

```bash
#make sure you check out the bedtools subtract page because it provides a nice visual illustration of what we are doing below
bedtools subtract -A -a ./pg_o91_con_rep1_cjun.sorted.bed -b ./pg_o91_con_rep1_igg.sorted.bed > ./pg_o91_con_rep1_cjun_no_igg.bed
bedtools subtract -A -a ./pg_o91_con_rep2_cjun.sorted.bed -b ./pg_o91_con_rep1_igg.sorted.bed > ./pg_o91_con_rep2_cjun_no_igg.bed

#let's check the number of peaks in our c-Jun replicates after subtraction
#remember, due to the nature of the peaks, subtraction using bedtools rarely results in mathematically exact operations with the number of peaks
wc -l ./pg_o91_con_rep1_cjun_no_igg.bed #6470
wc -l ./pg_o91_con_rep1_cjun_no_igg.bed #6365
```

Now we can generate a consensus peakset for our c-Jun CUR&RUN by using bedtools intersect. Since we might have replicate-specific noise in each individual peakset that leads to spurious peak calls, we need to make sure that the consensus peakset that we will use for our analyses only contains reproducible peaks. By reproducible peaks we mean peaks that were *reliably* called in both CUT&RUN experimental replicates. To create this consensus peakset we will use the bedtools intersect function (https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) that identify peaks/regions that are present in both of our replicates.

```bash
#below we are running the bedtools intersect function
bedtools intersect -wa -a pg_o91_con_rep1_cjun_no_igg.bed -b pg_o91_con_rep2_cjun_no_igg.bed | uniq > pg_o91_con_cjun_CONSENSUS.bed

#our consensus peakset contains about 4000 peaks (3956 peaks to be exact) and we retain about 60% of the peaks from each replicate
#this is an expected number of peaks for a transcription factor such as c-Jun which exhibits occupancy at specific regions of the genome
```

Now that we have generated the c-Jun consensus peakset, let us perform motif enrichment analysis and learn how to make a nice plot to display the results. But what is motif enrichment analysis? A you learned in the lectures, transcription factors usually bind DNA at specific sequences. This is how they achieve specificity in gene regulation. These specific sequences are called motifs. The motif for the transcription factor Sox9, for example, is AGAACAATGG. When we have a list of genomic regions (such as the peaks in our various peaksets), we might be interested in understanding if certain motifs (sequences) are over-represented or enriched among them. This could help us determine if these regions might be regulated by specific transcription factors. This is why we perform an analysis called a motif enrichment analysis. The details of the statistical models associated with this type of analysis are beyond the scope of our workshop, however it is important understand the working principle. Essentially, a motif enrichment analysis software takes all the sequences in your peakset and searches for motifs that have been associated with transcription factors in the literature. For each motif, it asks the following question: do I find this specific motif among the sequences in this peakset (for example, 4000 peaks) more often than expected when compared with a matched background set of genomic sequences? If the answer is yes, then the motif is considered statistically enriched in the peakset. This is an oversimplification of what actually happens during motif enrichment analysis, but it hopefully helps convey the basic idea.

Searching for enriched motifs in the peaksets of transcription factor CUT&RUNs is useful for several reasons. First, this analysis serves as an important quality control step. This is because, for a good quality transcription factor CUT&RUN, we expect the most enriched motifs to be associated with the transcription factor that we are analyzing. For example, in our c-Jun CUT&RUN we expect to see a significant enrichment of motifs that are associated with the AP-1 transcription factor complex (of which c-Jun is a member). We might see motifs for c-Jun as well as it's binding partners such as Fos, Fosl, Fra, etc (remember that c-Jun binds DNA in a homodimer with other partners). Second, motif enrichment can tell us about other transcription factors and transcription factor complexes that could be binding DNA alongside ouor transcription factor of interest. This might add mroe information to our understanding of how specific genes are regulated.

The software we will be using to perform our motif enrichment analysis is called "Hypergeometric Optimization of Motif EnRichment" or HOMER for short (like Homer Simpson). http://homer.ucsd.edu/homer/motif/

```bash
#don't run! 
#findMotifsGenome.pl ./pg_o91_con_cjun_CONSENSUS.bed /Data/Fjodor/MRE_PROJECT/mouse_genome/mm39_genome_helena/mm39.fa ./HOMER_pg_o91_con_cjun_CONSENSUS -size given
```

Let's make a nice plot that displays our motif enrichment data.

```r
setwd("~/cutrun")
motifs <- read_delim("./HOMER_pg_o91_con_cjun_CONSENSUS/knownResults.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)

nrow(motifs)

#below we are adding a column from 1 - 1049 (the total number of motifs in our dataset) that ranks the motifs
motifs$rank <- c(1:1049)

ggplot(motifs, mapping = aes(x=rank, y=-log10(`P-value`))) + 
  geom_point(size=4, color = "gray") + 
  theme_classic() + 
  theme(axis.text.y = element_text(size=18), 
        axis.text.x = element_text(size=18), 
        axis.title.x = element_text(size = 20, margin = unit(c(5, 0, 0, 0), "mm")),
        axis.title.y = element_text(size = 20, margin = unit(c(0, 5, 0, 0), "mm"))) +
  geom_point(data=subset(motifs, `q-value (Benjamini)` <= 0.05), aes(x=rank, y=-log10(`P-value`)), color="#B878BF", size=5) + 
  theme(axis.ticks.length=unit(.2, "cm"), 
        axis.text = element_text(family = "Helvetica", size = 14), 
        axis.title = element_text(family = "Helvetica", size=16),
        axis.text.x = element_text(colour = "black", margin = margin(t = 10)),
        axis.text.y = element_text(colour = "black", margin = margin(t = 10)),
        axis.title.x = element_text(margin = margin(t = 15)),
        axis.title.y = element_text(margin = margin(r = 15))) +
  labs(y = bquote(-log[10]~"(p-value)"),
       x = "Motif Rank") +
  geom_text_repel(data = subset(motifs, `Motif Name` == "Fosl2(bZIP)/3T3L1-Fosl2-ChIP-Seq(GSE56872)/Homer"),
                  aes(x=rank, y=-log10(`P-value`), label = "Fosl2 (77.35% of Peaks)"), box.padding = 14) +
  geom_text_repel(data = subset(motifs, `Motif Name` == "Fra2(bZIP)/Striatum-Fra2-ChIP-Seq(GSE43429)/Homer"),
                  aes(x=rank, y=-log10(`P-value`), label = "Fra2 (83.75% of Peaks)"), box.padding = 16) +
  geom_text_repel(data = subset(motifs, `Motif Name` == "Fos(bZIP)/TSC-Fos-ChIP-Seq(GSE110950)/Homer"),
                  aes(x=rank, y=-log10(`P-value`), label = "Fos (87.59% of Peaks)"), box.padding = 8) +
  geom_text_repel(data = subset(motifs, `Motif Name` == "Jun-AP1(bZIP)/K562-cJun-ChIP-Seq(GSE31477)/Homer"),
                  aes(x=rank, y=-log10(`P-value`), label = "Jun-AP1 (69.29% of Peaks)"), box.padding = 12)
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutrun/Rplot01.png" width="500">
</div>

</details>

Let us also make a representative cJun bigwig by averaging the bigwigs associated with the c-Jun experimental replicates.

```bash
#bigwigAverage -b ./pg_o91_con_rep1_cjun_cst.mm39_nodups.bam.bw ./pg_o91_con_rep2_cjun_cst.mm39_nodups.bam.bw -o ./pg_o91_con_cjun_AVERAGE.bw
ln -s /home/course/cutrun/pg_o91_con_cjun_AVERAGE.bw ~/cutrun
```

Let us make some tornado and profile plots to explore the relationship between our c-Jun CUT&RUn and ATAC-seq data genome-wide.

```bash
#computeMatrix reference-point -S ./pg_o91_control_TMM_SCALE_AVERAGE.bw ./pg_o91_day3_osteo_TMM_SCALE_AVERAGE.bw ./pg_o91_day6_osteo_TMM_SCALE_AVERAGE.bw -R ./pg_o91_con_cjun_CONSENSUS.bed -a 500 -b 500 -bs 20 -o pg_o91_con_cjun_CONSENSUS_ATAC.mat.gz --missingDataAsZero --referencePoint center -p10

ln -s /home/course/cutrun/pg_o91_con_cjun_CONSENSUS_ATAC.mat.gz ~/cutrun

cd ~/cutrun

plotHeatmap -m pg_o91_con_cjun_CONSENSUS_ATAC.mat.gz -out pg_o91_cjun_peaks_at_timecourse_atac.png --sortUsingSamples 1 --colorMap BuPu --samplesLabel "Day 0 Control" "Day 3 Osteo" "Day 6 Osteo" --whatToShow 'heatmap and colorbar' --plotTitle "ATAC-seq Signal at c-Jun Peaks" --regionsLabel "c-Jun CUT&RUN Peaks" --plotFileFormat png

plotProfile -m pg_o91_con_cjun_CONSENSUS_ATAC.mat.gz --perGroup -out pg_o91_cjun_peaks_at_timecourse_atac_PROFILE.png --samplesLabel "Day 0 Control" "Day 3 Osteo" "Day 6 Osteo" --plotTitle "ATAC-seq Signal at c-Jun Peaks" --regionsLabel "" --plotFileFormat png
```

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutrun/pg_o91_cjun_peaks_at_timecourse_atac.png" width="300">
</div>

</details>

<details>
  <summary>Click to reveal figure</summary>
  
<div align="center">
<img src="/genomics-course-usp2026/assets/images/cutrun/pg_o91_cjun_peaks_at_timecourse_atac_PROFILE.png" width="500">
</div>

</details>
