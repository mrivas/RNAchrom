.. _exampleCase:

Example case
============

Here, we'll use RNAchrom tools to analyze RNA-chromatin data. In their raw form, RNA-chromatin libraries are paired-end sequencing libraries. Each read pair correspond to an RNA-chromatin link, where one mate comes from the RNA aimer and the other from its chromatin target. For brevity, we'll call these as RNA-mate and DNA-mate, respectively. In the same fashion, we'll call the gene sources of the RNA-mates as aimer genes.

Data summary
------------

We generated two RNA-chromatin libraries on mouse E14: TwoStep1024, and TwoStep1219 (see Table 1). For each library, read pairs were divided into two sets: 

* **Aware pairs**, where at least one mate contained the linker sequence. Based on the direction of the linker sequence, each mate could be resolved as RNA-end or DNA-end.

* **Blind pairs**, where linker sequences were absent in both mates. Therefore, there was not direct way to map mates into RNA-end and DNA-ends.

As the RNA-mates were expected to span splice junctions, we used STAR --as splice tolerant aligner :cite:`Dobin2013`-- to map the reads against mm9. For all sets, we allowed up to 2 mismatches, restricted the intron size to less than 10k nt (Figure :num:`#alignment-diagram`). For both libraries most blind pairs could be aligned (>75%). Aware pairs, on the other hand, were less mappable (alignment rates ~49%) mainly because after removing the linker sequence we ended up with shorter mates sequences. 

.. _alignment-diagram:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/alignment_diagram.svg
   :width: 90%
   
   Construction of blind and aware-links. Paired-end fastq files (fastq1,2) were aligned using STAR (Dobin 2013) allowing for up to 2 mismatches. On one hand, reads not containing the linker sequence resulted in aligned (mates’ distance <10k) and chimeric (mates’s distance >600k, or on different chromosomes) alignments, both of them used to build blind-links. On the other hand, reads containing the linker sequence in were splitted into DNA and RNA fragments (determined by the removed linker sequence orientation) and then aligned to build aware-links. 
   

After the alignment, pairs whose both mates were uniquely mapped were considered **RNA-chromatin links**. Aware or bind pairs with one or two mates having multi-hits were discarded.  To avoid PCR-amplification biases, we removed duplicates. This resulted in ~16M links in both libraries, among them were 6k and 4k aware links, respectively (see Table1).

.. csv-table:: Table1: Number of links per library. To compute the % of links overlapping exons, we requiered at least one mate on a pair should overlap an exon.
   :header: "Sample", "# aware pairs (% uniquely aligned)", "# blind pairs (% uniquely aligned)", "# aware links (% overlapping exons)", "blind links (% overlapping exons)"

   TwoSteps1024, "105,657 (48.88)", "29,505,340 (75.01)", "5,996 (9.9)", "15,752,812 (9.52)" 
   TwoSteps1219, "71,981 (48.53)",  "31,765,882 (78.84)",  "4,357 (7.64)", "16,522,999 (7.98)"

**Inter-chromosomal links** were found only among-blind links: 509 and 368 for libraries TwoSteps1024 an TwoSteps1219, respectively.

Interesting cases: RNA-end over DNA-end enriched genes
------------------------------------------------------

We used aware-links to determine the genes interacting with distal genomic regions. For the genomic region of each gene, interesting cases were judged based on their RNA-ends over DNA-ends ratio. The basic idea is that compared to DNA-ends, RNA-ends should be more concentrated on aimer genes as these contain both:

#. links reaching other genes (targets of the current gene) and,
#. links of self-ligation (DNA of the current gene binding its owns nascent mRNA).

Exceptions to this would be genes that are themselves targets of other aimers. To test whether this is a generalized case we computed the distributions of RNA and DNA-ends only using aware-links overlapping coding regions. The histograms of RNA and DNA-ends per gene (Figure :num:`#mates-over-length`) are not only clearly different (negative skewness Figure :num:`#ratios`), but the mean value of RNA-ends per gene is bigger than its equivalent on DNA, meaning that most genes contain more RNA than DNA ends. Thus, there are clear cases where the number of RNA-ends is bigger than the number of DNA-ends. We focused on those as candidates. To stress the RNA over DNA enrichment we computed their ratio per genomic region (Figure :num:`#ratios`) and select the top 25%.

.. _mates-over-length:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/matesOverLength.svg
   :width: 90%
   
   Distribution of DNA and RNA mates per gene. All values were normalized by gene length.   

.. _ratios:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/ratios.svg
   :width: 90%

   Distribution of RNA over DNA reads per gene. All values were normalzied by gene length. Red line correspond to 75% quantile.

This resulted in 2 candidates genes that are shared among replicates (see Table 2).

.. csv-table:: Table 2: Candidate genes shared among replicates. DNA and RNA-mates correspond to the number of supporting mates on the genomic region of the respective gene (including introns).
   :header: "Gene ensemble ID","gene name", "Length (nt)", "Biotype", "# DNA mates", "# RNA mates", "Enrichment"
   
   ENSMUSG00000078942,Naip6, "36,554", protein coding, 0, 1, 0.693
   ENSMUSG00000018906,P4ha2, "31,571", protein coding, 0, 1, 0.693

Inferred aware-links
--------------------

Since the aware-links are a small fraction of the blind-links, and the previous enriched cases are statistically insignificant (low number of supporting links), we resorted to an inference method of aware-links. Among blind-links, we looked for mates spanning known splicing junctions. When one end of a pair meets this criteria it was inferred that it corresponds to the RNA-end. By asscociation, the other mate was considered as the DNA-mate. As example, see GM108000's links at the `WashU genome browser <http://epigenomegateway.wustl.edu/browser/?genome=mm9&session=KNZWb8e2Mq&statusId=178461723>`_. There, it can be seen that this gene contains only 1 aware-link (from library TwoSteps1219), still, many blind-links spanning known splicing junctions (track \*allSJ) support the connection of this gene with its upstream neighbor Gm10801-201. What's more, this inferred links seems to be supported by the remaining blind-links. 

.. To do
   Redo enrichment analysis without removal of duplicates (this were already removed at the sequence level, thus no need to do it at the alignment level).
   Inferred aware-links and use them to call enriched genes.


Bibliography
------------

.. bibliography:: Mendeley.bib
   :style: plain
   
