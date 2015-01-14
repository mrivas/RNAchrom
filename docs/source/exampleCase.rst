.. _exampleCase:

Example case
============

Data summary
------------

TwoStep1024 and TwoStep1219

Reads Alignment
---------------

We aligned the paired-end reads with STAR :cite:`Dobin2013`, allowing up to 2 mismatches and restricting the intron size (for mates spanning a spliced junction) to less than 10k nt. Reads were divided into two sets: one where at at least one mate of each pair contains the linker sequence, and another where both mates of each pair doesn’t contain the linker sequence (see Figure :num:`#alignment_diagram`). The first set was further split into RNA and DNA mates (based on the direction of the linker sequence; reads were the linker sequence was on the flanking regions were discarded), and paired-end aligned to produce aware-links. The second set was paired-end aligned to produce blind-links. In all cases duplicates were removed. For libraries TwoStep1024 and TwoStep1219, this resulted in 18M links for both libraries, from where 7k and 5k were aware links (see Table1) 

.. _alignment_diagram:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/alignment_diagram.svg
   :width: 90%

   Construction of blind and aware-links. Paired-end fastq files (fastq1,2) were aligned using STAR (Dobin 2013) allowing for up to 2 mismatches. On one hand, reads not containing the linker sequence resulted in aligned (mates’ distance <600k) and chimeric (mates’s distance >600k, or on different chromosomes) alignments, both of them used to build blind-links. On the other hand, reads containing the linker sequence in were splitted into DNA and RNA fragments (determined by the removed linker sequence orientation) and then aligned to build aware-links. 

.. csv-table:: Table1: Number of links per library. In parenthesis the % of links overlapping exons in at least one mate for all links and on the RNA mate of aware-links.
   :header: "Sample", "# all-links", "# aware-links"

   TwoSteps1024, 18\,053\,131 (9.53%), 6\,898 (6.07%)
   TwoSteps1219, 18\,048\,508 (8.30%), 5\,015 (4.62%)

Interesting cases: RNA over DNA enriched genes
----------------------------------------------

We use aware-links to determine candidate cases for aimers, i.e., genes that interact with distal genomic regions. For the genomic region of each gene, interesting cases were judged based on their RNA over DNA-ends ratio. Compared to DNA-mates, RNA-mates should be more concentrated on aimer genes as these contain both:

#. links reaching other genes (targets of the current gene) and,
#. links of self-ligation (DNA of the current gene binding its owns nascent mRNA).

Exceptions to this would be DNA regions that are targets of other aimers. To test whether this is the case we computed the distributions of RNA and DNA ends only using aware-links overlapping coding regions. The histograms of RNA and DNA-ends per gene (Figure :num:`#matesOverLength`) are not only clearly different (negative skewness Figure :num:`#ratios`), but the mean value of RNA-ends per gene is bigger than its equivalent on DNA, meaning that most genes contain more RNA than DNA ends. Thus, there are clear cases where the number of RNA is bigger than the number of DNA-mates. We focus on those for candidates. To stress the RNA over DNA enrichment we computed their ratio per genomic region (Figure :num:`#ratios`) and select the top 25%.

.. _matesOverLength:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/matesOverLength.svg
   :width: 90%
   
   Distribution of DNA and RNA mates per gene. All values were normalized by gene length.   

.. _rations:

.. figure:: https://132.239.135.28/public/RNAchrom/files/exampleCase/ratios.svg
   :width: 90%

   Distribution of RNA over DNA reads per gene. All values were normalzied by gene length. Red line correspond to 75% quantile.

This resulted in 4 candidates protein coding genes that are shared among replicates (see Table2).

.. csv-table:: Table 2: Candidate genes shared among replicates.
   :header: "Gene", "Length (nt)", "Biotype", "# DNA mates", "# RNA mates", "Enrichment"
   
   ENSMUSG00000058537, 3\,408,  protein coding, 0, 1, 0.693
   ENSMUSG00000075014, 755,     protein coding, 0, 4, 1.61
   ENSMUSG00000078942, 36\,554, protein coding, 0, 1, 0.693
   ENSMUSG00000018906, 31\,571, protein coding, 0, 1, 0.693
   ENSMUSG00000064358, 784,     protein coding, 0, 2, 1.1
