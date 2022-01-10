# Figure 2, Gruhl et al. 2021

## Figure 2A, Figure 2-Figure supplement 1A
### Scripts
- ```mcl.sh``` and ```mcl.rb``` + manual production with Illustrator
- detailed description on how numbers were retrieved can be found in the <b>Materials and methods</b> section of the paper (subsection <i>Identification of shared circRNA loci between species</i>), in short (for clarification):
 - parental genes and circRNAs as defined in the supplied GTF-files (e.g. ```mmCircRNAs.gtf```) were lifted from species A to species B;
 - lifted circRNA loci from species A in species B were then intersected (```bedtools intersect```, minimal overlap option <i>-f</i> set to 0.9) with the annotated circRNA locus from species B to identify shared circRNA loci;
 - using several helper scripts, all loci which were liftable and overlapped with a circRNA in species A and B were transformed into the abc-format for the MCL analysis (see subsection <i>Clustering of circRNA loci between species</i> in <b>Materials and methods</b>);
 - all *.abc-files were merged into one single file on which MCL was run (see ```mcl.sh```), the corresponding data files are now added to the provided data: ```parent_correct.abc``` (level 1: parental gene), ```orthologs.abc``` (level 2: circRNA locus) and ```exact.abc``` (level 3: start/stop exon);
 - using the custom-script ```mcl.rb```, the MCL output was summarised to receive the different "orthology clusters";
 - numbers shown in the figure are taken from the output files ```parent_correct.txt``` (level 1: parental gene), ```orthologs.txt``` (level 2: circRNA locus) and ```exact.txt``` (level 3: start/stop exon). 

### Files ###
```
.
|-- sharedLoci
	|-- exact.abc
	|-- exact.txt
	|-- orthologs.abc
	|-- orthologs.txt
	|-- parent_correct.abc
	|-- parent_correct.txt
``` 
 
## Figure 2B
### Scripts
```phastcons.r```

### Files
```
.
|-- mouse
	|-- phastcons
		|-- phastcons.bed
|-- rat
	|-- phastcons
		|-- phastcons.bed
|-- human
	|-- phastcons
		|-- phastcons.bed
```

## Figure 2C
### Scripts
```constitutiveExons.r```

### Files
```
.
|-- opossum
	|-- exonCounts
		|-- exonCounts.txt
		|-- mapped_exons.txt
|-- mouse
	|-- exonCounts
		|-- exonCounts.txt
		|-- mapped_exons.txt
|-- rat
	|-- exonCounts
		|-- exonCounts.txt
		|-- mapped_exons.txt
|-- rhesus
	|-- exonCounts
		|-- exonCounts.txt
		|-- mapped_exons.txt
|-- human
	|-- exonCounts
		|-- exonCounts.txt
		|-- mapped_exons.txt
```

## Figure 2D
### Scripts
```amplitude.r```

### Files
```
.
|-- opossum
	|-- introns
		|-- codingIntrons_gc.txt
	|-- parentalGenes
		|-- parentSummary.txt
		|-- exontable_gc.txt
|-- mouse
	|-- introns
		|-- codingIntrons_gc.txt
	|-- parentalGenes
		|-- parentSummary.txt
		|-- exontable_gc.txt
|-- rat
	|-- introns
		|-- codingIntrons_gc.txt
	|-- parentalGenes
		|-- parentSummary.txt
		|-- exontable_gc.txt
|-- rhesus
	|-- introns
		|-- codingIntrons_gc.txt
	|-- parentalGenes
		|-- parentSummary.txt
		|-- exontable_gc.txt
|-- human
	|-- introns
		|-- codingIntrons_gc.txt
	|-- parentalGenes
		|-- parentSummary.txt
		|-- exontable_gc.txt

```
## Figure 2E
### Scripts
- manual production with Illustrator