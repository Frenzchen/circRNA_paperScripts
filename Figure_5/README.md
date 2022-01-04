# Figure 5, Gruhl et al. 2021


## Figure 5A, Figure 5-Figure supplement 1
### Scripts
```dimer_orthology_plot.r```

Script produces one plot per species, which where then combined manually for the respective figure. Top-5 dimer names were added by hand-

### Files
```
.
|-- opossum
	|-- repeats
		|-- dimers_by_sharedness_v2.txt
		|-- repeatGroups.txt
		|-- overestimation_by_age.txt
|-- mouse
	|-- repeats
		|-- dimers_by_sharedness_v2.txt
		|-- repeatGroups.txt
		|-- overestimation_by_age.txt
|-- rat
	|-- repeats
		|-- dimers_by_sharedness_v2.txt
		|-- repeatGroups.txt
		|-- overestimation_by_age.txt
|-- rhesus
	|-- repeats
		|-- dimers_by_sharedness_v2.txt
		|-- repeatGroups.txt
		|-- overestimation_by_age.txt
|-- human
	|-- repeats
		|-- dimers_by_sharedness_v2.txt
		|-- repeatGroups.txt
		|-- overestimation_by_age.txt
```

## Figure 5B
Manual production of figure with Illustrator. UCSC genome browser views for circRNA loci in the Akt3 gene in human, mouse and opossum were downloaded together with all overlapping repeats present in the top-5 dimers of each species. Repeats from the same family were then counted and collapsed.


## Figure 5C, Figure 5-Figure supplement 3
### Scripts
```TopRepeats_milliDiv.r```
```plot_mfe.r```

Both scripts produce one plot (one for milliDiv, one for MFE) per species, which where then combined manually for the respective figures.

### Files
```
.
|-- opossum
	|-- repeats
		|-- dimers_rnacofold_opossum.txt
		|-- RNAcofold_pairwise_opossum.txt
		|-- usedDimers_milliDiv_v2.txt
		|-- usedDimer_repeats_v2.txt
	|-- GLMs_paper
		|-- geneTable_opossum.txt
|-- mouse
	|-- repeats
		|-- dimers_rnacofold_mouse.txt
		|-- RNAcofold_pairwise_mouse.txt
		|-- usedDimers_milliDiv_v2.txt
		|-- usedDimer_repeats_v2.txt
	|-- GLMs_paper
		|-- geneTable_mouse.txt
|-- rat
	|-- repeats
		|-- dimers_rnacofold_rat.txt
		|-- RNAcofold_pairwise_rat.txt
		|-- usedDimers_milliDiv_v2.txt
		|-- usedDimer_repeats_v2.txt
	|-- GLMs_paper
		|-- geneTable_rat.txt
|-- rhesus
	|-- repeats
		|-- dimers_rnacofold_rhesus.txt
		|-- RNAcofold_pairwise_rhesus.txt
		|-- usedDimers_milliDiv_v2.txt
		|-- usedDimer_repeats_v2.txt
	|-- GLMs_paper
		|-- geneTable_rhesus.txt
|-- human
	|-- repeats
		|-- dimers_rnacofold_human.txt
		|-- RNAcofold_pairwise_human.txt
		|-- usedDimers_milliDiv_v2.txt
		|-- usedDimer_repeats_v2.txt
	|-- GLMs_paper
		|-- geneTable_human.txt
```