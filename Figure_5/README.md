# Figure 5, Gruhl et al. 2019

## Figure 5A, Supplementary Figure 13A
### Scripts
```TopDimers_milliDiv.r```

Script produces one plot per species, which where then combined manually for the respective figure.

### Files
```
.
|-- opossum
	|-- repeats
		|-- usedDimers_milliDiv.txt
	|-- GLMs_paper
		|-- geneTable_opossum.txt
|-- mouse
	|-- repeats
		|-- usedDimers_milliDiv.txt
	|-- GLMs_paper
		|-- geneTable_mouse.txt
|-- rat
	|-- repeats
		|-- usedDimers_milliDiv.txt
	|-- GLMs_paper
		|-- geneTable_rat.txt
|-- rhesus
	|-- repeats
		|-- usedDimers_milliDiv.txt
	|-- GLMs_paper
		|-- geneTable_rhesus.txt
|-- human
	|-- repeats
		|-- usedDimers_milliDiv.txt
	|-- GLMs_paper
		|-- geneTable_human.txt
```

## Figure 5B, Supplementary Figure 13B
### Scripts
```dimer_orthology_plot.r```

Script produces one plot per species, which where then combined manually for the respective figure.

### Files
```
.
|-- opossum
	|-- repeats
		|-- dimers_by_sharedness.txt
		| -- repeatGroups.txt
|-- mouse
	|-- repeats
		|-- dimers_by_sharedness.txt
		| -- repeatGroups.txt
|-- rat
	|-- repeats
		|-- dimers_by_sharedness.txt
		| -- repeatGroups.txt
|-- rhesus
	|-- repeats
		|-- dimers_by_sharedness.txt
		| -- repeatGroups.txt
|-- human
	|-- repeats
		|-- dimers_by_sharedness.txt
		| -- repeatGroups.txt
```

### Figure 5C
Manual production of figure with Illustrator. UCSC genome brwoser views for circRNA loci in the Akt3 gene in human, mouse and opossum were downloaded together with all overlapping repeats present in the top-5 dimers of each species. Repeats from the same family were then counted and collapsed.