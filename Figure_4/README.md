# Figure 4, Gruhl et al. 2021

## Figure 4A, Figure 4-Figure supplement 1
### Scripts
```TEenrichment.r```

The script produces the basic plots. For the figure itself, plots for mouse and rhesus macaque were vertically reflected and aligned with the rat/human plot. Red numbers representing the top-3 TEs were added by hand, as well as p-values to bottom left/right of each plot (value caculated and printed by script). Plots were than arranged by hand to fit the page width.

### Files
```
.
|-- opossum
	|-- repeats
		|-- backgroundIntrons.bed
		|-- flankingIntrons.bed
		|-- repeat_backgroundIntron.intersection
		|-- repeat_flankingIntron.intersection
|-- mouse
	|-- repeats
		|-- backgroundIntrons.bed
		|-- flankingIntrons.bed
		|-- repeat_backgroundIntron.intersection
		|-- repeat_flankingIntron.intersection
|-- rat
	|-- repeats
		|-- backgroundIntrons.bed
		|-- flankingIntrons.bed
|-- rhesus
	|-- repeats
		|-- backgroundIntrons_lifted.bed
		|-- flankingIntrons_lifted.bed
		|-- repeat_backgroundIntron.intersection
		|-- repeat_flankingIntron.intersection
|-- human
	|-- repeats
		|-- backgroundIntrons.bed
		|-- flankingIntrons.bed
		|-- repeat_backgroundIntron.intersection
		|-- repeat_flankingIntron.intersection
```

## Figure 4B, D and E
### Scripts
```figure4BDE.r```

### Files
```
.
|-- opossum
	|-- repeats
		|-- overestimation_by_gene.txt
|-- mouse
	|-- repeats
		|-- dimers_mouse.txt
		|-- mouse_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_rodents.txt
		|-- overestimation_by_gene.txt
|-- rat
	|-- repeats
		|-- overestimation_by_gene.txt
|-- rhesus
	|-- repeats
		|-- overestimation_by_gene.txt	
|-- human
	|-- repeats
		|-- overestimation_by_gene.txt
```

## Figure 4C, Figure 4-Figure supplement 2
### Scripts
```phylogenetic_tree.r```

The script produces the basic plot. It was manually modified to (1) highlight the top-5 most recent dimers (purple text and box), (2) change linetype (dotted) of branches to seperate TE families from each other and to (3) add addtional information on the TE origin (red and black text).

### Files
```
.
|-- opossum
	|-- repeats
		|-- opossum_distmx.txt
|-- mouse
	|-- repeats
		|-- mouse_distmx.txt
|-- rat
	|-- repeats
		|-- rat_distmx.txt
|-- rhesus
	|-- repeats
		|-- rhesus_distmx.txt
|-- human
	|-- repeats
		|-- human_distmx.txt
```

## Figure 4F, Figure 4-Figure supplement 2
### Scripts
```dimer_age_pca.r```

Dimer labels do sometimes overlap if the PCA points are too close. Labels were separated manually in the figure itself for readability. Cumulative explained variance was printed from within script (see comment line 92) and added manually to the top-left of each plot.

### Files
```
.
|-- opossum
	|-- repeats
		|-- dimers_opossum.txt
		|-- opossum_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_opossum.txt
|-- mouse
	|-- repeats
		|-- dimers_mouse.txt
		|-- mouse_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_rodents.txt
|-- rat
	|-- repeats
		|-- dimers_rat.txt
		|-- rat_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_rodents.txt
|-- rhesus
	|-- repeats
		|-- dimers_rhesus.txt
		|-- rhesus_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_rhesus.txt
|-- human
	|-- repeats
		|-- dimers_human.txt
		|-- human_distmx.txt
		|-- repeat_names.txt
		|-- RNAcofold_pairwise_human.txt
```