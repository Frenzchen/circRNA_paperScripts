# Table 1, Gruhl et al. 2021

## Table 1, Supplementary Table 6
### Scripts
```glm_parental_opossum.r```
```glm_parental_mouse.r```
```glm_parental_rat.r```
```glm_parental_rhesus.r```
```glm_parental_human.r```

Each script reads the species-specific data table and runs regression models on the data to find the parental gene predictors. Individual output files are created and summarised by the following scripts:

```summaryTables_parental_glm.r``` (tables were then combined into a single one for the paper)


### Files
```
.
|-- opossum
	|-- geneTable_opossum.txt
|-- mouse
	|-- geneTable_mouse.txt
|-- rat
	|-- geneTable_rat.txt
|-- rhesus
	|-- geneTable_rhesus.txt
|-- human
	|-- geneTable_human.txt
```

## Supplementary Table 7
### Scripts
```glm_AgeHostpot.r```

### Files
```
.
|-- opossum
	|-- geneTable_opossum.txt
|-- mouse
	|-- geneTable_mouse.txt
|-- rat
	|-- geneTable_rat.txt
|-- rhesus
	|-- geneTable_rhesus.txt
|-- human
	|-- geneTable_human.txt
```