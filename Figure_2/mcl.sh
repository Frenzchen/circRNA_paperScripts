#!/bin/bash

# mcl analysis

# modules
module add SequenceAnalysis/StructurePrediction/mcl/14.137

# variables
base="parent_correct.abc"
name="parent_correct"
local=/scratch/local/monthly/fgruhl/orthology
scripts=~/bin/conservation/orthologous_circRNAs

# prep
cat opossum/conservation/opossum_$base mouse/conservation/mouse_$base rat/conservation/rat_$base rhesus/conservation/rhesus_$base human/conservation/human_$base >> $local/"$name".abc

cd $local

# mcl
mcxload -abc "$name".abc --stream-mirror -write-tab "$name".tab -o "$name".mci 
mcl "$name".mci -I 2
mcxdump -icl out."$name".mci.I20 -tabr "$name".tab -o "$name".dump

ruby $scripts/mcl.rb "$name".dump "$name".txt