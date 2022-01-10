#!/usr/bin/env ruby

# ruby 2.0
# 16/03/16
#
# Find different subclusters
#
# ARGV[0]: dump.cluster
# ARGV[1]: output file
#
# Cluster sorted according to age.


unless ARGV.size == 2
	puts 'wrong number of arguments'
	puts 'usage: ruby mcl.rb dump.species.mci.I20 cluster.txt'
	exit
end

# methods
##########################################################################################

def rna_age(species)
	s = species.values
	
	therian = s[0].any? && (s[1..2].all? {|x| x.any?} || s[3..-1].all? {|x| x.any?}) 
	eutherian1 = s[0].empty? && s[1..2].all? {|x| x.any?} && s[3..-1].any? {|x| x.any?}
	eutherian2 = s[0].empty? && s[1..2].any? {|x| x.any?} && s[3..-1].all? {|x| x.any?}	
 	rodents = s[0].empty? && s[1..2].all? {|x| x.any?} && s[3..-1].all? {|x| x.empty?}
 	primates = s[0..2].all? {|x| x.empty?} && s[3..-1].all? {|x| x.any?}

 	if therian
 		age = 'therian'
 	elsif eutherian1 || eutherian2
 		age = 'eutherian'
 	elsif rodents
 		age = 'rodents'
 	elsif primates
 	 	age = 'primates'
 	else
		i = s.index {|x| x.any?}
 		age = species.keys[i].to_s
 	end
end

# run
##########################################################################################

cluster = {}
cluster_id = 0

File.open(ARGV[0], 'r').readlines.each do |line|
	line = line.strip.split("\t")
	md = line.select {|gene| gene.match(/ENSMODG/)}
	mm = line.select {|gene| gene.match(/ENSMUSG/)}
	rn = line.select {|gene| gene.match(/ENSRNOG/)}
	rm = line.select {|gene| gene.match(/ENSMMUG/)}
	hs = line.select {|gene| gene.match(/ENSG/)}

	#md = line.select {|gene| gene.match(/mdCircRNA/)}
	#mm = line.select {|gene| gene.match(/mmCircRNA/)}
	#rn = line.select {|gene| gene.match(/rnCircRNA/)}
	#rm = line.select {|gene| gene.match(/rmCircRNA/)}
	#hs = line.select {|gene| gene.match(/hsCircRNA/)}
	
	id = "orthoCl-#{cluster_id += 1}"
	cluster[id] = {:opossum => md, :mouse => mm, :rat => rn, :rhesus => rm, :human => hs}
end

File.open(ARGV[1], 'w') do |output|
	cluster.each do |id, species|
		next if species.values.all? {|v| v.empty?}
 		rnas = species.values.collect {|x| x.any? ? x.join('|') : '-'}
		age = rna_age(species)
 		cl_length = species.values.collect {|x| x.length}
		output.puts [id, rnas, cl_length, age].join("\t")
	end
end

