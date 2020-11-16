#!/usr/bin/env ruby


`bedtools intersect -a codingExons_all.bed -b UTR_exons.bed -s -f 0.9 -wa | sort -k1,1 -k2,2n | bedtools merge -i - -nms -scores collapse > UTRexons_collapsed.txt`

`bedtools intersect -a codingExons_all.bed -b UTR_exons.bed -s -f 0.9 -v | sort -k1,1 -k2,2n | bedtools merge -i - -nms -scores collapse > nonUTRexons_collapsed.txt`

def exon_type(x, exon_type)
	t = nil
	
	if exon_type == 'utr'
		x.include?('c') ? t = 'c' : t = 'u'
	elsif exon_type == 'non-utr'
		x.include?('c') ? t = 'c' : t = 'o'
	else
		'wrong utr type'	
	end

end

def exon_type2(x)
	x.include?('c') ? 'c' : 'o'
end

#
##########################################################################################

File.open('exons_collapsed.bed', 'w') do |output|
	i = 1

	# read UTR exons
	utrs = File.open('UTRexons_collapsed.txt', 'r').readlines.each do |line|
		chr, start, stop, ids, types = line.strip.split("\t")
	
		types = types.split(',').uniq
		next if types.include?('u')
	
		gene_id = ids.split(';').collect {|id| id.split('|')[0]}.uniq.join(';')
		
		type = exon_type(types, 'utr')
		name = [gene_id, type, i].join('|')
		
		output.puts [chr, start, stop, name, 0].join("\t")
		i += 1
	end

	# read non UTR exons
	utrs = File.open('nonUTRexons_collapsed.txt', 'r').readlines.each do |line|
		chr, start, stop, ids, types = line.strip.split("\t")
	
		types = types.split(',').uniq
		next if types.include?('u')
	
		gene_id = ids.split(';').collect {|id| id.split('|')[0]}.uniq.join(';')
		
		type = exon_type(types, 'non-utr')
		name = [gene_id, type, i].join('|')
		
		output.puts [chr, start, stop, name, 0].join("\t")
		i += 1
	end
end

#`~/bin/bigWigAverageOverBed /archive/cig/kaessmann/fgruhl/genomes/mouse/mm10/mm10_GRCm38/mm10.60way.phastCons60wayPlacental.bw exons_collapsed.bed phastcons.bed`

#`bedtools intersect -a exons_collapsed.bed -b ../parentalGenes/constitutiveExons.bed -f 0.9 -wa | cut -f4 | uniq > constitutiveExons_ids.txt`


#`bedtools intersect -a spliceSite_exons.bed -b UTR_exons.bed -s -f 0.9 -v > nonUTRspliced_exons.bed`

# File.open('AccDon_exons.bed', 'w') do |output|
# 	File.open('nonUTRspliced_exons.bed', 'r').readlines.each do |line|
# 		chr, start, stop, id, type, strand = line.strip.split("\t")	
# 	
# 		acc, don = nil, nil
# 	
# 		if strand == '+'
# 			acc = [start.to_i-99, start.to_i+99]
# 			don = [stop.to_i-99, stop.to_i+99]
# 		else
# 			acc = [stop.to_i-99, stop.to_i+99]
# 			don = [start.to_i-99, start.to_i+99]
# 		end
# 
# 		output.puts [chr, acc, "#{id}|acc"].join("\t")
# 		output.puts [chr, don, "#{id}|don"].join("\t")
# 	end
# end

# scores = {}
# genes = {}
# 
# File.open('whole.bed', 'r').readlines.each do |line|
# 	chr, start, stop, score = line.strip.split("\t")
# 	pos = [chr, start.to_i, stop.to_i]
# 	scores[pos] = score
# end
# 
# puts "finished whole bed"
# 
# File.open('AccDon_exons.bed', 'r').readlines.each do |line|
# 	chr, start, stop, id = line.strip.split("\t")
# 	
# 	genes[id] = {} if !genes.has_key?(id)
# 
# 	start.to_i.upto(stop.to_i) do |n| 
# 	pos = [chr, n, n+1]
# 	scores.has_key?(pos) ? genes[id][pos] = scores[pos] : genes[id][pos] = nil
# 	end
# end
# 
# puts "finished scores"
# 
# File.open('test.txt', 'w') do |output|
# 	genes.each do |gene_id, v|
# 		i = 1
# 		v.each { |pos, score| output.puts [pos, gene_id, score, i].join("\t"); i +=1 }
# 	end
# end


