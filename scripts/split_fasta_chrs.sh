cd genome_chrs
csplit -s -z ../genome.fa '/>/' '{*}'
for i in xx* ; do \
	n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
	mv "$i" "chr$n.fa" ; \
done
