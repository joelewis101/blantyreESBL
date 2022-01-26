get comtifgs cith extract_fastas: ls *.txt | xargs -n1 -I% extract_fastas ../all_dassim_contigs.fa %
map with minimap2 using minimal_all
remove secomdary matches with samtools
