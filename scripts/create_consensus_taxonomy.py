#!/usr/bin/env python

# USAGE
# python create_consensus_taxonomy.py X Y Z A
# where X is the taxonomy mapping file for all NR seqs, Y is the representative
# file (i.e. one of the rep_set/ files with the 119 release), Z is the OTU 
# mapping file created from running pick_otus.py, and A is the output 
# consensus mapping file

from sys import argv
from cogent.parse.fasta import MinimalFastaParser

silva_taxa = open(argv[1], "U")

id_to_taxa = {}

for line in silva_taxa:
    curr_id = line.split()[0].strip()
    curr_taxa = " ".join(line.split()[1:]).strip()
    id_to_taxa[curr_id] = curr_taxa

rep_set_fasta = open(argv[2], "U")

ordered_ids = []
for label,seq in MinimalFastaParser(rep_set_fasta):
    ordered_ids.append(label)

ordered_ids = set(ordered_ids)

otu_mapping = open(argv[3], "U")

matched_ids = {}

for line in otu_mapping:
    if len(line.strip()) == 0:
        continue
    curr_line = line.strip().split('\t')
    curr_otu = curr_line[0]
    all_seqs = curr_line[1:]
    for seq in all_seqs:
        if seq in ordered_ids:
            matched_ids[seq] = all_seqs

output_taxa_mapping = open(argv[4], "w")

for id in ordered_ids:
    taxa_data = []
    final_taxa_string = []
    for curr_seq in matched_ids[id]:
        # Used example from http://stackoverflow.com/questions/3844801/check-if-all-elements-in-a-list-are-identical for comparing equality of elements
        # Have to step backwards through the bottom taxa up to the top taxa
        taxa_data.append(id_to_taxa[curr_seq].split(';')[::-1])
    levels = len(taxa_data[0])
    final_taxa_string = []
    for n in range(levels):
        curr_taxa_strings = []
        for tax in taxa_data:
            curr_taxa_strings.append(tax[n])
        match = all(x==curr_taxa_strings[0] for x in curr_taxa_strings)
        if match:
            final_taxa_string.append(curr_taxa_strings[0])
        else:
            final_taxa_string.append("Ambiguous_taxa")
    fixed_taxa_string = ";".join(final_taxa_string[::-1])
    output_taxa_mapping.write("%s\t%s\n" % (id, fixed_taxa_string))





    
