#! /usr/bin/python

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import sys

parsed_out = open(sys.argv[2], "w")
parsed_out.write("Query_name\t" + "Query_length\t" + "Hit_name\t" + "Query_start\t" + "Query_end\t" + "Frame\t" + "Hit_start\t" + "Hit_end\t" + "ID\n")

gapped_hits = open(sys.argv[3], "w")
gapped_hits.write("Query_name\t" + "Query_length\t" + "Hit_name\t" + "Query_start\t" + "Query_end\t" + "Frame\t" + "Hit_start\t" + "Hit_end\t" + "ID\n")

blast_records = NCBIXML.parse(open(sys.argv[1]))
rpoC1_seqs = []

nohits = open(sys.argv[4], "w")
gapped_seqs = []

for record in blast_records:
	hit_counter = 0
	query_name = record.query.split()[0]
	if record.alignments == []:
		nohits.write(query_name + "\n")
	else:
		for desc in record.descriptions:
			if desc.num_alignments < 2:
				for aln in record.alignments:
					hit_name = aln.hit_def
					for hsp in aln.hsps:
						if hit_counter == 0:
							if len(hsp.query) > 200:
							    query_len = hsp.align_length
							    query_start = hsp.query_start
							    query_end = hsp.query_end
							    hit_start = hsp.sbjct_start
							    hit_end = hsp.sbjct_end
							    frame = hsp.frame[1]
							    hit_id = float(hsp.identities) / float(query_len)
							    evalue = hsp.expect
							    my_id = "|".join([query_name, hit_name, str(query_len), str(evalue), str(hit_id)])
							    if not '-' in hsp.query and not '-' in hsp.sbjct:
							        if hsp.frame[1] > 0:
								    rpoC1_seq = Seq(hsp.query)
								elif hsp.frame[1] < 0:
								    rpoC1_seq = Seq(hsp.query).reverse_complement()
								rpoC1_seqs.append(SeqRecord(rpoC1_seq, id = my_id))
								parsed_out.write(query_name + "\t" + str(query_len) + "\t" + hit_name + "\t" + str(query_start) + "\t" + str(query_end) + "\t" + str(frame) + "\t" + str(hit_start) + "\t" + str(hit_end) + "\t" + str(hit_id) + "\n")
							    else:
								if hsp.frame[1] > 0:
								    gapped_seq = Seq(hsp.query)
								elif hsp.frame[1] < 0:
								    gapped_seq = Seq(hsp.query).reverse_complement()
								gapped_seqs.append(SeqRecord(gapped_seq, id = my_id))
								gapped_hits.write(query_name + "\t" + str(query_len) + "\t" + hit_name + "\t" + str(query_start) + "\t" + str(query_end) + "\t" + str(frame) + "\t" + str(hit_start) + "\t" + str(hit_end) + "\t" + str(hit_id) + "\n")
							    hit_counter += 1

SeqIO.write(rpoC1_seqs, sys.argv[5], "fasta")
SeqIO.write(gapped_seqs, sys.argv[6], "fasta")
nohits.close()
parsed_out.close()
gapped_hits.close()

