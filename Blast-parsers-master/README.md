# Blast-parsers
These are a variety of scripts, written in Python, to parse Blast outputs in XML format

blast_parser.py returns a tab separated table of the top hits in the Blast output

parse_blast_for_seqs.py parses Blast output and acts as a filter, separating hits without gaps and hits that have gaps in wither the reference or query sequence. It returns multiple files (a tab separated table of ungapped hits, a tab separated table of gapped hits, a fasta file of the ungapped query sequences, a fasta file of the gapped query sequences, and a list of query sequences that did not hit to the reference database.
