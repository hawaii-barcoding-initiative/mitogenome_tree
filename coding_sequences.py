#!/usr/bin/env python3

from Bio import SeqIO
import sys, os

def eprint(*args, **kwargs):
  print(*args, file=sys.stderr, **kwargs)

def main():
  seq_file = sys.argv[1]
  records = SeqIO.parse(seq_file,"genbank");
  print("accession\tspecies\tgene\tdescription\tnucleotides\tstart\tend")
  gene_desc = {}
  for seq in records:
    species = seq.annotations['organism']
    #  nucleotides = str(seq.seq)
    try:
      accession = seq.annotations['accessions'][0]
    except:
      accession = seq.name
    for f in seq.features:
      if f.type == "CDS":
        cds = f.location.extract(seq.seq)
        start = int(f.location.start)
        end = int(f.location.end)
        description = f.qualifiers['product'][0]
        try:
          gene = f.qualifiers['gene'][0]
          if description not in gene_desc:
            gene_desc[description] = set()
          gene_desc[description].add(gene)
        except:
          if description in gene_desc and len(gene_desc[description]) == 1:
            gene = gene_desc[description].pop()
          else:
            gene="NA"
        print(f"{accession}\t{species}\t{gene}\t{description}\t{cds}\t{start}\t{end}")


if __name__ == "__main__":
  main()