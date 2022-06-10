#!/usr/bin/env python3


# this script tries to extract all protein-coding genes and 16S and 12S from a mitogenome in genbank format

from Bio import SeqIO
import sys, os

# list all possible names of 12S and 16S rRNAs (list according to MIDORI folks)
rna_12s = ["12 ribosomal RNA","12 S ribosomal RNA","12S","12S (SSU) ribosomal RNA","12S mitochondrial ribosomal RNA","12S mitochondrial ribosomal RNA","12S ribosmal RNA","12S ribosoma RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA","12S ribosomal RNA gene","12S ribosomal RNA subunit","12S ribosomal RNA, small subunit","12S ribosomar RNA","12S ribosome RNA","12S ribosoml RNA","12S ribosonal RNA","12S ribosormal RNA","12S ribosormal RNA","12S ribsomal RNA","12S ribsomal RNA","12S ribsosmal RNA","12S RNA","12S small ribosomal RNA","12S small ribosomal RNA subunit","12S small subunit ribosomal RNA","12S small subunit ribosomal RNA","12S small subunit ribosomal RNA","12S small subunit ribosomal RNA","12S small subunit RNA","12S-rRNA","12SrRNA","12SS ribosomal RNA","2S ribosomal RNA","mitochondrial SSU ribosomal RNA","ribosomal RNA small subunit","ribosomal RNA small subunit","rns","rRNA small subunit","rRNA small subunit","rRNA-12S","rrnS","s-RNA","s-rRNA","s-rRNA","small ribosomal RNA","small ribosomal RNA","small ribosomal RNA","small ribosomal RNA","small ribosomal RNA","small ribosomal RNA subunit","small ribosomal RNA subunit","small ribosomal RNA subunit","small ribosomal RNA subunit RNA","small subinut ribosomal RNA","small subnit ribosomal RNA","small subumit ribosomal RNA","small subumit ribosomal RNA","small subunit","small subunit 12S ribosomal RNA","small subunit of 12S ribosomal RNA","small subunit of ribosomal RNA","small subunit ribosmal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA","small subunit ribosomal RNA gene","small-subunit ribosomal RNA subunit","srRNA","subunit ribosomal RNA"]
rna_16s = ["l-ribosomal RNA","16 ribosomal RNA","16 S ribosomal RNA","16S","16S (LSU) ribosomal RNA","16S large ribosomal RNA","16S large subunit ribosomal RNA","16S large subunit ribosomal RNA","16S large subunit ribosomal RNA","16S large subunit ribsomal RNA","16S mitochondrial ribosomal RNA","16S riboaomal RNA","16S ribosamal RNA","16S ribosmal RNA","16S ribosmal RNA","16S ribosoaml RNA","16S ribosomail RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA","16S ribosomal RNA gene","16S ribosomal RNA gene","16S ribosomal RNA homolog","16S ribosomal RNA protein","16S ribosomal RNA subunit","16S ribosomal RNA subunit RNA","16S ribosomal RNA, large subunit","16S ribosome RNA","16S ribosoml RNA","16S ribosommal RNA","16S ribososmal RNA","16S ribososomal RNA","16S ribsomal RNA","16S risbosomal RNA","16S rivbosomal RNA","16S-rRNA","l-rRNA","l-rRNA","l-rRNA","l6S ribosomal RNA","l6S ribosomal RNA","large 16S ribosomal RNA","large ribosomal RNA","large ribosomal RNA","large ribosomal RNA","large ribosomal RNA","large ribosomal RNA","large ribosomal RNA subnuit","large ribosomal RNA subunit","large ribosomal RNA subunit","large ribosomal RNA subunit","large ribosomal RNA subunit","large ribosomal RNA subunit","large ribosomal RNA subunit RNA","large subunit","large subunit 16S ribosomal RNA","large subunit of 16S ribosomal RNA","large subunit of ribosomal RNA","large subunit of ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA","large subunit ribosomal RNA (lrRNA)","large subunit ribosomal RNA gene","LrRNA","lrRNA","lsu ribosomal RNA","mitochondrial large ribosomal RNA subunit RNA","q6S ribosomal RNA","ribosomal RNA large subunit","ribosomal RNA large subunit","ribosomal RNA large subunit","ribosomal RNA, large subunit","rnl","rRNA large subunit","rRNA large subunit","rrnL"]

def main():
  # load the sequence file, it should have a bunch of mitogenomes
  seq_file = sys.argv[1]
  records = SeqIO.parse(seq_file,"genbank");
  # output tab-separated data to stdout
  print("accession\tspecies\tgene\tdescription\ttype\tnucleotides\tstart\tend")
  gene_desc = {}
  for seq in records:
    # pull species
    species = seq.annotations['organism']
    # get accession
    try:
      accession = seq.annotations['accessions'][0]
    except:
      accession = seq.name
    # look for coding regions or rRNAs
    for f in seq.features:
      if f.type == "CDS" or f.type == "rRNA":
        # pull out the info we care about
        cds = f.location.extract(seq.seq)
        start = int(f.location.start)
        end = int(f.location.end)
        description = f.qualifiers['product'][0]
        try:
          if f.type == "rRNA":
            raise Exception("rRNA")
          gene = f.qualifiers['gene'][0]
          if description not in gene_desc:
            gene_desc[description] = set()
          gene_desc[description].add(gene)
        except:
          if f.type == "CDS":
            if description in gene_desc and len(gene_desc[description]) == 1:
              gene = gene_desc[description].pop()
            else:
              gene="NA"
          # we'll trigger this exception because rRNAs don't have a 'gene' feature
          elif f.type == "rRNA": 
            if description in rna_12s:
              gene="12S"
            elif description in rna_16s:
              gene="16S"
            else:
              gene="NA"
        print(f"{accession}\t{species}\t{gene}\t{description}\t{f.type}\t{cds}\t{start}\t{end}")


if __name__ == "__main__":
  main()