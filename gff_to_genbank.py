#!/usr/bin/env python3
"""gff_to_genbank.py
Parses GFF (Prodigal/GeneMark), BLAST tabular outputs (SwissProt/MIL-1), tRNAscan-SE output and rRNA BLAST results,
and writes an annotated GenBank file (GenBank format via Biopython).

Usage example:
  python3 /mnt/data/scripts/gff_to_genbank.py \
    --genome data/scaffolds.fasta \
    --gff results/prodigal.gff \
    --blast_sp results/scaffolds.hits_from_SwissProt.txt \
    --blast_mil results/scaffolds.hits_from_MIL_1.txt \
    --trna results/trnascan.out \
    --rrna results/rrna_blast.tsv \
    --out results/GENOME.gbk
"""

import argparse, os, sys, re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def parse_args():
    p = argparse.ArgumentParser(description='Convert GFF+BLAST+tRNA+rRNA -> GenBank')
    p.add_argument('--genome', required=True, help='Input genome FASTA (contigs)')
    p.add_argument('--gff', default=None, help='GFF file with gene/CDS features (Prodigal/GeneMark)')
    p.add_argument('--blast_sp', default=None, help='BLASTp vs SwissProt (tabular)')
    p.add_argument('--blast_mil', default=None, help='BLASTp vs MIL-1 (tabular)')
    p.add_argument('--trna', default=None, help='tRNAscan-SE output file')
    p.add_argument('--rrna', default=None, help='rRNA BLAST TSV')
    p.add_argument('--out', required=True, help='Output GenBank file path')
    return p.parse_args()

def read_blast_tabular(path):
    hits = {}
    if not path or not os.path.exists(path):
        return hits
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 2:
                continue
            q = cols[0]; s = cols[1]
            # try to parse known columns; allow variable length
            pident = float(cols[2]) if len(cols) > 2 and cols[2] != '' else None
            length = int(cols[3]) if len(cols) > 3 and cols[3] != '' else None
            evalue = cols[4] if len(cols) > 4 else ''
            bits = cols[5] if len(cols) > 5 else ''
            stitle = cols[6] if len(cols) > 6 else ''
            entry = dict(sseqid=s, pident=pident, length=length, evalue=evalue, bitscore=bits, stitle=stitle)
            # keep best by bitscore if numeric
            if q not in hits:
                hits[q] = entry
            else:
                try:
                    if float(bits) > float(hits[q].get('bitscore', 0)):
                        hits[q] = entry
                except:
                    pass
    return hits

def parse_gff(gff_path):
    feats = []
    if not gff_path or not os.path.exists(gff_path):
        return feats
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9:
                continue
            contig, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype.lower() not in ('cds','cdregion','gene'):
                continue
            attrd = {}
            for a in attrs.split(';'):
                if '=' in a:
                    k,v = a.split('=',1)
                    attrd[k]=v
            feats.append(dict(contig=contig, type=ftype, start=int(start), end=int(end), strand=strand, attrs=attrd))
    return feats

def parse_trnascan(path):
    trnas = []
    if not path or not os.path.exists(path):
        return trnas
    with open(path) as fh:
        for line in fh:
            if line.strip()=='' or line.startswith('Seq') or line.startswith('Name') or line.startswith('--------'):
                continue
            cols = line.split()
            if len(cols) < 5:
                continue
            # common layout fallback: name, #, begin, end, type, score...
            try:
                name = cols[0]; begin = int(cols[2]); end = int(cols[3]); ttype = cols[4]
                trnas.append(dict(seqname=name, start=begin, end=end, strand='+', type=ttype))
            except:
                # try other positions
                try:
                    name = cols[0]; begin = int(cols[1]); end = int(cols[2]); ttype = cols[3]
                    trnas.append(dict(seqname=name, start=begin, end=end, strand='+', type=ttype))
                except:
                    continue
    return trnas

def parse_rrna(path):
    rr = []
    if not path or not os.path.exists(path):
        return rr
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 9:
                continue
            q,s,pident,length,qstart,qend,sstart,send,evalue = cols[:9]
            try:
                sstart=int(sstart); send=int(send); pidentf = float(pident)
            except:
                continue
            rr.append(dict(qseqid=q, sseqid=s, pident=pidentf, sstart=sstart, send=send))
    return rr

def write_gbk(out, genome_records, cds, trna, rrna, blast_sp, blast_mil):
    """
    Build SeqRecord objects with features and write GenBank.
    This version normalizes coordinates (start <= end) and sets strand correctly
    even if GFF provides start > end or BLAST returns reversed coordinates.
    """
    from Bio import SeqIO
    recs = []
    cds_by_contig = {}
    for c in cds:
        cds_by_contig.setdefault(c['contig'], []).append(c)

    for cid, seqrec in genome_records.items():
        rec = SeqRecord(Seq(str(seqrec.seq)), id=seqrec.id, name=seqrec.id, description=seqrec.description)
        feats = []

        # CDS features
        for c in cds_by_contig.get(cid, []):
            # read raw coordinates
            raw_start = int(c['start'])
            raw_end = int(c['end'])
            # normalize: start <= end in Python indexing
            start0 = min(raw_start, raw_end) - 1   # convert to 0-based
            end0 = max(raw_start, raw_end)         # end is exclusive in FeatureLocation
            # determine strand: prefer explicit GFF strand if valid, else infer from coords
            strand_symbol = c.get('strand', '+')
            if strand_symbol == '+':
                strand = 1
            elif strand_symbol == '-':
                strand = -1
            else:
                # infer from coordinates (if start > end -> negative)
                strand = 1 if raw_end >= raw_start else -1

            # build qualifiers
            qual = {}
            locus = c['attrs'].get('ID') or c['attrs'].get('locus_tag') or c['attrs'].get('Name')
            if locus: qual['locus_tag'] = [locus]
            protid = c['attrs'].get('protein_id')
            if protid: qual['protein_id'] = [protid]

            # try to attach annotation from BLAST (SwissProt or MIL)
            qid = locus or protid
            prod = None
            if qid and qid in blast_sp:
                prod = blast_sp[qid].get('stitle') or blast_sp[qid].get('sseqid')
            elif qid and qid in blast_mil:
                prod = blast_mil[qid].get('stitle') or blast_mil[qid].get('sseqid')
            if prod:
                qual['product'] = [prod if len(prod) < 250 else prod[:247] + '...']

            # add CDS feature (FeatureLocation requires start <= end)
            loc = FeatureLocation(start0, end0, strand=strand)
            feats.append(SeqFeature(location=loc, type='CDS', qualifiers=qual))

        # tRNAs
        for t in trna:
            # map by seqname if provided
            if t.get('seqname') and t.get('seqname') not in (cid, '.', '', None):
                continue
            ts = int(t['start']); te = int(t['end'])
            tstart = min(ts, te) - 1
            tend = max(ts, te)
            tstrand = 1 if t.get('strand', '+') == '+' else -1
            loc = FeatureLocation(tstart, tend, strand=tstrand)
            qual = {'product': [t.get('type', 'tRNA')], 'note': [t.get('type', 'tRNA')]}
            feats.append(SeqFeature(location=loc, type='tRNA', qualifiers=qual))

        # rRNAs
        for r in rrna:
            if r['sseqid'] != cid:
                continue
            rs = int(r['sstart']); re = int(r['send'])
            rstart = min(rs, re) - 1
            rend = max(rs, re)
            rstrand = 1 if re - rs >= 0 else -1
            loc = FeatureLocation(rstart, rend, strand=rstrand)
            qual = {'product': [r['qseqid']], 'note': [f"pident={r['pident']}"]}
            feats.append(SeqFeature(location=loc, type='rRNA', qualifiers=qual))

        rec.features = feats
        recs.append(rec)
    for rec in recs:
        # обязательный минимум для Biopython GenBank writer
        rec.annotations.setdefault('molecule_type', 'DNA')   # 'DNA' или 'RNA' по смыслу
        rec.annotations.setdefault('topology', 'linear')     # 'linear' или 'circular'
        # полезные дополнительные поля (можно заполнить конкретными значениями позже)
        rec.annotations.setdefault('organism', 'unknown')
        rec.annotations.setdefault('source', rec.id)
        rec.annotations.setdefault('data_file_division', 'PLN')
        rec.annotations.setdefault('date', '01-JAN-2000')

    # write output
    with open(out, 'w') as fh:
        SeqIO.write(recs, fh, 'genbank')
    return out


def main():
    args = parse_args()
    if not os.path.exists(args.genome):
        sys.exit('Genome FASTA not found: %s' % args.genome)
    genome_records = {rec.id:rec for rec in SeqIO.parse(args.genome, 'fasta')}
    cds = parse_gff(args.gff) if args.gff else []
    blast_sp = read_blast_tabular(args.blast_sp) if args.blast_sp else {}
    blast_mil = read_blast_tabular(args.blast_mil) if args.blast_mil else {}
    trna = parse_trnascan(args.trna) if args.trna else []
    rrna = parse_rrna(args.rrna) if args.rrna else []
    out = write_gbk(args.out, genome_records, cds, trna, rrna, blast_sp, blast_mil)
    print('Wrote GenBank to', out)

if __name__ == '__main__':
    main()
