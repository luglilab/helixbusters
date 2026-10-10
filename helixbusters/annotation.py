"""Strand-aware GTF feature partitions and exact half-open interval annotation."""
from bisect import bisect_right
from collections import Counter, defaultdict
import gzip
import hashlib
from pathlib import Path
import re

FEATURES = ('promoter', 'exon', 'intron', 'intergenic')


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return merged


class GTFAnnotation:
    def __init__(self, path, header, upstream=2000, downstream=500, genome=None):
        if upstream < 0 or downstream < 1:
            raise ValueError('Promoter upstream must be >=0 and downstream >=1')
        self.header = dict(header)
        self.genes = {}
        self.segments = {}
        self.starts = {}
        self.gene_candidates = {}
        aliases = {}
        for chrom in self.header:
            aliases[chrom] = chrom
            if re.fullmatch(r'(chr)?([0-9]+|X|Y)', chrom):
                aliases[chrom.removeprefix('chr')] = chrom
                aliases['chr' + chrom.removeprefix('chr')] = chrom
        opener = gzip.open if str(path).endswith('.gz') else open
        ignored = 0
        header_comments = []
        with opener(path, 'rt') as handle:
            for number, line in enumerate(handle, 1):
                if line.startswith('#'):
                    if len(header_comments) < 20:
                        header_comments.append(line.rstrip())
                    if genome:
                        from helixbusters.genomes import normalize_genome
                        builds = re.findall(r'\b(GRCh37|GRCh38|GRCm38|GRCm39|hg19|hg38|mm10|mm39)\b', line)
                        if any(normalize_genome(build) != normalize_genome(genome) for build in builds):
                            raise ValueError('GTF header genome build disagrees with declared genome')
                    continue
                if not line.strip():
                    continue
                fields = line.rstrip('\n').split('\t')
                if len(fields) != 9:
                    raise ValueError(f'GTF row {number} must have nine fields')
                original, _, feature, raw_start, raw_end, _, strand, _, attributes = fields
                if original not in aliases:
                    ignored += 1
                    continue
                if feature not in {'gene', 'transcript', 'exon'}:
                    continue
                chrom = aliases[original]
                start, end = int(raw_start) - 1, int(raw_end)
                if strand not in {'+', '-'} or not 0 <= start < end <= self.header[chrom]:
                    raise ValueError(f'Invalid or reference-incompatible GTF coordinates at row {number}')
                attrs = dict(re.findall(r'(\w+)\s+"([^"]*)"', attributes))
                gene_id = attrs.get('gene_id')
                if not gene_id:
                    raise ValueError(f'Missing gene_id at GTF row {number}')
                key = (chrom, gene_id)
                gene = self.genes.setdefault(key, {'start': start, 'end': end, 'strand': strand,
                    'name': attrs.get('gene_name', gene_id), 'biotype': attrs.get('gene_type', attrs.get('gene_biotype', 'unknown')),
                    'exons': [], 'tss': set()})
                if gene['strand'] != strand:
                    raise ValueError(f'Gene has inconsistent strand: {gene_id}')
                gene['start'], gene['end'] = min(gene['start'], start), max(gene['end'], end)
                if feature == 'exon':
                    gene['exons'].append((start, end))
                if feature == 'transcript':
                    gene['tss'].add(start if strand == '+' else end - 1)
        if not self.genes:
            raise ValueError('No GTF genes match the reference chromosomes')
        events = defaultdict(list)
        for (chrom, gene_id), gene in self.genes.items():
            if not gene['exons']:
                raise ValueError(f'Gene {gene_id} has no exon annotation; use a complete gene/transcript/exon GTF')
            gene['exons'] = merge_intervals(gene['exons'])
            if not gene['tss']:
                gene['tss'].add(gene['start'] if gene['strand'] == '+' else gene['end'] - 1)
            intervals = [('body', [(gene['start'], gene['end'])]), ('exon', gene['exons'])]
            promoters = []
            for tss in gene['tss']:
                start, end = ((tss - upstream, tss + downstream) if gene['strand'] == '+'
                              else (tss - downstream + 1, tss + upstream + 1))
                start, end = max(0, start), min(self.header[chrom], end)
                if start < end:
                    promoters.append((start, end))
            intervals.append(('promoter', merge_intervals(promoters)))
            for kind, rows in intervals:
                for start, end in rows:
                    events[chrom].extend([(start, kind, gene_id, 1), (end, kind, gene_id, -1)])
        for chrom, length in header:
            active = {kind: Counter() for kind in ('promoter', 'exon', 'body')}
            grouped = defaultdict(list)
            for position, kind, gene_id, delta in events[chrom]:
                grouped[position].append((kind, gene_id, delta))
            previous = 0
            rows = []
            candidates = []
            for position in sorted(set(grouped) | {length}):
                if previous < position:
                    category = next((kind for kind in ('promoter', 'exon', 'body') if active[kind]), None)
                    candidates.append(tuple(sorted(active[category])) if category else ())
                    category = 'intron' if category == 'body' else category or 'intergenic'
                    genes = tuple(sorted(set().union(*(set(values) for values in active.values()))))
                    rows.append((previous, position, category, genes))
                for kind, gene_id, delta in grouped[position]:
                    active[kind][gene_id] += delta
                    if active[kind][gene_id] == 0:
                        del active[kind][gene_id]
                previous = position
            self.segments[chrom] = rows
            self.starts[chrom] = [row[0] for row in rows]
            self.gene_candidates[chrom] = candidates
        hasher = hashlib.sha256()
        with Path(path).open('rb') as handle:
            for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
                hasher.update(block)
        self.provenance = {'gtf': str(path), 'gtf_sha256': hasher.hexdigest(), 'genome': genome,
            'gtf_header_comments': header_comments,
            'genes': len(self.genes), 'ignored_rows_unknown_contigs': ignored,
            'promoter_upstream': upstream, 'promoter_downstream': downstream,
            'promoter_policy': 'all annotated transcript TSS; gene-span TSS fallback if transcripts absent; strand-aware relative interval [-upstream, downstream)',
            'priority': list(FEATURES), 'gene_body_definition': 'union of gene spans; exon is union across all isoforms; intron is residual body outside promoters and any exon',
            'association_policy': 'all overlapping gene-body or promoter gene IDs; associations are not proof of regulation',
            'coordinates': 'input GTF 1-based inclusive converted to 0-based half-open'}

    def overlaps(self, chrom, start, end):
        if chrom not in self.header or not 0 <= start < end <= self.header[chrom]:
            raise ValueError('Query lies outside the annotation reference')
        rows = self.segments[chrom]
        index = max(0, bisect_right(self.starts[chrom], start) - 1)
        while index < len(rows) and rows[index][0] < end:
            left, right, category, genes = rows[index]
            if right > start:
                yield min(end, right) - max(start, left), category, genes
            index += 1

    def annotate(self, chrom, start, end):
        coverage = dict.fromkeys(FEATURES, 0)
        genes = set()
        for length, category, gene_ids in self.overlaps(chrom, start, end):
            coverage[category] += length
            genes.update(gene_ids)
        if sum(coverage.values()) != end - start:
            raise ValueError('Feature partition does not conserve query length')
        dominant = max(FEATURES, key=lambda kind: coverage[kind])
        return dominant, coverage, sorted(genes)

    def gene_fields(self, chrom, ids):
        return (';'.join(ids), ';'.join(self.genes[(chrom, gene)]['name'] for gene in ids),
                ';'.join(sorted({self.genes[(chrom, gene)]['biotype'] for gene in ids})))

    def site_gene_candidates(self, chrom, position):
        """Genes in the highest-priority feature at a 1-bp DSB endpoint."""
        if chrom not in self.header or not 0 <= position < self.header[chrom]:
            raise ValueError('Site lies outside the annotation reference')
        index = bisect_right(self.starts[chrom], position) - 1
        return self.segments[chrom][index][2], self.gene_candidates[chrom][index]
