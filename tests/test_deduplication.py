import itertools
import json
from pathlib import Path
import random
import tempfile
import unittest

from helixbusters.deduplication import deduplicate_bam, five_prime_position, group_umis
from tests.bam_helpers import HEADER, read, write_bam


class TestCoordinates(unittest.TestCase):
    def test_positive_negative_and_zero_boundary(self):
        self.assertEqual(five_prime_position(read(start=0)), 0)
        self.assertEqual(five_prime_position(read(start=0, flag=16)), 19)
        self.assertEqual(five_prime_position(read(start=100, flag=16)), 119)

    def test_reference_span_handles_internal_indels(self):
        self.assertEqual(five_prime_position(read(flag=16, cigar="10M3D10M")), 122)
        self.assertEqual(five_prime_position(read(flag=16, cigar="10M3I10M")), 119)
        self.assertEqual(five_prime_position(read(flag=16, cigar="10=2X8=")), 119)

    def test_strict_5prime_and_permitted_3prime_clipping(self):
        for cigar, flag in [("3S20M", 0), ("20M3S", 16), ("3H20M", 0),
                            ("20M3H", 16), ("2I20M", 0), ("20M2D", 16)]:
            with self.subTest(cigar=cigar, flag=flag):
                self.assertIsNone(five_prime_position(read(cigar=cigar, flag=flag)))
        self.assertEqual(five_prime_position(read(cigar="20M3S")), 100)
        self.assertEqual(five_prime_position(read(cigar="3S20M", flag=16)), 119)

    def test_aligned_policy_does_not_invent_unclipped_coordinate(self):
        self.assertEqual(five_prime_position(read(cigar="3S20M"), "aligned"), 100)
        self.assertEqual(five_prime_position(read(cigar="20M3S", flag=16), "aligned"), 119)
        self.assertIsNone(five_prime_position(read(cigar="10M100N10M"), "aligned"))
        self.assertIsNone(five_prime_position(read(flag=4)))


class TestUmiFamilies(unittest.TestCase):
    def test_exact_preserves_distinct_umis(self):
        families = group_umis({"AAAAAAAA": 10, "AAAAAAAC": 1})
        self.assertEqual(len(families), 2)
        self.assertEqual(sum(f.read_count for f in families), 11)

    def test_directional_corrects_errors_but_preserves_balanced_umis(self):
        families = group_umis({"AAAAAAAA": 10, "AAAAAAAC": 1}, "directional")
        self.assertEqual(len(families), 1)
        self.assertEqual(families[0].representative, "AAAAAAAA")
        self.assertEqual(families[0].read_count, 11)
        self.assertEqual(len(group_umis({"AAAAAAAA": 10, "AAAAAAAC": 9}, "directional")), 2)

    def test_directional_transitive_chain(self):
        families = group_umis({"AAAAAAAA": 10, "AAAAAAAC": 4, "AAAAAACC": 1}, "directional")
        self.assertEqual(len(families), 1)
        self.assertEqual(families[0].read_count, 15)

    def test_singleton_ties_and_input_order_are_deterministic(self):
        items = [("AAAAAAAA", 1), ("AAAAAAAC", 1), ("CCCCCCCC", 1)]
        expected = group_umis(dict(items), "directional")
        self.assertEqual(len(expected), 2)
        for permutation in itertools.permutations(items):
            self.assertEqual(group_umis(dict(permutation), "directional"), expected)

    def test_distance_zero_one_two(self):
        counts = {"AAAAAAAA": 10, "AAAAAACC": 1}
        self.assertEqual(len(group_umis(counts, "directional", 0)), 2)
        self.assertEqual(len(group_umis(counts, "directional", 1)), 2)
        self.assertEqual(len(group_umis(counts, "directional", 2)), 1)

    def test_against_bruteforce_directed_graph(self):
        # An independent all-pairs Hamming oracle validates neighbor enumeration
        # and shared descendants with dense, randomly generated graphs.
        rng = random.Random(2026)
        pool = ["".join(x) for x in itertools.product("ACGT", repeat=3)]
        for distance in (1, 2):
            for _ in range(20):
                counts = {u: rng.randint(1, 12) for u in rng.sample(pool, 25)}
                edges = {u: {v for v in counts if u != v
                             and sum(a != b for a, b in zip(u, v)) <= distance
                             and counts[u] >= 2 * counts[v] - 1} for u in counts}
                assigned, expected = set(), []
                for root in sorted(counts, key=lambda u: (-counts[u], u)):
                    if root in assigned:
                        continue
                    reachable = {root}
                    while True:
                        expanded = reachable | set().union(*(edges[u] for u in reachable))
                        if expanded == reachable:
                            break
                        reachable = expanded
                    members = reachable - assigned
                    expected.append((root, tuple(sorted(members)), sum(counts[u] for u in members)))
                    assigned |= members
                actual = [(f.representative, f.members, f.read_count)
                          for f in group_umis(counts, "directional", distance)]
                self.assertEqual(actual, expected)

    def test_invalid_counts_and_options(self):
        for counts in ({"": 1}, {"AAAAAAAN": 1}, {"AAAA": 1, "AAA": 1},
                       {"AAAA": 0}, {"AAAA": 1.5}):
            with self.subTest(counts=counts), self.assertRaises(ValueError):
                group_umis(counts)
        with self.assertRaises(ValueError):
            group_umis({}, "unknown")
        with self.assertRaises(ValueError):
            group_umis({}, max_distance=3)


class TestBamDeduplication(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def run_dedup(self, records, **options):
        bam = write_bam(self.root / "input.bam", records)
        return deduplicate_bam(bam, self.root / "families.tsv", self.root / "counts.bed",
                               output_molecules=self.root / "molecules.bed",
                               output_sites=self.root / "sites.tsv",
                               output_qc=self.root / "qc.json", **options)

    def test_strands_contigs_and_nearby_positions_stay_separate(self):
        records = [read(start=119), read(start=119), read(start=119, umi="CCCCCCCC"),
                   read(start=100, flag=16), read(start=120),
                   read(start=119, ref=1), read(start=119, ref=2),
                   read(start=119, ref=3), read(start=119, ref=4)]
        stats = self.run_dedup(records, method="directional")
        self.assertEqual(stats["accepted_reads"], 9)
        self.assertEqual(stats["deduplicated_molecules"], 8)
        self.assertEqual(stats["duplicate_reads"], 1)
        self.assertEqual(stats["strand_loci"], 7)
        self.assertEqual((self.root / "counts.bed").read_text().splitlines(), [
            "chr1\t119\t120\t3", "chr1\t120\t121\t1", "chrX\t119\t120\t1",
            "chrY\t119\t120\t1", "chrM\t119\t120\t1", "scaffold_1\t119\t120\t1"])
        bed = [line.split("\t") for line in (self.root / "molecules.bed").read_text().splitlines()]
        self.assertEqual(len(bed), 8)
        self.assertTrue(all(len(row) == 6 and row[4] == "0" for row in bed))
        self.assertEqual({r[5] for r in bed if r[0] == "chr1" and r[1] == "119"}, {"+", "-"})
        self.assertEqual(json.loads((self.root / "qc.json").read_text()), stats)

    def test_reverse_ends_join_even_with_different_alignment_starts(self):
        stats = self.run_dedup([read(start=100, flag=16, cigar="120M"),
                                read(start=150), read(start=200, flag=16)])
        self.assertEqual(stats["deduplicated_molecules"], 2)
        self.assertEqual((self.root / "counts.bed").read_text(), "chr1\t150\t151\t1\nchr1\t219\t220\t1\n")

    def test_directional_family_support_and_reproducible_outputs(self):
        records = [read()] * 10 + [read(umi="AAAAAAAC")] + [read(umi="CCCCCCCC")] * 5
        stats = self.run_dedup(records, method="directional")
        self.assertEqual(stats["exact_molecules"], 3)
        self.assertEqual(stats["deduplicated_molecules"], 2)
        self.assertEqual(stats["umi_groups_merged"], 1)
        before = (self.root / "families.tsv").read_bytes()
        self.assertIn(b"AAAAAAAA\t11\n", before)
        self.run_dedup(list(reversed(records)), method="directional")
        self.assertEqual((self.root / "families.tsv").read_bytes(), before)

    def test_filters_umi_parsing_and_duplicate_flags(self):
        records = [read(name="many_underscores_in_original_name_AAAAAAAA"),
                   read(name="no_suffix", rx="CCCCCCCC"),
                   read(flag=1024, umi="GGGGGGGG"),
                   read(flag=256), read(flag=2048), read(flag=512), read(mapq=19),
                   read(mapq=255), read(cigar="2S20M"), read(name="missing"),
                   read(umi="AAAAAAAN"), read(umi="AAAA"),
                   read(umi="TTTTTTTT", rx="NNNNNNNN"),
                   read(start=-1, ref=-1, flag=4, cigar=None)]
        stats = self.run_dedup(records, min_mapq=20)
        self.assertEqual(stats["accepted_reads"], 3)
        self.assertEqual(stats["skipped"], {
            "secondary": 1, "supplementary": 1, "qc_fail": 1, "low_mapq": 1,
            "unknown_mapq": 1, "ambiguous_five_prime": 1, "missing_umi": 1,
            "invalid_umi": 3, "unmapped": 1})
        self.assertEqual(stats["input_reads"], stats["accepted_reads"] + sum(stats["skipped"].values()))

    def test_paired_selection_is_explicit(self):
        records = [read(flag=65), read(flag=145)]
        with self.assertRaisesRegex(ValueError, "explicit read_selection"):
            self.run_dedup(records)
        first = self.run_dedup(records, read_selection="read1")
        self.assertEqual(first["accepted_reads"], 1)
        self.assertEqual(first["skipped"], {"other_mate": 1})
        self.assertIn("\t100\t101\t", (self.root / "counts.bed").read_text())
        self.run_dedup(records, read_selection="read2")
        self.assertIn("\t119\t120\t", (self.root / "counts.bed").read_text())
        with self.assertRaisesRegex(ValueError, "paired-read flags"):
            self.run_dedup([read()], read_selection="read1")

    def test_all_missing_umis_fail_without_replacing_existing_outputs(self):
        (self.root / "counts.bed").write_text("previous successful result\n")
        with self.assertRaisesRegex(ValueError, "No usable UMIs"):
            self.run_dedup([read(name="no_umi_suffix")])
        self.assertEqual((self.root / "counts.bed").read_text(), "previous successful result\n")
        self.assertEqual(list(self.root.glob(".*.tmp")), [])
        self.assertFalse((self.root / "qc.json").exists())

    def test_empty_bam_writes_zero_qc(self):
        stats = self.run_dedup([])
        self.assertEqual(stats["input_reads"], 0)
        self.assertEqual(stats["deduplicated_molecules"], 0)
        self.assertEqual((self.root / "molecules.bed").read_text(), "")

    def test_unsorted_bam_rejected_even_if_header_claims_sorted(self):
        for records in ([read(start=200), read(start=100)],
                        [read(ref=1), read(ref=0)],
                        [read(start=-1, ref=-1, flag=4, cigar=None), read()]):
            bam = write_bam(self.root / "unsorted.bam", records, sort=False)
            with self.assertRaisesRegex(ValueError, "coordinate-sorted"):
                deduplicate_bam(bam, self.root / "a.tsv", self.root / "b.bed")
            self.assertFalse((self.root / "a.tsv").exists())

    def test_input_and_output_aliases_rejected(self):
        bam = write_bam(self.root / "input.bam", [read()])
        original = bam.read_bytes()
        with self.assertRaisesRegex(ValueError, "distinct"):
            deduplicate_bam(bam, bam, self.root / "b.bed")
        self.assertEqual(bam.read_bytes(), original)

    def test_multiple_samples_rejected(self):
        header = dict(HEADER, RG=[{"ID": "a", "SM": "sample_a"}, {"ID": "b", "SM": "sample_b"}])
        bam = write_bam(self.root / "input.bam", [read()], header=header)
        with self.assertRaisesRegex(ValueError, "multiple SM"):
            deduplicate_bam(bam, self.root / "a.tsv", self.root / "b.bed")


if __name__ == "__main__":
    unittest.main()
