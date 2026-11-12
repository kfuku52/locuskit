# LOCUSKIT neighborhood tables

LOCUSKIT uses explicit gene-to-species mappings and reports what was observed
around each requested gene. It accepts precomputed homology groups or DIAMOND hits and can run DIAMOND
on a supplied protein FASTA. Coordinate-only extraction uses the Python standard library and the
inherited FractBias annotation reader.

## Generate and validate

A genome manifest has the tab-separated columns species_id, assembly_id and
gff. Paths resolve relative to the manifest. Optional feature and attribute
columns default to gene and ID. Use one assembly per species. The locus table
includes annotated genes without sequence data. If a tree uses representative
transcript IDs, supply an alias table with species_id, locus_id and gene_id to
map physical loci to those exact tree labels. Collisions and unknown source loci
are errors.

~~~bash
locuskit loci --manifest_tsv genomes.tsv --aliases_tsv aliases.tsv --out_file loci.tsv
locuskit neighbors --focal_tsv mapping.tsv --loci_tsv loci.tsv --groups_tsv homology.tsv \
  --window 20 --out_file family.context.tsv
locuskit validate --in_file family.context.tsv --focal_tsv mapping.tsv
locuskit export-legacy --in_file family.context.tsv --out_file family_synteny.tsv
~~~

For multiple families using the same annotation, pass a TSV manifest with
`focal_tsv` and `out_file` columns. Both paths resolve relative to the manifest.

~~~bash
locuskit neighbors-batch --manifest_tsv families.tsv --loci_tsv loci.tsv \
  --protein_fasta proteins.faa --window 20
~~~

The manifest command reads the locus table once and uses the same neighborhood
and homology options as `neighbors`. Each family uses its own set of proteins
for DIAMOND. With a shared FASTA, the batch first streams the proteome into a
temporary FASTA containing only proteins found in the requested windows.
LOCUSKIT rejects overlapping input and output paths and repeated
outputs before writing. Each output is replaced atomically. If a later job
fails, earlier completed tables remain available.

The focal mapping uses gene and species as in aioi. The locus table has
species_id, assembly_id, gene_id, seqid, start, end and strand. Existing
normalized coordinates can be passed directly. Use positive, one-based,
inclusive intervals and one row per gene locus. Strand is +, -, . or ?.

The optional homology table uses species_id, gene_id and group_id. A group
membership is scoped to its input table. To use existing DIAMOND hits instead,
pass --hits_tsv with headerless qseqid, sseqid and evalue columns. The default
cutoff is 0.01. All selected neighbor sequence IDs must be unique across
species. LOCUSKIT takes connected components of admitted hits among observed
neighbors. A gene that occurs in the hit file without a retained link gets a
singleton group. A gene with no hit record keeps its coordinate and an empty
group. To run the search, use --protein_fasta and optionally --threads. LOCUSKIT
selects the neighborhood proteins, runs DIAMOND from PATH in temporary
storage, and retains their groups. Genes missing from the protein FASTA keep
their coordinates and an empty group. These groups describe sequence similarity. Locus conservation, event
assignment and homology probabilities need further evidence.

Plain UTF-8 TSV and gzip-compressed TSV are accepted. Named outputs replace
previous outputs atomically after validation. Output paths may not alias
inputs.

## Shared context v1

All rows begin with schema_version 1. A record_type of focal records the
requested ID and species, annotation status, requested window and the observed
number of annotated positions on each side when known. A record_type of
neighbor adds the neighbor ID and signed rank relative to the transcription
direction of its focal gene. Rank -1 is the nearest upstream gene and +1 is the
nearest downstream gene. Positions are ordered by start, end and gene ID on
each contig. Intervals describe the gene on the row. The assembly and contig are
shared between a focal gene and its neighbors. Validation checks that every
appearance of a locus has the same coordinates and that ranks follow the
focal strand and the genomic order. Coverage counts reflect supplied
annotation and can be zero. A truncated window can reflect the available
annotation near a contig boundary.

| Field | Meaning |
|---|---|
| focal_id, species_id | Exact identifiers from the focal mapping. |
| assembly_id, seqid, start, end, strand | Coordinates of the focal or neighbor in its assembly. All five are blank if unavailable. |
| window, observed_upstream, observed_downstream | Requested and observed annotated ranks. Focal rows only. |
| status | observed, annotation_missing, locus_missing, orientation_unknown or legacy_unknown. Focal rows only. |
| offset, neighbor_id | Signed rank and gene identifier. Neighbor rows only. |
| group_id, group_size | Optional sequence similarity group and its unique gene count. |
| evalue_cutoff, group_method | Optional cutoff and source. Methods are provided_groups, similarity_hits, diamond and legacy. |

annotation_missing means the focal species had no rows in the supplied locus
table. locus_missing means the species was represented but the focal ID was
absent. orientation_unknown retains the known interval, with no signed ranks.
legacy_unknown preserves old offsets and group memberships while leaving
coordinates, window size and coverage blank. Missing context is unknown
evidence for downstream callers.

## Use existing GeneGalleon output

~~~bash
locuskit import-legacy --in_file family_synteny.tsv --focal_tsv mapping.tsv \
  --out_file family.context.tsv
locuskit export-legacy --in_file family.context.tsv --out_file family_tree_plot.tsv
~~~

Both seven-column legacy output and the variant with evalue_cutoff can be
imported. The original producer already attaches species names. The optional
mapping covers older tables without a species column and adds missing focal
records with legacy_unknown status. The legacy export contains only grouped
neighbors, which the existing tree_plot reader can use directly.

## FractBias compatibility

The LOCUSKIT distribution includes the kffractbias Python package and legacy
command. Use locuskit fractbias with the existing subcommands and options.
Applications that inspect the old distribution metadata can install the
metadata-only project under compat/kffractbias alongside LOCUSKIT. In a migration
environment, remove the previous full kffractbias distribution before installing
LOCUSKIT and its metadata shim so that a single wheel owns the Python files.
The [FractBias user guide](fractbias.md) describes the analyses and citation.
