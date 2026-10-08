# Micro-C: FASTQ to `.cool`, with QC and resolution selection

Companion to `microc_fastq_to_cool.sh`, version 2, 8 October 2026.

## Scope and interpretation

The script processes one **bulk, genome-wide, paired-end Micro-C library without UMIs**, starting with `R1.fastq.gz`, `R2.fastq.gz`, and an uncompressed reference FASTA. It produces filtered contact pairs, raw counts and balancing weights in a standalone `.cool`, a multiresolution `.mcool`, and QC reports. It is not a single-cell, capture Micro-C, or UMI-aware workflow.

**A finished script is not equivalent to a passed experiment.** Wet-lab digestion, sample identity, biological plausibility, and replicate agreement require information beyond one FASTQ pair. The script lists these as `NOT_REVIEWED` in `qc/manual_review.tsv`; complete that review before interpreting the maps.

Three kinds of numbers appear below:

- **Hard requirements:** file integrity and count conservation; violations stop processing.
- **Processing settings:** reproducible choices that determine which reads or bins are retained. They are not universal sample acceptance criteria.
- **Contextual benchmarks:** published values for particular protocols, input amounts, and depths. Falling outside them prompts investigation rather than automatic rejection.

All percentages must be compared with the **same denominator, genome, chromosome set, filters, and sequencing depth**. Several metrics have no scientifically defensible universal cutoff; this guide states that explicitly.

## Install and run

Linux/Bash, substantial temporary disk space, and a conda/mamba environment are expected. Sorting, pairs, and the finest matrix can consume substantial storage. `THREADS` is per major tool, so concurrent processes can use more than that total; `SORT_MEMORY` is the sort buffer, not total RAM.

```bash
mamba create -n microc -c conda-forge -c bioconda --strict-channel-priority \
  python=3.11 'numpy<2' 'pandas<3' pairtools=1.1.2 cooler=0.10.4 \
  bwa samtools htslib fastqc multiqc matplotlib cutadapt
conda activate microc

SAMPLE=sample1 OUTDIR=sample1_microc THREADS=8 RESOLUTION=10000 \
  bash microc_fastq_to_cool.sh /path/to/genome.fa R1.fastq.gz R2.fastq.gz
```

Use a **new output directory** for each run. Run biological replicates separately. For multiple sequencing lanes from the same physical library, concatenate all R1 lane files in one order and all corresponding R2 files in exactly the same order before processing; cross-lane duplicates then get removed together. Do not mix independently prepared libraries before deduplication.

For human/mouse and other references containing alternative haplotypes, supply a text file of primary nuclear chromosomes, one FASTA-matching name per line:

```bash
CHROMS_FILE=primary_nuclear_chromosomes.txt \
EXCLUDE_CHROMS=chrM,MT,M SAMPLE=sample1 OUTDIR=sample1_primary \
  bash microc_fastq_to_cool.sh /path/to/genome.fa R1.fastq.gz R2.fastq.gz
```

Adapt organelle names for your organism, including chloroplast contigs when relevant. Alignment still uses the full supplied FASTA. Existing BWA index files are reused; if absent, the script builds an index under the output directory. When setting `BWA_INDEX`, ensure it was built from this exact assembly.

## Step 0 — Configuration and reproducibility

The script checks dependencies and parameter formats and records settings, tool versions, Python packages, and logs. The requested `RESOLUTION` is included in `RESOLUTIONS`; every requested bin size must be a positive integer multiple of the smallest size.

**Keep:** valid inputs and consistent reference/annotation builds. A mismatched chromosome name or assembly is an input error, not a low-quality biological sample. Check reference provenance manually; matching names alone cannot prove matching assemblies.

## Step 1 — FASTQ integrity, sequence QC, and optional trimming

Both complete compressed files are read to check record structure, mate identities, read counts, and gzip integrity. FastQC runs on both mates and, if trimming is enabled, on both trimmed mates. Every FastQC module status is preserved in `fastqc_modules.tsv`.

| Metric | Value to keep / trigger | Why this value; meaning outside it | Action |
|---|---|---|---|
| Paired-file integrity | **Zero** malformed records or mismatched mate IDs; equal, nonzero pair counts | Pair reconstruction assumes corresponding reads occur in the same order. Any violation can join unrelated molecules or lose data. | Stop; repair or regenerate paired files. Do not independently discard reads from one mate. |
| Base quality | Most bases at **Q30 or higher** is desirable, with no universal required percentage | Q30 represents nominal 0.1% base-call error; Q20 represents 1%. Low-quality tails reduce confident mapping. These are sequencing metrics, not Micro-C-specific grades. | Inspect cycle profiles and post-trim mapping yield. |
| FastQC per-base quality | WARN: median **<25** or lower quartile **<10**; FAIL: median **<20** or lower quartile **<5**, at any position [1] | These are FastQC's module defaults. The lower quartile detects a substantial subset of poor reads. A flag does not by itself invalidate the contact library. | Review affected cycles and mates; consider justified tail trimming. |
| Adapter content | FastQC WARN **>5%**, FAIL **>10%**, for a detected adapter at any position [2] | Short Micro-C inserts can cause sequencing into adapters. A high value can reflect insert length rather than failed ligation. | Check the library trace and actual adapter sequence; optionally trim and compare yield. |
| Per-base N content | FastQC WARN **>5%**, FAIL **>20%** at any position [3] | Ambiguous calls reduce usable sequence. These are positional thresholds, not thresholds for the aggregate `pct_N` column. | Investigate cycle/sequencing problems and mapping loss. |
| GC, sequence content, duplication, overrepresented sequences | **No universal pass range** | MNase bias, genome composition, recurrent contact ends, adapters, and contamination can all affect these modules. FastQC sequence duplication is not the coordinate duplicate fraction used later. | Identify flagged sequences and compare with matched libraries; investigate unexpected organisms if indicated. |

Trimming is **off by default**, since BWA performs local alignment. Enable it only with the actual kit's 3′ sequencing adapters:

```bash
TRIM_ADAPTERS=1 ADAPTER_R1=ACTUAL_R1_ADAPTER ADAPTER_R2=ACTUAL_R2_ADAPTER \
  OUTDIR=sample1_trimmed bash microc_fastq_to_cool.sh genome.fa R1.fastq.gz R2.fastq.gz
```

Replace the adapter placeholders before running. Cutadapt uses paired processing, Q20 tail trimming, and a default minimum retained length of **30 nt**; it removes the pair if either mate becomes too short. These are adjustable processing settings, not established Micro-C acceptance thresholds. Higher stringency sacrifices more reads; lower stringency retains more potentially uninformative sequence. Review `cutadapt.json`, its log, post-trim FastQC, and the fraction of raw pairs entering alignment. Do not trim genomic ligation junctions as if they were sequencing adapters. [4]

### Essential wet-lab QC accompanying Step 1

Inspect the **pre-ligation MNase digestion trace**. A published bulk Micro-C-XL protocol targets roughly **70–90% mononucleosomes**, with little subnucleosomal smear. Mononucleosomal DNA is approximately **150 bp before ligation**. These are protocol-specific benchmarks, not a target for the final adapter-containing library trace. [5]

- Too many di-/oligonucleosomes suggests insufficient digestion, which can produce unligated dinucleosome contamination and reduce effective spatial resolution.
- Excess short subnucleosomal material suggests overdigestion and possible loss or bias of informative chromatin fragments. A mononucleosome fraction above 90% alone does not prove overdigestion: inspect the full trace and follow the kit's specification.
- Check the final library for adapter dimers and the protocol's expected insert distribution. Compare digestion across experimental conditions so MNase differences are not mistaken for biology.

FASTQ read length and genomic separation of paired ends cannot substitute for these traces: ligated ends may originate far apart on a chromosome or on different chromosomes.

## Step 2 — Reference and chromosome selection

FASTA lengths define both parsing and matrix coordinates. By default, contigs named `chrM`, `MT`, or `M` are excluded from the matrix; other contigs remain unless `CHROMS_FILE` restricts them.

**Keep:** chromosome names and sizes that exactly match the alignment reference. Organelle contamination has **no universal acceptable percentage** across organisms and preparation methods. A high fraction means sequencing effort was spent outside the intended nuclear assay and can suggest poor nuclear enrichment. The report counts pairs touching explicitly excluded chromosomes before exclusion. Additional contigs omitted through `CHROMS_FILE` also reduce nuclear yield, but are not included in that particular organelle/excluded-name metric.

## Step 3 — Alignment and contact reconstruction

`bwa mem -SP` aligns the mates without normal paired-end rescue assumptions. Supplementary alignments are preserved for `pairtools parse` to interpret ligation-spanning reads. Contacts are then sorted for deduplication; an ordinary coordinate-sorted BAM is not used as parser input.

| Setting / metric | Value | Explanation and interpretation |
|---|---|---|
| Mapping quality | **MAPQ ≥30 at both ends** | Recommended stringent unique-mapping threshold in pairtools protocols. Nominally corresponds to ≤0.1% mapping error per alignment, subject to BWA's calibration and reference completeness; it is not a measured false-contact rate. Lower values admit more ambiguous locations; raising the threshold loses more contacts. [6] |
| Complex read walks | **`5unique`** | Reports one contact per sequenced molecule, selecting outer unique alignments when needed. Some complex-walk contacts may be indirect. Use `WALKS_POLICY=mask` for a more conservative comparison, and keep the policy identical across samples. [6] |
| Rescue molecule size | **750 bp** | Pairtools' default bound for rescuing physical read-through molecules. This is not a maximum genomic distance between contact ends. Change only with knowledge of the physical insert distribution. [6] |
| Mapping fraction, single-sided fraction, usable yield | **No universal threshold** | Depend on genome, read length, assembly, quality, protocol, and MAPQ. Unexpected loss relative to matched controls suggests an assembly mismatch, contamination, short/poor reads, or complex/ambiguous alignments. | 

`parse.stats` retains pair-type information. Later filtering retains unique and rescued-unique `UU`, `UR`, and `RU` contacts; rejecting every non-`UU` pair would unnecessarily discard rescued Micro-C molecules. [7]

## Step 4 — Duplicate removal and library complexity

Duplicates are defined by both mapped endpoints and strands. The default is **exact endpoint matching (`DEDUP_MISMATCH=0`)**. This is a conservative starting choice for dense Micro-C, not a universal optimum. Pairtools defaults to a 3 bp tolerance; different published workflows use other tolerances. Broader tolerances remove near-identical alignments but can also merge independent molecules. Without UMIs, exact independent endpoint collisions cannot be reliably distinguished from duplicates.

The reported duplicate fraction is **duplicates / mapped pairs before deduplication**, not duplicates / all FASTQ pairs.

For a **5–10 million read-pair pilot from at least 5 million input cells**, the Micro-C-XL protocol describes **below approximately 5–10% duplicates** as desirable. A higher pilot fraction can indicate limited complexity, excessive amplification, low input, or sequencer clustering artifacts. This range is **not a failure cutoff for a deeply sequenced or low-input library**: duplication naturally increases as an existing molecule pool is sampled repeatedly. [5]

Inspect unique yield and `complexity_and_convergence.tsv` together. The complexity estimate assumes random sampling from a finite pool; unstable or infinite estimates in shallow data do not prove unlimited complexity. If original standard Illumina identifiers are available, `RUN_BYTILE=1` also generates tile-level duplicate statistics to help investigate clustering effects. It does not turn this into a UMI-aware workflow. [8]

## Step 5 — Nuclear contact QC, before distance filtering

`contact_qc.tsv` gives a count, denominator, and percentage for every metric. Cis, trans, and the cis distance bands use **unique, high-quality nuclear pairs before any minimum-distance cut**. The cis distance bands sum to the cis fraction; cis plus trans sums to 100%.

| Metric | Value to keep / benchmark | What a departure means and what to do |
|---|---|---|
| Cis fraction | **≥70% desirable** in the cited bulk Micro-C-XL protocol [5] | Below 70% can indicate excess random interchromosomal ligation, preparation problems, or different biology. High cis is not sufficient for passing QC: unligated local fragments also inflate it. Evaluate distance and orientation profiles together. |
| Trans fraction | Corresponding **≤30%** benchmark, using the same denominator | Elevated trans can reflect background ligation, but true chromosome rearrangements or organism-specific architecture can also contribute. Inspect the whole-genome map and matched controls. |
| Cis <1 kb, 1–10 kb, ≥10 kb | **No universal required percentage** | Micro-C genuinely enriches local contacts. Excess short cis plus strong strand bias suggests contamination; a reduced long-range fraction weakens distal analyses. Compare protocols and sample conditions at matched depth. |
| Strand orientation at sufficient distance | Four orientations approach approximately **25% each**, not a strict tolerance band [5,8] | Persistent bias at well-sampled longer distances suggests an artifact or processing problem. At very long, sparse distances, noisy proportions are expected. There is no universal distance at which convergence must occur. |
| Inward-facing (`+-`) enrichment near ~300–350 bp | **Diagnostic pattern, no universal fraction cutoff** [5] | Can indicate unligated/undigested dinucleosomes. Inspect the digestion trace and whether the excess persists across short distances. Short-range strand bias is not proof that every affected contact is false. |
| Contact decay, P(s) | Broad decay with distance; **no universal acceptable slope or plateau position** | Early flattening or unusual short-range peaks can indicate increased background or fragment artifacts; cell state and crosslinking also change the curve. Compare matched controls and orientation-specific curves. [8] |

`orientation_and_scaling.png` shows orientation proportions and raw distance scaling. The scaling calculation corrects distance-bin width and finite chromosome length, then normalizes by total cis contacts. It does **not** correct mappability or coverage biases and is not a balanced or observed/expected matrix. Do not treat its absolute values as a universal assay score.

## Step 6 — Optional short-distance exclusion

The default **`MIN_CIS_BP=0` retains all cis distances**. This preserves real nucleosome-scale contacts. A threshold such as **500 or 1,000 bp** may be appropriate for a specific downstream long-range analysis if strand profiles show local contamination, but it also deletes real short-range interactions. Choose it from the observed QC and scientific question, not as a compulsory Micro-C cleanup step.

`distance_filter_qc.tsv` records contacts retained and removed. All pre-filter nuclear pairs and their QC remain available. No restriction-fragment assignment, restriction-site motif filter, or conventional restriction-fragment “dangling-end” rule is applied: MNase lacks a fixed restriction-site map.

## Step 7 — Raw contact matrices

`cooler cload pairs` bins retained contacts, and `cooler zoomify` sums raw counts at coarser resolutions. Pairtools coordinates are **1-based**; the script correctly omits Cooler's `--zero-based` flag. [9]

**Hard QC:** every resolution must have exactly the same total stored contact count as `matrix_input.stats`. A mismatch indicates data loss or an implementation/coordinate problem; the script stops. Individual matrix pixels are counts of contacts between two bins, not the number of independently validated interactions.

## Step 8 — Balance each resolution

Balancing creates bin weights independently at every resolution. It never replaces or deletes the stored raw contact counts.

| Parameter / check | Script value | Why; meaning outside the retained range |
|---|---|---|
| Nonzero pixels per bin | **`--min-nnz 10`** | Cooler masks bins with too few nonzero entries after applicable filters. This improves numerical stability; ten entries are not sufficient evidence for a biological feature. |
| Low-coverage outliers | **`--mad-max 5`** | Filters very low marginals using a log-coverage median-absolute-deviation rule. Low bins can be gaps, unmappable regions, or genuinely undercovered sequence. This is not “five times lower than mean coverage.” |
| Near-diagonal exclusion during balancing | **`--ignore-diags 2`**, **`--ignore-dist 1000`** | Ignores the main/adjacent diagonals and applies the distance-derived exclusion when larger. This limits domination by local contacts. The effective distance depends on bin size; excluded pixels remain in the raw map. |
| Convergence tolerance | **`--tol 1e-5`** | Numerical convergence criterion, not an error rate or acceptable contamination percentage. |
| Iteration limit | **500** | Computational stopping limit, not a biological threshold. Failure to converge calls for reviewing coverage, filters, and structure; more iterations alone may not fix it. |
| Usable weights | **Converged, finite, positive weights** in the regions analyzed | Invalid/masked weights mark unusable bins for balanced analysis. There is no universal allowed percentage of masked bins. Compare across resolutions and inspect which chromosomes/regions are affected. |

These settings follow Cooler's documented options. [9] Optional `BLACKLIST_BED` masks matching bins **for balancing only**; it neither removes pairs nor changes the coverage-QC denominator. A failed normalization stores NaN weights rather than presenting unconverged weights as successful. A `.cool` with valid raw counts can still be created in that case; it is not ready for balanced analyses.

## Step 9 — Select a defensible resolution

The default outputs are **1, 2, 5, 10, 20, 50, and 100 kb**, with **10 kb** selected for the standalone file. These are requested bin sizes, not promises of supported biological resolution. Extremely fine maps can be built from shallow data and still contain little information.

1. Start with `resolution_qc.tsv`. A historical Hi-C heuristic is the finest resolution at which **at least 80% of bins have at least 1,000 contacts**. The script reports a practical adaptation using symmetric raw row sums, with diagonal counts counted once, over **all bins of the chosen chromosomes**, including empty/gap bins. This is a coverage heuristic, not a statistical guarantee of loop detection. [10,11]
2. Also inspect **cis, off-diagonal coverage**, because trans contacts and intense local/diagonal Micro-C counts can inflate all-contact coverage. No validated universal cutoff is asserted for this companion metric. `resolution_by_chromosome.tsv` reveals chromosomes hidden by a genome-wide average.
3. Require usable balancing weights and acceptable map appearance in the regions of interest. A global pass does not establish adequate coverage at every locus.
4. Match resolution to the question: compartments usually tolerate coarser maps; domains and loops require finer, deeper maps; nucleosome-level interpretation requires much stronger local support. Select the coarsest bin size that answers the question reliably, and confirm stability with neighboring resolutions and replicates.
5. Compare samples at the same resolution and assess comparable usable contact depth. Use depth-matched subsampling when evaluating replicate/map quality. Do not assign one universal bin size solely from total read count.

The 80%/1,000 heuristic originated in Hi-C. It is useful for consistency but does not establish a Micro-C-specific performance boundary. Falling below it means most bins lack that chosen coverage level, not that the entire experiment failed. Passing it does not guarantee reliable single-pixel calls, loops, or differential contacts. Do not discard empty bins simply to make the metric pass.

If the preferred resolution already exists in the `.mcool`, export it without realigning:

```bash
cooler cp sample1_microc/maps/sample1.mcool::/resolutions/20000 \
  sample1_microc/maps/sample1.20000.cool
```

For another set of resolutions, change `RESOLUTIONS` and `RESOLUTION` before running; every bin size must be an integer multiple of the finest one. Use a manageable finest grid for the genome and depth. Creating a finer grid never creates additional contacts.

## Step 10 — Map inspection, replicates, and release of results

Inspect `map_preview.png`, `whole_genome_preview.png`, and the manual checklist. The local preview defaults to a central region of the longest chromosome; specify `VIEW_REGION=chr1:10000000-12000000` for a relevant locus, keeping the region within 1,000 bins. The overview is a coarse raw map; it is not a copy-number-corrected map.

Check for empty stripes, unusually strong interchromosomal blocks, chromosome-wide coverage differences, local contamination, and whether balancing improves visibility without introducing misleading artifacts. Chromosome gaps and sex-chromosome dosage may be expected; CNVs/rearrangements may be real biology. The raw and balanced previews use independent logarithmic color scales: their colors cannot be compared numerically.

For biological replicates, process separately, compare digestion and QC metrics, then compare distance-stratified contact agreement and maps at matched depth and resolution. A method such as **HiCRep** accounts for contact-distance effects better than an unstratified genome-wide Pearson correlation. There is **no universal acceptable SCC/correlation cutoff** across depths, resolutions, and biological conditions; compare with relevant replicate and nonreplicate controls. [12] The supplied one-library script does not compute replicate concordance or automatically verify sample identity.

Only pool libraries after this review. Sum accepted libraries' **raw counts**, then rebalance the pooled matrix. Do not average balanced matrices as a substitute for pooling counts, and retain independent replicate results for downstream statistical testing.

## Output checklist

| Output under `OUTDIR` | What to inspect |
|---|---|
| `run_settings.sh`, `logs/` | Exact settings, versions, mapping/filtering/balancing logs, failures |
| `qc/fastq_integrity.tsv`, `qc/fastqc_*/`, `qc/fastqc_modules.tsv` | Read integrity and all FastQC module results |
| `qc/cutadapt.json` and trimmed FastQC, if enabled | Adapter trimming and read-loss effects |
| `qc/parse.stats`, `qc/dedup.stats`, optional `qc/duplicates_by_tile.tsv` | Mapping classes, duplicates, complexity |
| `qc/contact_qc.tsv` | Explicit-denominator yield, excluded chromosomes, cis/trans, distance bands |
| `qc/complexity_and_convergence.tsv` | Pairtools complexity and orientation convergence statistics |
| `qc/orientation_by_distance.tsv`, `qc/orientation_and_scaling.png` | Short-range contamination and raw P(s) |
| `qc/nuclear_before_distance_filter.stats`, `qc/matrix_input.stats`, `qc/distance_filter_qc.tsv` | Effect of minimum-distance filtering |
| `qc/resolution_qc.tsv`, `qc/resolution_by_chromosome.tsv` | Coverage, masked bins, convergence, count conservation |
| `qc/map_preview.png`, `qc/whole_genome_preview.png`, `qc/manual_review.tsv` | Visual checks and outstanding manual assessments |
| `qc/multiqc/`, if MultiQC succeeds | Convenient aggregate report; dedicated QC outputs remain authoritative |
| `pairs/` | Sorted, nonduplicate, nuclear, and final filtered contacts for reanalysis |
| `maps/SAMPLE.RESOLUTION.cool`, `maps/SAMPLE.mcool` | Final standalone and multiresolution matrices, raw counts plus weights |

The script intentionally retains intermediate pairs and the finest raw matrix. Archive what is needed for reproducibility before deleting large intermediates.

## Validation and limits

The workflow is configured for pairtools 1.1.2 and cooler 0.10.4. Its downstream contact-processing, filtering, balancing, QC, and export logic is checked with synthetic alignments containing known duplicates, short-range contacts, organelle contacts, low-MAPQ reads, and unmapped reads. FASTQ-integrity checks are tested separately. Synthetic checks establish software behavior, not assay performance. No user FASTQs or reference genome were supplied or processed; native BWA alignment, FastQC, optional Cutadapt, and wet-lab/replicate QC still require validation in the actual analysis environment.

## Sources

1. [FastQC: per-base sequence quality](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html).
2. [FastQC: adapter content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/10%20Adapter%20Content.html).
3. [FastQC: per-base N content](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html).
4. [Cutadapt user guide](https://cutadapt.readthedocs.io/en/stable/guide.html).
5. Metelova et al. (2023), [Micro-C-XL protocol, JoVE, DOI: 10.3791/64579](https://doi.org/10.3791/64579); [protocol PDF](https://api.jove.com/api/article/pdf/64579). The digestion, pilot duplication, and cis benchmarks above are contextualized to this protocol.
6. [Pairtools protocols and pipelines](https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html); [CLI reference](https://pairtools.readthedocs.io/en/latest/cli_tools.html).
7. [Pairtools formats and pair types](https://pairtools.readthedocs.io/en/latest/formats.html).
8. [Pairtools statistics documentation](https://pairtools.readthedocs.io/en/latest/stats.html); Open2C et al. (2024), [Pairtools: from sequencing data to chromosome contacts, DOI: 10.1371/journal.pcbi.1012164](https://doi.org/10.1371/journal.pcbi.1012164).
9. [Cooler CLI reference](https://cooler.readthedocs.io/en/latest/cli.html) and [API reference](https://cooler.readthedocs.io/en/latest/api.html).
10. Rao et al. (2014), [A 3D Map of the Human Genome at Kilobase Resolution, DOI: 10.1016/j.cell.2014.11.021](https://doi.org/10.1016/j.cell.2014.11.021).
11. [HiCRes: a tool to estimate and predict the genomic resolution of Hi-C libraries, DOI: 10.1093/nar/gkab1235](https://doi.org/10.1093/nar/gkab1235).
12. Yang et al. (2017), [HiCRep: assessing the reproducibility of Hi-C data using a stratum-adjusted correlation coefficient, DOI: 10.1101/gr.220640.117](https://doi.org/10.1101/gr.220640.117).
