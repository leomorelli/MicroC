#!/usr/bin/env bash
# Micro-C FASTQ -> contact pairs -> QC -> .mcool and a standalone .cool
# Version 2, 2026-10-08. Companion: MicroC_pipeline_QC_guide.md.
# For bulk, genome-wide, paired-end Micro-C without UMIs (Linux/Bash).
# Run each independently prepared library separately. Combine sequencing lanes
# from the SAME library before this workflow so duplicates across lanes are removed.
# Do not pool biological replicates before checking their reproducibility.
#
# SETUP (run once, outside this script; mamba can be replaced with conda):
# mamba create -n microc -c conda-forge -c bioconda --strict-channel-priority \
#   python=3.11 'numpy<2' 'pandas<3' pairtools=1.1.2 cooler=0.10.4 \
#   bwa samtools htslib fastqc multiqc matplotlib cutadapt
# conda activate microc
#
# RUN:
# bash microc_fastq_to_cool.sh /path/to/genome.fa R1.fastq.gz R2.fastq.gz
# RESOLUTION=5000 OUTDIR=sample2 SAMPLE=sample2 bash microc_fastq_to_cool.sh \
#   /path/to/genome.fa R1.fastq.gz R2.fastq.gz
#
# Only a reference FASTA is additionally required; indexes can be built below.
# Parameters can be changed here or supplied as environment variables.
# References (QC ranges are contextual benchmarks, NOT universal pass/fail rules):
# https://pairtools.readthedocs.io/en/latest/protocols_pipelines.html
# https://pairtools.readthedocs.io/en/latest/formats.html
# https://cooler.readthedocs.io/en/latest/cli.html
# https://doi.org/10.3791/64579                Micro-C-XL protocol and QC
# https://doi.org/10.1371/journal.pcbi.1012164  Orientation-dependent QC
# https://doi.org/10.1016/j.cell.2014.11.021    Coverage-based map resolution
# https://doi.org/10.1093/nar/gkab1235          HiCRes; resolution limitations

set -Eeuo pipefail
export LC_ALL=C
export MPLBACKEND=Agg
# Avoid hidden numerical-library thread pools in addition to pipeline threads.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
trap 'echo "ERROR at line $LINENO. Check logs; this run is incomplete." >&2' ERR

if [[ $# -lt 1 || $# -gt 3 || ${1:-} == --help ]]; then
    echo "Usage: bash $0 reference.fa [R1.fastq.gz] [R2.fastq.gz]"
    echo "Settings: SAMPLE, OUTDIR, THREADS, RESOLUTION, RESOLUTIONS, CHROMS_FILE, MIN_CIS_BP"
    [[ ${1:-} == --help ]] && exit 0
    exit 2
fi
REF=$(realpath "$1")
R1=$(realpath "${2:-R1.fastq.gz}")
R2=$(realpath "${3:-R2.fastq.gz}")
SAMPLE=${SAMPLE:-microc}
OUTDIR=${OUTDIR:-microc_output}
THREADS=${THREADS:-8}       # Per major tool; mapping/parsing/sorting overlap.
SORT_MEMORY=${SORT_MEMORY:-4G}  # GNU sort buffer, NOT total pipeline RAM.
RESOLUTION=${RESOLUTION:-10000} # Requested standalone .cool bin size, in bp.
RESOLUTIONS=${RESOLUTIONS:-1000,2000,5000,10000,20000,50000,100000}
MAPQ=${MAPQ:-30}           # BWA MAPQ >=30 at BOTH ends: stringent mapping filter.
WALKS_POLICY=${WALKS_POLICY:-5unique}
# 5unique reports one contact per molecule; complex walks may be indirect.
# Use mask for a more conservative policy. Keep this identical across samples.
MAX_MOLECULE_SIZE=${MAX_MOLECULE_SIZE:-750}
# 750 bp is pairtools' rescue default, NOT a genomic contact-distance cutoff.
# Adjust to the physical insert-size distribution if your protocol requires it.
DEDUP_MISMATCH=${DEDUP_MISMATCH:-0}
# 0 = exact endpoint/strand duplicates; an explicit conservative starting choice
# for dense Micro-C data, NOT an established universal optimum. pairtools defaults
# to 3 bp; some published pipelines use 1 bp. Compare 0/1/3 in a pilot if needed.
# Without UMIs, independent molecules sharing endpoints cannot be distinguished.
MIN_CIS_BP=${MIN_CIS_BP:-0}
# 0 retains ALL genomic distances in the output map. This is intentional:
# inspect short-range orientation QC BEFORE imposing e.g. 500 or 1000 bp.
# Setting 1000 removes real local nucleosome contacts as well as contaminants.
CHROMS_FILE=${CHROMS_FILE:-}
# Optional text file: one analysis chromosome name per line, matching the FASTA.
# For human/mouse whole-genome analyses, a PRIMARY nuclear chromosome list is
# recommended, to exclude alternative haplotypes/decoys from matrix axes and QC.
# Alignment still uses the full supplied reference. No organism is assumed here.
EXCLUDE_CHROMS=${EXCLUDE_CHROMS:-chrM,MT,M}
# Check names for your organism; add chloroplast/other organelle contigs if needed.
BLACKLIST_BED=${BLACKLIST_BED:-} # Optional same-assembly BED, for balancing only.
BALANCE_IGNORE_BP=${BALANCE_IGNORE_BP:-1000}
# At least the main and immediately adjacent diagonals are ignored in balancing.
# This DOES NOT delete those contacts from the raw .cool/.mcool.
BWA_INDEX=${BWA_INDEX:-$REF}    # Existing BWA-MEM index prefix, if available.
TRIM_ADAPTERS=${TRIM_ADAPTERS:-0} # 1 enables paired Cutadapt + post-trim FastQC.
ADAPTER_R1=${ADAPTER_R1:-}       # Actual library 3' sequencing adapter, when trimming.
ADAPTER_R2=${ADAPTER_R2:-}       # Supply from the kit; do NOT trim genomic ligation motifs.
TRIM_QUALITY=${TRIM_QUALITY:-20} # Tail trimming Q20 is nominal 1% base-error probability.
MIN_READ_LENGTH=${MIN_READ_LENGTH:-30} # Practical minimum AFTER trimming, not a QC standard.
RUN_BYTILE=${RUN_BYTILE:-0}     # 1 only with original standard Illumina read identifiers.
VIEW_REGION=${VIEW_REGION:-}   # Optional, e.g. chr1:10000000-12000000; <=1000 bins.

# STEP 0: Check settings/dependencies and create a NEW output directory.
# Refuse an existing output directory to avoid overwrites or appending old stats.
[[ $SAMPLE =~ ^[A-Za-z0-9_.-]+$ ]] || { echo 'Invalid SAMPLE name' >&2; exit 2; }
for value in "$THREADS" "$RESOLUTION" "$MAPQ" "$MAX_MOLECULE_SIZE" "$DEDUP_MISMATCH" "$MIN_CIS_BP" "$BALANCE_IGNORE_BP" "$TRIM_QUALITY" "$MIN_READ_LENGTH"; do
    [[ $value =~ ^[0-9]+$ ]] || { echo "Expected nonnegative integer: $value" >&2; exit 2; }
done
(( THREADS > 0 && RESOLUTION > 0 && MAX_MOLECULE_SIZE > 0 ))
[[ $TRIM_ADAPTERS =~ ^[01]$ && $RUN_BYTILE =~ ^[01]$ ]] || { echo 'Boolean settings must be 0 or 1' >&2; exit 2; }
for file in "$CHROMS_FILE" "$BLACKLIST_BED"; do
    [[ -z $file || -s $file ]] || { echo "Missing/empty optional file: $file" >&2; exit 2; }
done
for file in "$REF" "$R1" "$R2"; do [[ -s $file ]] || { echo "Missing/empty: $file" >&2; exit 2; }; done
[[ $R1 != "$R2" ]] || { echo 'R1 and R2 must be different files' >&2; exit 2; }
for program in bwa samtools bgzip fastqc pairtools cooler python; do
    command -v "$program" >/dev/null || { echo "Missing program: $program" >&2; exit 2; }
done
python -c 'import cooler, numpy, pandas, matplotlib, h5py'
if (( TRIM_ADAPTERS )); then
    command -v cutadapt >/dev/null
    [[ -n $ADAPTER_R1 && -n $ADAPTER_R2 ]] || { echo 'Supply ADAPTER_R1 and ADAPTER_R2 for trimming' >&2; exit 2; }
fi
RESOLUTIONS=$(python - "$RESOLUTIONS" "$RESOLUTION" <<'PY'
import sys
r = sorted(set(map(int, sys.argv[1].split(','))) | {int(sys.argv[2])})
if r[0] <= 0 or any(x % r[0] for x in r):
    raise SystemExit('All resolutions must be positive integer multiples of the finest bin size.')
print(','.join(map(str, r)))
PY
)
BASE_RES=${RESOLUTIONS%%,*}
mkdir "$OUTDIR"
OUTDIR=$(realpath "$OUTDIR")
mkdir -p "$OUTDIR"/{reference,pairs,maps,reads,qc/fastqc_R1,qc/fastqc_R2,logs,tmp}
exec > >(tee "$OUTDIR/logs/pipeline.log") 2>&1
printf 'Started %s\n' "$(date -Is)"
declare -p REF R1 R2 SAMPLE THREADS SORT_MEMORY RESOLUTION RESOLUTIONS MAPQ \
    WALKS_POLICY MAX_MOLECULE_SIZE DEDUP_MISMATCH MIN_CIS_BP CHROMS_FILE \
    EXCLUDE_CHROMS BLACKLIST_BED BALANCE_IGNORE_BP BWA_INDEX TRIM_ADAPTERS \
    ADAPTER_R1 ADAPTER_R2 TRIM_QUALITY MIN_READ_LENGTH RUN_BYTILE VIEW_REGION > "$OUTDIR/run_settings.sh"
{ pairtools --version; cooler --version; fastqc --version; samtools --version; python --version; } > "$OUTDIR/logs/versions.txt"
# BWA writes its version in the alignment log below. Record the full environment.
python -m pip freeze > "$OUTDIR/logs/python_packages.txt"

# STEP 1: Raw-read QC. Both FASTQs must have corresponding reads in the same order.
# Inspect FastQC HTML: base qualities, adapters, read lengths, overrepresented reads.
# Practical sequencing benchmark: most bases Q30 or above is desirable; this is
# not a Micro-C-specific hard cutoff. FastQC duplication estimates are NOT the
# coordinate-based PCR duplicate metric used below. There is no universal
# acceptable adapter percentage; substantial read-through motivates trimming.
# Trimming is off by default: BWA local alignment handles clipped/chimeric reads.
# If adapter read-through justifies trimming, set TRIM_ADAPTERS=1 and supply both
# adapter sequences; paired trimming and post-trim FastQC then run below.
# Do not require a restriction-site motif or trim genomic ligation junctions.
echo 'STEP 1: FASTQ quality control'
# Hard integrity requirements: zero malformed records, zero mate-ID mismatches,
# equal nonzero read counts. Stream the entire gzip inputs to detect truncation.
# This extra pass is deliberate and adds I/O time on large libraries.
python - "$R1" "$R2" "$OUTDIR/qc/fastq_integrity.tsv" <<'PY'
import csv, gzip, sys
from collections import Counter
from contextlib import ExitStack
hist = [Counter(), Counter()]; lengths = [Counter(), Counter()]; bases = [0,0]; ns = [0,0]
def record(f):
    h = f.readline()
    if not h: return None
    s, p, q = (f.readline().rstrip(b'\r\n') for _ in range(3))
    if not h.startswith(b'@') or not p.startswith(b'+') or not s or len(s)!=len(q):
        raise SystemExit('Malformed/empty FASTQ record or sequence/quality length mismatch.')
    if min(q)<33 or max(q)>126: raise SystemExit('Invalid Phred+33 quality character.')
    name=h[1:].split()[0]
    if name.endswith((b'/1',b'/2')): name=name[:-2]
    return name,s,q
count=0
with ExitStack() as stack:
    f1=stack.enter_context(gzip.open(sys.argv[1],'rb'))
    f2=stack.enter_context(gzip.open(sys.argv[2],'rb'))
    while True:
        a,b=record(f1),record(f2)
        if a is None and b is None: break
        if a is None or b is None: raise SystemExit('R1/R2 read counts differ.')
        if a[0]!=b[0]: raise SystemExit(f'R1/R2 IDs differ at pair {count+1}.')
        for j,(_,s,q) in enumerate((a,b)):
            hist[j].update(q); lengths[j][len(s)]+=1; bases[j]+=len(s); ns[j]+=s.upper().count(b'N')
        count+=1
if count==0: raise SystemExit('FASTQ files have no reads.')
with open(sys.argv[3],'w') as f:
    w=csv.writer(f,delimiter='\t')
    w.writerow(['mate','reads','bases','min_length','max_length','pct_bases_Q20','pct_bases_Q30','pct_N','integrity'])
    for j in range(2):
        w.writerow([f'R{j+1}',count,bases[j],min(lengths[j]),max(lengths[j]),
            100*sum(n for q,n in hist[j].items() if q>=53)/bases[j],
            100*sum(n for q,n in hist[j].items() if q>=63)/bases[j],100*ns[j]/bases[j],'PASS'])
PY
fastqc --threads 1 --outdir "$OUTDIR/qc/fastqc_R1" "$R1"
fastqc --threads 1 --outdir "$OUTDIR/qc/fastqc_R2" "$R2"
# FastQC defaults: per-base quality WARN if median<25 or lower quartile<10;
# FAIL if median<20 or lower quartile<5. Adapter content WARN>5%, FAIL>10%.
# Per-base N content WARN>5%, FAIL>20%. These trigger review, not sample rejection.
if (( TRIM_ADAPTERS )); then
    mkdir -p "$OUTDIR/qc/fastqc_trimmed_R1" "$OUTDIR/qc/fastqc_trimmed_R2"
    cutadapt -j "$THREADS" -a "$ADAPTER_R1" -A "$ADAPTER_R2" \
        -q "$TRIM_QUALITY" --minimum-length "$MIN_READ_LENGTH" --pair-filter=any \
        --json "$OUTDIR/qc/cutadapt.json" \
        -o "$OUTDIR/reads/trimmed.R1.fastq.gz" -p "$OUTDIR/reads/trimmed.R2.fastq.gz" \
        "$R1" "$R2" > "$OUTDIR/logs/cutadapt.log"
    R1="$OUTDIR/reads/trimmed.R1.fastq.gz"; R2="$OUTDIR/reads/trimmed.R2.fastq.gz"
    fastqc --threads 1 --outdir "$OUTDIR/qc/fastqc_trimmed_R1" "$R1"
    fastqc --threads 1 --outdir "$OUTDIR/qc/fastqc_trimmed_R2" "$R2"
fi
# Retain all native FastQC statuses in an easily comparable table.
python - "$OUTDIR/qc" <<'PY'
import csv,sys,zipfile
from pathlib import Path
qc=Path(sys.argv[1])
with (qc/'fastqc_modules.tsv').open('w') as out:
    w=csv.writer(out,delimiter='\t'); w.writerow(['report','status','module','filename'])
    for path in sorted(qc.glob('fastqc_*/*.zip')):
        with zipfile.ZipFile(path) as z:
            key=next(k for k in z.namelist() if k.endswith('/summary.txt'))
            for line in z.read(key).decode().splitlines(): w.writerow([str(path.relative_to(qc)),*line.split('\t')])
PY

# WET-LAB QC (not recoverable conclusively from FASTQ): inspect pre-ligation MNase
# digestion and final library traces. Several academic bulk protocols target
# ~70-90% mononucleosomes with little subnucleosomal material; other kits differ.
# Mononucleosomal DNA is ~150 bp BEFORE ligation, not the final library insert.
# Excess di-/oligonucleosomes suggests underdigestion; a substantial short smear
# suggests overdigestion. Compare digestion levels across biological conditions.

# STEP 2: Prepare reference lengths and alignment index, keeping input files intact.
# Use an uncompressed FASTA; names/sizes are derived from the SAME assembly.
echo 'STEP 2: Reference preparation'
ln -s "$REF" "$OUTDIR/reference/genome.fa"
samtools faidx "$OUTDIR/reference/genome.fa"
cut -f1,2 "$OUTDIR/reference/genome.fa.fai" > "$OUTDIR/reference/all.chrom.sizes"
python - "$OUTDIR/reference/all.chrom.sizes" "$CHROMS_FILE" "$EXCLUDE_CHROMS" \
    "$OUTDIR/reference/analysis.chrom.sizes" <<'PY'
import sys
from pathlib import Path
src, keepfile, excluded, dst = sys.argv[1:]
sizes = dict(line.split()[:2] for line in Path(src).read_text().splitlines())
keep = set(sizes)
if keepfile:
    keep = {l.split()[0] for l in Path(keepfile).read_text().splitlines() if l.strip() and not l.startswith('#')}
    if keep - sizes.keys():
        raise SystemExit(f'Chromosomes absent from FASTA: {sorted(keep - sizes.keys())}')
keep -= set(excluded.split(','))
if not keep:
    raise SystemExit('No analysis chromosomes remain.')
Path(dst).write_text(''.join(f'{c}\t{n}\n' for c, n in sizes.items() if c in keep))
PY
index_complete=1
for suffix in amb ann bwt pac sa; do [[ -s ${BWA_INDEX}.${suffix} ]] || index_complete=0; done
if (( ! index_complete )); then
    BWA_INDEX="$OUTDIR/reference/bwa"
    bwa index -p "$BWA_INDEX" "$REF" 2> "$OUTDIR/logs/bwa_index.log"
fi
printf 'Effective BWA index: %s\n' "$BWA_INDEX"

# STEP 3: Align -> reconstruct contact pairs -> sort for deduplication.
# -SP disables ordinary paired-end pairing/mate rescue, inappropriate for Hi-C.
# Preserve supplementary alignments; pairtools interprets chimeric molecules.
# MAPQ >=30 defines unique alignments here. No universal mapping-rate percentage:
# assess relative to matched libraries, read length, genome, and this MAPQ cutoff.
# 5unique yields one pair per sequenced molecule, enabling interpretable counts.
# --drop-sam saves space; retained .pairs include read IDs, strands and MAPQs.
# Sorting follows parsing, not coordinate sorting of BAM (which splits read IDs).
echo 'STEP 3: BWA-MEM, contact parsing and sorting'
bwa mem -t "$THREADS" -SP "$BWA_INDEX" "$R1" "$R2" 2> "$OUTDIR/logs/bwa_mem.log" | \
    pairtools parse --chroms-path "$OUTDIR/reference/all.chrom.sizes" \
        --min-mapq "$MAPQ" --walks-policy "$WALKS_POLICY" \
        --max-molecule-size "$MAX_MOLECULE_SIZE" --add-columns mapq --drop-sam \
        --output-stats "$OUTDIR/qc/parse.stats" 2> "$OUTDIR/logs/parse.log" | \
    pairtools sort --nproc "$THREADS" --memory "$SORT_MEMORY" \
        --tmpdir "$OUTDIR/tmp" --nproc-in 1 --nproc-out 2 \
        --output "$OUTDIR/pairs/${SAMPLE}.sorted.pairs.gz" 2> "$OUTDIR/logs/sort.log"

# STEP 4: Remove PCR/optical duplicates using BOTH endpoints and strands.
# The duplicate fraction denominator is mapped pairs BEFORE deduplication.
# Published example: <~5-10% duplication in a 5-10 million read-pair pilot from
# >=5 million input cells. This is NOT a limit for deep or low-input libraries.
# Sequencer clustering duplicates can inflate pilot estimates; inspect library
# complexity and unique yield, rather than applying a fixed duplicate cutoff.
echo 'STEP 4: Deduplication'
dedup_extra=()
if (( RUN_BYTILE )); then
    dedup_extra+=(--output-bytile-stats "$OUTDIR/qc/duplicates_by_tile.tsv")
fi
pairtools dedup --max-mismatch "$DEDUP_MISMATCH" --nproc-in 1 --nproc-out 2 \
    "${dedup_extra[@]}" \
    --output "$OUTDIR/pairs/${SAMPLE}.nodup.pairs.gz" \
    --output-stats "$OUTDIR/qc/dedup.stats" \
    "$OUTDIR/pairs/${SAMPLE}.sorted.pairs.gz" 2> "$OUTDIR/logs/dedup.log"

# STEP 5: Retain high-quality contacts on the requested nuclear chromosomes.
# Keep rescued unique pairs (UR/RU) as well as UU; otherwise short Micro-C
# fragments spanning a ligation junction can be lost unnecessarily.
# These are mapping-filtered pairs, not proof that every short contact is genuine.
# QC is computed BEFORE any optional genomic-distance exclusion.
echo 'STEP 5: Mapping/chromosome filtering and Micro-C QC'
pairtools select --chrom-subset "$OUTDIR/reference/analysis.chrom.sizes" \
    --nproc-in 1 --nproc-out 2 \
    --output "$OUTDIR/pairs/${SAMPLE}.nuclear.pairs.gz" \
    "((pair_type == 'UU') or (pair_type == 'UR') or (pair_type == 'RU')) and (mapq1 >= $MAPQ) and (mapq2 >= $MAPQ)" \
    "$OUTDIR/pairs/${SAMPLE}.nodup.pairs.gz"
pairtools stats --no-yaml --output "$OUTDIR/qc/nuclear_before_distance_filter.stats" \
    "$OUTDIR/pairs/${SAMPLE}.nuclear.pairs.gz"

# Short-range QC: an excess of inward (+-, FR) pairs around ~300-350 bp may
# reflect undigested dinucleosomes. Do NOT require 25% per orientation at every
# short distance; at sufficiently long distances, the four curves approach ~25%.
# There is no universal convergence distance. Use the observed profiles.
# For mammalian bulk Micro-C, >=70% cis is a published desirable benchmark,
# NOT a universal cutoff; unligated fragments also inflate the cis fraction.
# Report cis <1kb, 1-10kb, >=10kb and trans with a common explicit denominator.
python - "$OUTDIR/qc" "$OUTDIR/reference/analysis.chrom.sizes" "$EXCLUDE_CHROMS" <<'PY'
import sys, re, csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
qc, chromfile = Path(sys.argv[1]), Path(sys.argv[2])
excluded = set(sys.argv[3].split(','))
def read_stats(path):
    out = {}
    for line in path.read_text().splitlines():
        fields = line.split('\t')
        if len(fields) == 2:
            try: out[fields[0]] = float(fields[1])
            except ValueError: pass
    return out
d = read_stats(qc / 'dedup.stats')
n = read_stats(qc / 'nuclear_before_distance_filter.stats')
total = n.get('total', 0)
if total == 0:
    raise SystemExit('No high-quality nuclear contacts remain; inspect reference, mapping and chromosome names.')
rows = [
    ('both_unmapped_or_below_MAPQ', d.get('total_unmapped', 0), d.get('total', 0), 'all pairs entering alignment'),
    ('one_end_mapped', d.get('total_single_sided_mapped', 0), d.get('total', 0), 'all pairs entering alignment'),
    ('both_mapped_before_dedup', d.get('total_mapped', 0), d.get('total', 0), 'all parsed read pairs'),
    ('coordinate_duplicates', d.get('total_dups', 0), d.get('total_mapped', 0), 'mapped pairs before deduplication'),
    ('nuclear_unique_yield', total, d.get('total', 0), 'all parsed read pairs'),
    ('cis', n.get('cis', 0), total, 'unique high-quality nuclear pairs'),
    ('trans', n.get('trans', 0), total, 'unique high-quality nuclear pairs'),
    ('cis_lt_1kb', n.get('cis', 0)-n.get('cis_1kb+', 0), total, 'unique high-quality nuclear pairs'),
    ('cis_1_to_10kb', n.get('cis_1kb+', 0)-n.get('cis_10kb+', 0), total, 'unique high-quality nuclear pairs'),
    ('cis_ge_10kb', n.get('cis_10kb+', 0), total, 'unique high-quality nuclear pairs'),
]
# Report organelle/explicitly excluded chromosomes BEFORE their removal.
excluded_pairs=0
for key,value in d.items():
    fields=key.split('/')
    if len(fields)==3 and fields[0]=='chrom_freq' and (fields[1] in excluded or fields[2] in excluded):
        excluded_pairs+=value
rows.append(('excluded_chromosome_pairs',excluded_pairs,d.get('total_nodups',0),'all mapped nonduplicate pairs before chromosome selection'))
integrity=qc/'fastq_integrity.tsv'
if integrity.exists():
    raw_count=int(next(csv.DictReader(integrity.open(),delimiter='\t'))['reads'])
    rows.extend([('pairs_entering_alignment',d.get('total',0),raw_count,'raw FASTQ pairs'),
                 ('nuclear_yield_vs_raw',total,raw_count,'raw FASTQ pairs')])
with (qc / 'contact_qc.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['metric','count','denominator_count','percent','denominator'])
    for name, count, denom, label in rows:
        w.writerow([name, int(count), int(denom), 100*count/denom if denom else 'NA', label])
with (qc/'complexity_and_convergence.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['source','metric','value'])
    for label,data in [('dedup',d),('nuclear_before_distance_filter',n)]:
        for key,value in data.items():
            if 'complexity' in key or 'dist_freq_convergence' in key: w.writerow([label,key,value])
# A naive complexity estimate assumes random sampling from a finite molecule pool.
# Infinite/unstable estimates in shallow data are NOT evidence of unlimited yield.
# Convert logarithmic distance-bin counts to raw opportunity-normalized scaling.
# Counts / sum_c sum_{s=lo}^{hi-1}(L_c-s) / total_cis. This corrects bin width
# and finite chromosome length, but NOT mappability/coverage biases; not ICE/OE.
lengths = [int(l.split()[1]) for l in chromfile.read_text().splitlines()]
groups = {}
for key, count in n.items():
    match = re.fullmatch(r'dist_freq/(\d+)-(\d+)/([+-]{2})', key)
    if match:
        a, b, orientation = match.groups()
        groups.setdefault((int(a), int(b)), {})[orientation] = count
orientations = ['+-','-+','++','--']
records = []
for (a,b), counts in sorted(groups.items()):
    if a < 1: continue  # Do not interpret coincident endpoints as genomic scaling.
    k = np.array([counts.get(o, 0) for o in orientations], dtype=float)
    if k.sum() == 0: continue
    opportunities = 0.0
    for L in lengths:
        hi = min(b, L)
        if hi > a:
            opportunities += (hi-a)*L - (a+hi-1)*(hi-a)/2
    if opportunities <= 0: continue
    records.append((a,b,k,k/k.sum(),k/(opportunities*max(n.get('cis',0),1))))
with (qc / 'orientation_by_distance.tsv').open('w') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(['start_bp','end_bp']+[f'count_{o}' for o in orientations]+[f'fraction_{o}' for o in orientations]+[f'raw_scaling_{o}' for o in orientations])
    for a,b,k,frac,ps in records: w.writerow([a,b,*k,*frac,*ps])
if records:
    x = np.array([(a*b)**0.5 for a,b,*_ in records])
    fig, axes = plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for j,o in enumerate(orientations):
        axes[0].semilogx(x,[r[3][j] for r in records],label=o)
        y = np.array([r[4][j] for r in records]); y[y == 0] = np.nan
        axes[1].loglog(x,y,label=o)
    axes[0].axhline(0.25,color='grey',ls='--',lw=1)
    axes[0].set(xlabel='Genomic separation (bp)',ylabel='Fraction within distance bin',ylim=(0,1))
    axes[1].set(xlabel='Genomic separation (bp)',ylabel='Relative raw P(s)')
    for ax in axes: ax.legend(title='Ordered strands'); ax.spines[['top','right']].set_visible(False)
    fig.savefig(qc / 'orientation_and_scaling.png',dpi=180); plt.close(fig)
PY

# STEP 6: Apply an OPTIONAL minimum cis separation to the matrix input.
# No restriction-fragment assignment/filter is used: MNase has no fixed site map.
# Always preserve the nuclear pairs and pre-filter QC for revisiting this choice.
pairtools select --nproc-in 1 --nproc-out 2 \
    --output "$OUTDIR/pairs/${SAMPLE}.filtered.pairs.gz" \
    "(chrom1 != chrom2) or (abs(pos1-pos2) >= $MIN_CIS_BP)" \
    "$OUTDIR/pairs/${SAMPLE}.nuclear.pairs.gz"
pairtools stats --no-yaml --output "$OUTDIR/qc/matrix_input.stats" \
    "$OUTDIR/pairs/${SAMPLE}.filtered.pairs.gz"
python - "$OUTDIR/qc" <<'PY'
import csv,sys
from pathlib import Path
qc=Path(sys.argv[1])
def count(path):
    for line in path.read_text().splitlines():
        if line.startswith('total\t'): return int(line.split('\t')[1])
    raise SystemExit('Missing total in pairtools statistics.')
before=count(qc/'nuclear_before_distance_filter.stats'); after=count(qc/'matrix_input.stats')
if after==0: raise SystemExit('No contacts remain after distance filtering. Reduce MIN_CIS_BP or inspect mapping.')
with (qc/'distance_filter_qc.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['before','retained','removed','pct_retained'])
    w.writerow([before,after,before-after,100*after/before])
PY

# STEP 7: Count contacts in fixed genomic bins and construct multiple resolutions.
# pairtools .pairs positions are 1-based; do NOT use cooler's --zero-based flag.
# Columns 2/3 and 4/5 contain the two chromosome/position pairs.
# All coarser maps are sums of RAW counts; never coarsen balanced pixel values.
echo 'STEP 7: Contact matrices'
RAW_COOL="$OUTDIR/maps/${SAMPLE}.${BASE_RES}.raw.cool"
MCOOL="$OUTDIR/maps/${SAMPLE}.mcool"
cooler cload pairs --temp-dir "$OUTDIR/tmp" -c1 2 -p1 3 -c2 4 -p2 5 \
    "$OUTDIR/reference/analysis.chrom.sizes:$BASE_RES" \
    "$OUTDIR/pairs/${SAMPLE}.filtered.pairs.gz" "$RAW_COOL"
cooler zoomify --nproc "$THREADS" --resolutions "$RESOLUTIONS" \
    --out "$MCOOL" "$RAW_COOL"

# STEP 8: Balance EACH resolution independently. Raw counts remain unchanged.
# --min-nnz 10: require >=10 nonzero matrix entries per bin (algorithm setting,
# NOT a scientific coverage threshold). --mad-max 5: low-coverage outlier filter.
# tol 1e-5: algorithm convergence tolerance, not a data-quality percentage.
# No universal acceptable fraction of masked bins; compare primary chromosomes
# across replicates/resolutions and investigate excessive masking.
# store_nan avoids presenting unconverged weights as a successful normalization.
# Balance failure may reflect sparse data or unusual genomes, not just bin size.
echo 'STEP 8: Matrix balancing'
balance_extra=()
if [[ -n $BLACKLIST_BED ]]; then balance_extra+=(--blacklist "$BLACKLIST_BED"); fi
IFS=',' read -r -a resolution_array <<< "$RESOLUTIONS"
for r in "${resolution_array[@]}"; do
    cooler balance --nproc "$THREADS" --ignore-diags 2 --ignore-dist "$BALANCE_IGNORE_BP" \
        --min-nnz 10 --mad-max 5 --tol 1e-5 --max-iters 500 \
        --convergence-policy store_nan "${balance_extra[@]}" \
        "$MCOOL::/resolutions/$r" 2> "$OUTDIR/logs/balance_${r}.log"
done

# STEP 9: Assess coverage at each resolution WITHOUT loading a dense matrix.
# Historical heuristic (Rao 2014): >=80% of bins with >=1000 contacts.
# Here pct_bins_ge1000_all uses raw symmetric row sums (diagonal counted ONCE),
# over ALL bins of the chosen chromosome set, including zero-coverage/gap bins.
# This is a practical adaptation, not a formal resolution estimate or guarantee
# of loop detection. The second measure excludes trans AND diagonal pixels:
# short-range/diagonal Micro-C enrichment can inflate the all-contact criterion.
# Do not discard empty bins to make the metric pass. Specify primary chromosomes.
# A blacklist masks balancing weights but does not alter these coverage denominators.
echo 'STEP 9: Coverage and resolution assessment'
python - "$MCOOL" "$OUTDIR/qc/resolution_qc.tsv" "$RESOLUTION" "$OUTDIR/qc/matrix_input.stats" <<'PY'
import sys, csv
import numpy as np
import cooler
mcool, output, selected, statsfile = sys.argv[1:]
expected=int(next(line.split('\t')[1] for line in open(statsfile) if line.startswith('total\t')))
perchrom=[]
with open(output,'w') as f:
    w = csv.writer(f,delimiter='\t')
    w.writerow(['bin_bp','n_bins','n_pairs','median_contacts_all','pct_bins_ge1000_all',
                'median_contacts_cis_offdiag','pct_bins_ge1000_cis_offdiag',
                'pct_bins_with_finite_weight','balance_converged'])
    for group in sorted(cooler.fileops.list_coolers(mcool),key=lambda g:int(g.split('/')[-1])):
        c = cooler.Cooler(f'{mcool}::{group}'); bins = c.bins()[:]
        chrom = bins['chrom'].cat.codes.to_numpy()
        n = len(bins); total = np.zeros(n); cis = np.zeros(n)
        for start in range(0,c.info['nnz'],1_000_000):
            p = c.pixels()[start:start+1_000_000]
            i=p['bin1_id'].to_numpy(); j=p['bin2_id'].to_numpy(); v=p['count'].to_numpy()
            np.add.at(total,i,v)
            other=i!=j; np.add.at(total,j[other],v[other])
            keep=other & (chrom[i]==chrom[j])
            np.add.at(cis,i[keep],v[keep]); np.add.at(cis,j[keep],v[keep])
        assert c.info['sum']==expected, f'Contact-count mismatch at {c.binsize} bp'
        weights=bins['weight'].to_numpy()
        finite = np.isfinite(weights) & (weights>0)
        with c.open('r') as root:
            converged=bool(np.all(root['bins/weight'].attrs.get('converged',False)))
        w.writerow([c.binsize,n,c.info['sum'],np.median(total),100*np.mean(total>=1000),
                    np.median(cis),100*np.mean(cis>=1000),100*np.mean(finite),converged])
        for name in c.chromnames:
            mask=bins['chrom'].eq(name).to_numpy()
            perchrom.append([c.binsize,name,int(mask.sum()),float(np.median(total[mask])),
                100*float(np.mean(total[mask]>=1000)),100*float(np.mean(finite[mask]))])
        if str(c.binsize)==selected and (not converged or not finite.any()):
            print(f'WARNING: {selected}-bp normalization failed or has no usable weights. '
                  'Use raw counts only, inspect QC, and/or choose coarser bins.',file=sys.stderr)
from pathlib import Path
with (Path(output).parent/'resolution_by_chromosome.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['bin_bp','chrom','n_bins','median_contacts','pct_bins_ge1000','pct_bins_with_finite_weight'])
    w.writerows(perchrom)
PY
# Export the requested resolution as a standalone .cool, including its weights.
cooler cp "$MCOOL::/resolutions/$RESOLUTION" "$OUTDIR/maps/${SAMPLE}.${RESOLUTION}.cool"

# STEP 10: Inspect maps at two scales and export a manual review checklist.
# No fixed loop/TAD count is a universal QC threshold: cell cycle/perturbations
# can alter real structure. Previews use logarithmic, independent color scales;
# use the same limits and matched contact depths for quantitative comparisons.
echo 'STEP 10: Map previews and final review checklist'
python - "$OUTDIR" "$SAMPLE" "$RESOLUTION" "$VIEW_REGION" <<'PY'
import csv,sys
from pathlib import Path
import numpy as np
import cooler
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
out,sample,res,region=Path(sys.argv[1]),sys.argv[2],int(sys.argv[3]),sys.argv[4]
c=cooler.Cooler(str(out/'maps'/f'{sample}.{res}.cool'))
chrom=c.chromsizes.idxmax()
if not region:
    length=int(c.chromsizes[chrom]); width=min(length,5_000_000,res*1000)
    start=max(0,(length-width)//2); region=f'{chrom}:{start}-{start+width}'
lo,hi=c.extent(region)
if hi-lo>1000: raise SystemExit('VIEW_REGION exceeds 1000 bins. Choose a smaller region.')
fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
for ax,balanced in zip(axes,[False,True]):
    a=np.asarray(c.matrix(balance=balanced).fetch(region),dtype=float)
    positive=a[np.isfinite(a)&(a>0)]
    if len(positive):
        vmin=max(float(np.percentile(positive,5)),np.finfo(float).tiny)
        vmax=max(float(np.percentile(positive,99)),vmin*1.01)
        im=ax.imshow(np.ma.masked_invalid(np.where(a>0,a,np.nan)),origin='lower',cmap='Reds',norm=LogNorm(vmin=vmin,vmax=vmax))
        fig.colorbar(im,ax=ax,shrink=0.7)
    else: ax.text(0.5,0.5,'No usable values',ha='center',transform=ax.transAxes)
    ax.set_title('Balanced' if balanced else 'Raw counts'); ax.set_xlabel('Bin'); ax.set_ylabel('Bin')
fig.suptitle(region); fig.savefig(out/'qc/map_preview.png',dpi=170); plt.close(fig)
# Whole-genome counts aggregated into ~400 plus chromosome-boundary bins,
# WITHOUT a dense fine-resolution matrix. Many tiny contigs can exceed this.
factor=max(1,int(np.ceil(c.info['nbins']/400)))
coarse=out/'qc'/f'{sample}.overview.cool'
if factor>1:
    cooler.coarsen_cooler(c.uri,str(coarse),factor=factor,chunksize=1_000_000)
    overview=cooler.Cooler(str(coarse))
else:
    overview=c
# Refuse unexpectedly large dense previews; the contact maps are already saved.
if overview.info['nbins']>3000:
    raise SystemExit('Too many contigs for a genome preview; use CHROMS_FILE with primary chromosomes. Maps are saved.')
a=overview.matrix(balance=False)[:].astype(float)
fig,ax=plt.subplots(figsize=(6,6),layout='constrained')
v=a[a>0]
if len(v):
    vmin=max(1,float(np.percentile(v,5))); vmax=max(vmin*1.01,float(np.percentile(v,99)))
    im=ax.imshow(np.ma.masked_where(a<=0,a),cmap='Reds',norm=LogNorm(vmin=vmin,vmax=vmax))
    fig.colorbar(im,ax=ax,label='Raw counts')
for name in overview.chromnames:
    begin,end=overview.extent(name)
    ax.axhline(begin-.5,color='grey',lw=.2); ax.axvline(begin-.5,color='grey',lw=.2)
ax.set(title='Whole-genome raw overview',xlabel='Genomic bin',ylabel='Genomic bin')
fig.savefig(out/'qc/whole_genome_preview.png',dpi=170); plt.close(fig)
with (out/'qc/manual_review.tsv').open('w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['check','status','notes'])
    for check in ['MNase digestion trace and protocol-specific target','Final library size/adapter-dimer trace',
                  'Sample/reference identity','Raw and post-trim FastQC flags',
                  'Short-range strand bias and distance cutoff','P(s) versus matched controls',
                  'Coverage gaps, trans blocks, CNV/rearrangements','Raw and balanced map appearance',
                  'Replicate concordance at matched depth (when replicates exist)',
                  'Resolution and pooling decision']:
        w.writerow([check,'NOT_REVIEWED',''])
PY
# MultiQC is optional; the dedicated TSV/PNG reports above are always produced.
if command -v multiqc >/dev/null; then
    multiqc "$OUTDIR/qc" -o "$OUTDIR/qc/multiqc" > "$OUTDIR/logs/multiqc.log" 2>&1 || \
        echo 'WARNING: MultiQC failed; inspect individual QC reports.'
fi
printf '\nCompleted. Inspect QC BEFORE using balanced matrices.\n'
printf 'Standalone matrix: %s\n' "$OUTDIR/maps/${SAMPLE}.${RESOLUTION}.cool"
printf 'Multiple resolutions: %s\n' "$MCOOL"
printf 'Contact QC: %s\n' "$OUTDIR/qc/contact_qc.tsv"
printf 'Orientation QC: %s\n' "$OUTDIR/qc/orientation_and_scaling.png"
printf 'Resolution QC: %s\n' "$OUTDIR/qc/resolution_qc.tsv"
printf 'Manual review: %s\n' "$OUTDIR/qc/manual_review.tsv"
# Keep independent replicate matrices for reproducibility/differential analyses.
# When pooling accepted libraries, SUM raw counts and re-balance the pooled map.
# Intermediate sorted/pairs files are retained and can be large; remove only when
# you have archived the necessary pairs/QC and no longer need to revisit filters.
