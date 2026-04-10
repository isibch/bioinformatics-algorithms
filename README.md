# Bioinformatics Algorithms

A collection of classical sequence alignment algorithms implemented from scratch in Python.

The goal of this project is to understand and implement the core algorithms used in bioinformatics for pairwise sequence alignment.

Implemented algorithms include global, local and affine-gap alignments with a shared API, CLI interface, visualization tools and benchmarking utilities.

---

# Implemented Algorithms

### Needleman–Wunsch
Global sequence alignment with linear gap penalties.

### Smith–Waterman
Local sequence alignment algorithm for detecting high-scoring subsequences.

### Gotoh Algorithm
Affine gap penalty alignment using three dynamic programming matrices.

Variants implemented:

- Gotoh Global
- Gotoh Local
- Gotoh Free-shift

---

# Installation

Clone the repository
```bash
git clone https://github.com/yourusername/bioinformatics-algorithms.git
cd bioinformatics-algorithms
```

Run the CLI from the project root:

```bash
PYTHONPATH=src python3 -m bioalgos.cli GATTACA GCATGCU --method nw
```

---

# Usage

Example:

```bash
PYTHONPATH=src python3 -m bioalgos.cli GATTACA GCATGCU --method nw
```

Example Output:
```
Method: nw
Score: 0
Checked score: 0

G-ATTACA
| | |.|.
GCA-TGCU

Stats:
  matches: 4
  mismatches: 2
  gaps: 2
  alignment_length: 8
  identity: 0.5
```

Supported methods:
```
nw
sw
gotoh-global
gotoh-local
gotoh-freeshift
```

---

# Benchmark

Run the built-in benchmark:

```bash
PYTHONPATH=src python3 scripts/benchmark_alignment.py
```
Example Outout:

```
Benchmarking alignment methods

Sequences: GATTACA vs GCATGCU
method             score    checked  avg_time_ms
----------------------------------------------------
nw                 0        0        0.02
sw                 2        2        0.02
gotoh-global       -1       -1       0.05
gotoh-local        2        2        0.05
gotoh-freeshift    1        1        0.05
```

---

# Project Structure

```
bioinformatics-algorithms
│
├── src
│   └── bioalgos
│       ├── alignment
│       │   ├── __init__.py
│       │   ├── api.py
│       │   ├── needleman_wunsch.py
│       │   ├── smith_waterman.py
│       │   ├── gotoh_global.py
│       │   ├── gotoh_local.py
│       │   ├── gotoh_freeshift.py
│       │   ├── utils.py
│       │   └── visualization.py
│       │
│       └── cli.py
│
├── scripts
│   └── benchmark_alignment.py
│
├── README.md
└── .gitignore
```

---

# Future Work

Planned extensions:
	•	substitution matrices (BLOSUM / PAM)
	•	FASTA file support
	•	multiple sequence alignment
	•	k-mer indexing
	•	BLAST-like heuristic search
	•	phylogenetic tree algorithms (UPGMA / Neighbor Joining)