# 🧬 Transcriptomic Profiling & End-to-End RNA-Seq Pipeline: Asthma vs. Control
[![Bioinformatics](https://img.shields.io/badge/Domain-Transcriptomics%20%26%20Bioinformatics-8A2BE2?style=for-the-badge&logo=dna)](https://github.com)
[![Pipeline](https://img.shields.io/badge/Pipeline-Bowtie2%20%7C%20StringTie%20%7C%20DESeq2-2E8B57?style=for-the-badge&logo=linux)](https://github.com)
[![R-Bioconductor](https://img.shields.io/badge/R-Bioconductor_v4.3-276DC3?style=for-the-badge&logo=r)](https://bioconductor.org)
[![Status](https://img.shields.io/badge/Status-Completed-success?style=for-the-badge)](https://github.com)

---

## 📌 Executive Summary

This repository houses an end-to-end, fully reproducible transcriptomic pipeline for analyzing high-throughput RNA-Seq data. Utilizing human **Chromosome 21 (GRCh38)** as a targeted model dataset, the project investigates differential gene expression dynamics and dysregulated molecular pathways comparing **Asthmatic** airway tissue against healthy **Control** groups.

### Key Analytical Achievements:
* **High-Precision Alignment    : ** Achieved an **overall alignment rate of 97.98%** using `Bowtie2`.
* **Advanced Feature Debugging  : ** Resolved a **95.38% NA mapping anomaly** in `AnnotationDbi` via custom regular expressions, elevating annotation efficiency to **93.08%**.
* **Biostatistical Variance Control:** Captured **84% of total data variance** (PC1: 50%, PC2: 34%) using Variance Stabilizing Transformation (`VST`) in `DESeq2`.
* **Biological Biomarker Discovery:** Identified **7 core Differentially Expressed Genes (DEGs)** pointing toward a dual molecular pathology: **anti-inflammatory immune suppression** (`IL10RB`, `MX1`) coupled with **glycolytic metabolic reprogramming** (`PFKL`, `SUMO3`).

---

## 📂 Repository Structure

```text
├── data/
│   ├── metadata.csv                      # Sample metadata (Asthma vs Control)
│   └── gene_count_matrix.csv             # Raw integer count matrix from prepDE.py
├── scripts/
│   ├── 01_upstream_processing.sh         # Bash script for FastQC, Trimmomatic, Bowtie2
│   └── 02_deseq2_analysis.R              # R script for DESeq2, VST, & Regex Annotation
├── results/
│   ├── Hasil_DESeq2_Annotated_Fix.csv    # Final annotated differential expression table
│   ├── figures/
│   │   ├── Rplot_plotPCA.png             # PCA plot (Group clustering)
│   │   ├── Rplot_VolcanoPlot.pdf         # Differential expression volcano plot
│   │   ├── Rplot_pheatmap.pdf            # Top 20 DEGs clustered heatmap
│   │   └── string_hires_image.png        # STRING DB Protein-Protein Interaction Network
│   └── laporan-metodologi-rna-seq.pdf    # Formal methodology & results research report
└── README.md                             # Project documentation
```

---

## ⚡ Pipeline Workflow Architecture

```text
  [ Raw FASTQ Reads ] 
           │
           ▼ (FastQC v0.11.9 & Trimmomatic v0.39)
  [ Cleaned Reads ] ───► 64.19% Retained (Phred > 30)
           │
           ▼ (Bowtie2 v2.4.5 - GRCh38 Chr21)
  [ SAM / Sorted BAM ] ───► 97.98% Overall Alignment Rate
           │
           ▼ (StringTie v2.2.1 & prepDE.py)
  [ Raw Count Matrix ] ───► gene_count_matrix.csv
           │
           ▼ (DESeq2 - Median-of-Ratios Normalization)
  [ Statistical Testing ] ───► Threshold: Padj ≤ 0.05 & |Log2FC| ≥ 1.0
           │
           ▼ (Regex Feature Cleaning)
  [ Annotated DEGs ] ───► 93.08% Successful Annotation Success
           │
           ▼ (STRING DB v12.0)
  [ PPI Network ] ───► Immune Suppression & Metabolic Reprogramming Axis
```

---

## 📊 Key Biological Findings (7 Core DEGs)

Filtering parameters: **Adjusted P-value (\\(P_{adj}\\)) \\(\le 0.05\\)** and **\\(|\log_2\text{Fold Change}| \ge 1.0\\)**.

| Official Symbol | Log2 Fold Change | Adjusted \\(P\\)-value (\\(P_{adj}\\)) | Expression Status | Primary Biological Function / Pathology |
| :--- | :---: | :---: | :---: | :--- |
| **`PFKL`** | **+8.73** | \\(7.76 \times 10^{-6}\\) | **Upregulated** | **Glikolisis (Rate-Limiting):** Mendorong *metabolic reprogramming* untuk pasokan energi sel peradangan. |
| **`SUMO3`** | **+24.07** | \\(8.42 \times 10^{-6}\\) | **Upregulated** | **Modifikasi Pasca-Translasi:** Sumoilasi protein untuk menstabilkan faktor transkripsi saat stres seluler. |
| **`IL10RB`** | **-11.85** | \\(1.69 \times 10^{-11}\\) | **Downregulated** | **Reseptor Anti-inflamasi:** Subunit reseptor IL-10/IL-22; hilangnya penekan peradangan saluran napas. |
| **`KCNJ15`** | **-25.39** | \\(4.83 \times 10^{-6}\\) | **Downregulated** | **Saluran Ion Kalium:** Mengatur homeostasis potensial membran sel epitel dan imunitas. |
| **`MX1`** | **-8.70** | \\(0.019\\) | **Downregulated** | **GTPase Antiviral:** Komponen utama imunitas bawaan (*innate immunity*) dan respons interferon. |
| **`MCM3AP`** | **-8.76** | \\(6.10 \times 10^{-5}\\) | **Downregulated** | **Regulasi Siklus Sel:** Berperan dalam transpor mRNA dan inisiasi replikasi DNA. |
| **`RPS26P4`** | **-2.08** | \\(1.34 \times 10^{-7}\\) | **Downregulated** | **Pseudogen Ribosom:** Komponen struktural RNA. |

---

## 🛠️ Data Engineering & Troubleshooting Log

<details>
<summary><b>🔍 1. Resolution of 95.38% Missing Values in Gene Annotation (Click to Expand)</b></summary>

<br>

* **Problem:** Direct execution of `AnnotationDbi::mapIds()` using raw count matrix row names yielded **95.38% NA values** (124/130 genes).
* **Root Cause Analysis:** StringTie output concats gene symbols with pipes (e.g., `ENSG00000141959|PFKL`) and Ensembl ID version suffixes (e.g., `.10`), which fail exact matching against `org.Hs.eg.db`.
* **Technical Solution:** Implemented a two-tiered cleaning approach in R:
  ```R
  # Remove pipe delimiters and version numbers simultaneously
  ensembl_clean <- sub("[|\\.].*", "", rownames(res))
  
  # Primary mapping via org.Hs.eg.db
  res$symbol <- mapIds(org.Hs.eg.db, keys=ensembl_clean, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  
  # Secondary fallback to StringTie extracted symbols
  stringtie_symbol <- sub(".*\\|", "", rownames(res))
  res$symbol <- ifelse(is.na(res$symbol), stringtie_symbol, res$symbol)
  ```
* **Impact:** Reduced annotation missingness from **95.38% down to 6.92%**, restoring high-confidence biological identity to the dataset.

</details>

<details>
<summary><b>🔍 2. DESeq2 Variance Stabilizing Transformation (VST) Sub-sampling Fix</b></summary>

<br>

* **Problem:** Standard `vst(dds)` failed with `Error: less than 'nsub' rows`.
* **Root Cause:** Default `vst()` expects \\(n \ge 1000\\) genes for random sub-sampling, whereas Chromosome 21 subsetting produced \\(< 1000\\) rows.
* **Technical Solution:** Invoked `varianceStabilizingTransformation(dds, blind = FALSE)` directly.
* **Impact:** Generated robust variance-stabilized values enabling accurate Principal Component Analysis (PCA).

</details>

---

## 🚀 Reproduction & Execution Guide

### Prerequisites
* Linux Ubuntu 22.04 LTS via WSL2
* R (v4.3+) & RStudio
* Installed CLI Tools: `FastQC`, `Trimmomatic`, `Bowtie2`, `Samtools`, `StringTie`

### Step 1: Upstream Processing (Linux CLI)
```bash
# Clone Repository
git clone https://github.com/your-username/RNASeq-Asthma-DGE-Analysis.git
cd RNASeq-Asthma-DGE-Analysis/scripts

# Execute Quality Control & Trimming
trimmomatic PE -phred33 input_R1.fq input_R2.fq \
  paired_R1.fq unpaired_R1.fq paired_R2.fq unpaired_R2.fq \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

# Alignment & BAM Management
bowtie2 -x chr21_index -1 paired_R1.fq -2 paired_R2.fq -S sample.sam
samtools view -bS sample.sam | samtools sort -o sample.sorted.bam
samtools index sample.sorted.bam

# Transcript Assembly & Count Matrix Generation
stringtie -p 4 -G chr21.gtf -o sample.gtf -e sample.sorted.bam
python prepDE.py -i sample_list.txt
```

### Step 2: Downstream Analysis (RStudio)
```R
source("scripts/02_deseq2_analysis.R")
# Outputs annotated CSV and publication-ready plots to results/
```

---

## 🕸️ Network Pathology (STRING DB Analysis)

Reconstruction of the Protein-Protein Interaction (PPI) network on **STRING DB v12.0** revealed a functional core connected component: **`IL10RB` — `MX1` — `PFKL` — `SUMO3`**:

```text
   [IL10RB] ───────── [MX1] ───────── [PFKL] ───────── [SUMO3]
 (Immune Brake)   (Antiviral)     (Glycolysis)    (Sumoylation)
```

* **Immune Dysregulation Axis:** Downregulation of `IL10RB` (anti-inflammatory signaling) and `MX1` (interferon innate response) compromises the tissue's capacity to resolve inflammation.
* **Metabolic Reprogramming Axis:** Marked upregulation of `PFKL` indicates an inflammatory metabolic switch toward hyper-glycolysis, supported by `SUMO3`-mediated post-translational stabilization.

---

## 📄 Documentation & References

* **Full Methodology Report:** See [`results/laporan-metodologi-rna-seq.pdf`](results/laporan-metodologi-rna-seq.pdf) for the publication-ready academic report.
* **Reference Genome:** Ensembl Human GRCh38 (Chromosome 21).

---

## 👤 Author & Contact

**Gracecielo Isaiah**  
*Aspirant Biological Data Analyst & Bioinformatics Consultant*  
* **Focus Areas:** Computational Transcriptomics, Biomarker Discovery, & Functional Omics  
* **Tools Mastered:** R/Bioconductor, Python, Linux CLI, DESeq2, STRING DB  
```

