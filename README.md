<img width="1192" height="261" alt="rnaseq_pipeline_diagram" src="https://github.com/user-attachments/assets/600fe0a6-2f59-4780-a494-3ad9b5d72d6b" />

# Bulk-RNA-seq-Analysis-Pipeline 🧬
Welcome to the **Bulk RNA-seq Pipeline** – a clean, reproducible, and beginner-friendly workflow for processing RNA-seq data, starting from raw **SRA IDs** and ending with a **gene counts matrix** ready for downstream **differential expression analysis**.
This repository contains a step-by-step pipeline for analyzing bulk RNA-seq data from the study by Guo et al., Nature Communications (2019) — GSE106305.
We focus on LNCaP and PC3 prostate cancer cell lines, comparing normoxia vs hypoxia conditions to identify differentially expressed genes (DEGs).
This pipeline is designed to help **bioinformaticians, biologists, and beginners** run RNA-seq analysis step by step, without being overwhelmed by scattered commands. Everything is modular, documented, and ready to extend.

---

## 🔍 Overview

The pipeline covers:

1. **Data Retrieval** → Download `.sra` files from NCBI SRA and convert to `.fastq.gz`
2. **Quality Control** → Run **FastQC** and summarize with **MultiQC**
3. **(Optional) Trimming** → Clean low-quality reads with **Trimmomatic**
4. **Alignment** → Map reads to the **human genome (GRCh38)** using **HISAT2**
5. **Post-processing** → Sort & index alignments with **Samtools**
6. **Quantification** → Count reads per gene using **featureCounts**
7. **QC of Aligned Reads** → Assess alignments with **Qualimap**
8. **Output** → A clean counts matrix for analysis in R (DESeq2/edgeR)

---

## ⚙️ Why these tools?

* **HISAT2** → memory-efficient and splice-aware aligner (great for human RNA-seq)
* **FastQC + MultiQC** → easy visualization of raw read quality
* **Trimmomatic (optional)** → only used if adapter contamination or poor bases are detected
* **featureCounts** → fast, multi-threaded gene quantification
* **Qualimap** → intuitive BAM-level QC (coverage, bias, duplication)

💡 **Tip:** Trimming is often unnecessary if FastQC reports look good. Avoid trimming unless needed to prevent data loss.

---

## 📂 Repository Structure

```
bulk-rnaseq-pipeline/
├── README.md                  <- This file
├── LICENSE                    <- License (MIT recommended)
├── requirements.txt           <- Python dependencies
├── environment.yml            <- Conda environment (optional)
├── config/
│   ├── sra_ids.txt            <- SRA IDs (one per line)
│   └── samples.csv            <- Metadata: sample, condition, fastq path
├── scripts/
│   ├── download_fastq.py       <- Download & convert SRA to FASTQ
│   ├── qc.sh                   <- FastQC + MultiQC
│   ├── trimming.sh             <- Trimmomatic (optional)
│   ├── alignment.sh            <- HISAT2 + Samtools
│   ├── counts_featurecounts.sh <- Quantification
│   └── qualimap_qc.sh          <- BAM QC
├── results/                   <- Outputs (not uploaded, large)
│   ├── fastqc/
│   ├── multiqc/
│   ├── aligned_bam/
│   └── counts/
└── R/
    ├── dge_analysis.R          <- DESeq2 analysis
    └── plots.R                 <- PCA, volcano, heatmaps
```

---

## 🔧 Installation

### Install required tools (Ubuntu / WSL)

```bash
sudo apt update
sudo apt install -y sra-toolkit fastqc multiqc hisat2 samtools
```

### Install Python dependencies

```bash
pip install -r requirements.txt
```

### (Optional) Conda environment

```bash
conda env create -f environment.yml
conda activate bulk-rnaseq
```

---

## 🚀 Usage

### 1️⃣ Download FASTQ files from SRA

Edit `config/sra_ids.txt` and run:

```bash
python scripts/download_fastq.py
```

### 2️⃣ Run Quality Control

```bash
bash scripts/qc.sh
```

### 3️⃣ (Optional) Trim Reads

```bash
bash scripts/trimming.sh
```

### 4️⃣ Align Reads

```bash
bash scripts/alignment.sh
```

### 5️⃣ Count Genes

```bash
bash scripts/counts_featurecounts.sh
```

### 6️⃣ QC of Alignments

```bash
bash scripts/qualimap_qc.sh
```

### 7️⃣ Differential Expression in R

```bash
Rscript R/dge_analysis.R
```

---

## 📊 Outputs

* **QC reports** → FastQC HTMLs + MultiQC summary
* **Aligned BAMs** → Sorted & indexed BAM files
* **Counts matrix** → Gene expression counts (ready for DE analysis)
* **Plots** → PCA, volcano, heatmap (from R)

---

## ✅ Best Practices

* Never upload raw FASTQ or BAM files to GitHub (too large) → only keep scripts + summary files.
* Use `.gitignore` to exclude `results/fastq/` and `results/aligned_bam/`.
* Record exact reference genome & GTF versions for reproducibility.
* Test pipeline on a **single sample** before scaling to all.

---

## 👩‍💻 Author

Created with ❤️ for bioinformatics enthusiasts and researchers starting with RNA-seq analysis.

📬 Contributions, suggestions, and improvements are welcome!

---

## 📜 License

This project is licensed under the **MIT License** – feel free to use, modify, and share.
