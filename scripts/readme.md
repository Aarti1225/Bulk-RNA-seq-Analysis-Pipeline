# Clean, Reproducible Version of Pipeline (rna_seq_pipeline.sh)

#!/bin/bash
set -euo pipefail

# ===============================
# RNA-Seq Pipeline (WSL/Ubuntu)
# ===============================

# ---- 0. Update system & install dependencies ----
sudo apt update
sudo apt install -y sra-toolkit hisat2 samtools fastqc multiqc default-jre

# ---- 1. Download SRA files ----
# Example: prefetch SRR7179504
# Or use python3 fastq_download.py for batch download
prefetch SRR7179504

# ---- 2. Convert SRA to FASTQ ----
fastq-dump --outdir fastq --gzip --skip-technical \
  --readids --read-filter pass --dumpbase --split-3 --clip SRR7179504.sra

# ---- 3. Quality control ----
fastqc fastq/*.fastq.gz -o fastqc_results/ --threads 8
multiqc fastqc_results/ -o multiqc_report/

# ---- 4. (Optional) Trimming ----
java -jar Trimmomatic-0.39/trimmomatic-0.39.jar SE \
  -threads 4 fastq/SRR7179504.fastq.gz \
  fastq/SRR7179504_trimmed.fastq.gz TRAILING:10 -phred33

# ---- 5. Concatenate & Rename FASTQ files ----
cat SRR7179504_pass.fastq.gz SRR7179505_pass.fastq.gz \
    SRR7179506_pass.fastq.gz SRR7179507_pass.fastq.gz > LNCAP_Normoxia_S1.fastq.gz

cat SRR7179508_pass.fastq.gz SRR7179509_pass.fastq.gz \
    SRR7179510_pass.fastq.gz SRR7179511_pass.fastq.gz > LNCAP_Normoxia_S2.fastq.gz

cat SRR7179520_pass.fastq.gz SRR7179521_pass.fastq.gz \
    SRR7179522_pass.fastq.gz SRR7179523_pass.fastq.gz > LNCAP_Hypoxia_S1.fastq.gz

cat SRR7179524_pass.fastq.gz SRR7179525_pass.fastq.gz \
    SRR7179526_pass.fastq.gz SRR7179527_pass.fastq.gz > LNCAP_Hypoxia_S2.fastq.gz

mv SRR7179536_pass.fastq.gz PC3_Normoxia_S1.fastq.gz
mv SRR7179537_pass.fastq.gz PC3_Normoxia_S2.fastq.gz
mv SRR7179540_pass.fastq.gz PC3_Hypoxia_S1.fastq.gz
mv SRR7179541_pass.fastq.gz PC3_Hypoxia_S2.fastq.gz

# ---- 6. Download genome index (GRCh38) ----
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -xvzf grch38_genome.tar.gz

# ---- 7. Align reads with HISAT2 ----
hisat2 -q -x grch38/genome -U LNCAP_Hypoxia_S1.fastq.gz | \
  samtools sort -o LNCAP_Hypoxia_S1.bam
samtools index LNCAP_Hypoxia_S1.bam

# (Repeat alignment for all samples; or make a loop)

# ---- 8. Count features with featureCounts ----
mkdir -p quants
featureCounts -S 2 -a Homo_sapiens.GRCh38.114.gtf \
  -o quants/featurecounts.txt *.bam

# ---- 9. QC of aligned reads ----
./qualimap_v2.3/qualimap rnaseq \
  -bam LNCAP_Hypoxia_S1.bam \
  -gtf Homo_sapiens.GRCh38.114.gtf \
  -outdir rnaseq_qc_results \
  --java-mem-size=8G

# ----10. Differential Expression in R (DESeq2)----
move to R
  

[fastq_download.py](https://github.com/user-attachments/files/22005205/fastq_download.py)

import time
import subprocess
sra_numbers = [
    "SRR7179504", "SRR7179505", "SRR7179506", "SRR7179507",
    "SRR7179508", "SRR7179509", "SRR7179510", "SRR7179511",
    "SRR7179520", "SRR7179521", "SRR7179522", "SRR7179523",
    "SRR7179524", "SRR7179525", "SRR7179526", "SRR7179527",
    "SRR7179536", "SRR7179537", "SRR7179540","SRR7179541"
    ]

for sra_id in sra_numbers:
    print("\n=== Downloading:", sra_id, "===")
    prefetch_cmd = f"prefetch {sra_id}"
    print("Command:", prefetch_cmd)

    start_time = time.time()
    subprocess.call(prefetch_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ Download time for {sra_id}: {elapsed_min:.2f} minutes")

for sra_id in sra_numbers:
    sra_path = f"/home/smriti_baskworkspace/bulkrnaseq_tutorial/{sra_id}/{sra_id}.sra"
    print("\n=== Generating FASTQ for:", sra_id, "===")
    fastq_dump_cmd = (
        f"fastq-dump --outdir fastq --gzip --skip-technical "
        f"--readids --read-filter pass --dumpbase --split-3 --clip {sra_path}"
    )
    print("Command:", fastq_dump_cmd)

    start_time = time.time()
    subprocess.call(fastq_dump_cmd, shell=True)
    end_time = time.time()

    elapsed_min = (end_time - start_time) / 60
    print(f"⏱ FASTQ generation time for {sra_id}: {elapsed_min:.2f} minutes")

[convert_sra.sh](https://github.com/user-attachments/files/22005262/convert_sra.sh)

#!/bin/bash

# Make sure fastq folder exists
mkdir -p ~/fastq

# Loop over all .sra files
for sra in $(find ~ -name "*.sra"); do
    base=$(basename $sra .sra)   # Extract SRR ID
    out=~/fastq/${base}.fastq.gz # Expected output file

    if [ -f "$out" ]; then
        echo "✅ $out already exists, skipping..."
    else
        echo "➡️ Converting $sra to $out ..."
        fastq-dump --outdir ~/fastq \
          --gzip --skip-technical --readids \
          --read-filter pass --dumpbase \
          --split-3 --clip "$sra"
    fi
done


[run_sra.sh](https://github.com/user-attachments/files/22005286/run_sra.sh)
[#!/bin/bash

# List of SRA IDs
sra_numbers=(
    SRR7179504 SRR7179505 SRR7179506 SRR7179507
    SRR7179508 SRR7179509 SRR7179510 SRR7179511
    SRR7179520 SRR7179521 SRR7179522 SRR7179523
    SRR7179524 SRR7179525 SRR7179526 SRR7179527
    SRR7179536 SRR7179537 SRR7179540 SRR7179541
)

for SRR in "${sra_numbers[@]}"; do
    echo -e "\n=== Processing: $SRR ==="

    # Ensure directory exists
    mkdir -p "$SRR"

    # Clean up any stale lock files
    if [ -f "$SRR/$SRR.sra.lock" ]; then
        echo "⚠️ Found lock file for $SRR — removing it..."
        rm -f "$SRR/$SRR.sra.lock"
    fi

    # Download SRA file if not already present
    if [ ! -f "$SRR/$SRR.sra" ]; then
        echo "Downloading $SRR ..."
        prefetch "$SRR" --output-directory "$SRR"
    else
        echo "$SRR already downloaded ✅"
    fi

    # Check if FASTQ already exists
    if ls fastq/${SRR}*.fastq.gz 1> /dev/null 2>&1; then
        echo "FASTQ for $SRR already exists ✅ — skipping conversion."
    else
        # Convert to FASTQ
        echo "Converting $SRR to FASTQ ..."
        fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass \
                   --dumpbase --split-3 --clip "$SRR/$SRR.sra"
    fi
done

Uploading run_sra.sh…]()


[run_fastqc.sh](https://github.com/user-attachments/files/22005231/run_fastqc.sh)


#!/bin/bash

# Input folder with fastq files
FASTQ_DIR=~/fastq

# Output folder for FastQC results
QC_DIR=~/qc_fastqc
mkdir -p $QC_DIR

# Loop through all fastq.gz files
for fq in $FASTQ_DIR/*.fastq.gz; do
    base=$(basename $fq .fastq.gz)
    out_zip=$QC_DIR/${base}_fastqc.zip
    out_html=$QC_DIR/${base}_fastqc.html

    if [[ -f "$out_zip" && -f "$out_html" ]]; then
        echo " FastQC already done for $base, skipping..."
    else
        echo " Running FastQC on $fq ..."
        fastqc -o $QC_DIR -f fastq $fq
    fi
done

echo " FastQC finished. Results in $QC_DIR"




[make_grch38.sh](https://github.com/user-attachments/files/22005300/make_grch38.sh)
#!/bin/sh

#
# Downloads sequence for the GRCh38 release 84 version of H. sapiens (human) from
# Ensembl.
#
# Note that Ensembl's GRCh38 build has three categories of compressed fasta
# files:
#
# The base files, named ??.fa.gz
#
# By default, this script builds and index for just the base files,
# since alignments to those sequences are the most useful.  To change
# which categories are built by this script, edit the CHRS_TO_INDEX
# variable below.
#

ENSEMBL_RELEASE=84
ENSEMBL_GRCh38_BASE=ftp://ftp.ensembl.org/pub/release-${ENSEMBL_RELEASE}/fasta/homo_sapiens/dna

get() {
	file=$1
	if ! wget --version >/dev/null 2>/dev/null ; then
		if ! curl --version >/dev/null 2>/dev/null ; then
			echo "Please install wget or curl somewhere in your PATH"
			exit 1
		fi
		curl -o `basename $1` $1
		return $?
	else
		wget $1
		return $?
	fi
}

HISAT2_BUILD_EXE=./hisat2-build
if [ ! -x "$HISAT2_BUILD_EXE" ] ; then
	if ! which hisat2-build ; then
		echo "Could not find hisat2-build in current directory or in PATH"
		exit 1
	else
		HISAT2_BUILD_EXE=`which hisat2-build`
	fi
fi

rm -f genome.fa
F=Homo_sapiens.GRCh38.dna.primary_assembly.fa
if [ ! -f $F ] ; then
	get ${ENSEMBL_GRCh38_BASE}/$F.gz || (echo "Error getting $F" && exit 1)
	gunzip $F.gz || (echo "Error unzipping $F" && exit 1)
	mv $F genome.fa
fi

CMD="${HISAT2_BUILD_EXE} -p 4 genome.fa genome"
echo Running $CMD
if $CMD ; then
	echo "genome index built; you may remove fasta files"
else
	echo "Index building failed; see error message"
fi



[hisat2alignment.sh](https://github.com/user-attachments/files/22005195/hisat2alignment.sh)
#!/bin/bash

# Directory containing FASTQ files
FASTQ_DIR="fastq"
GENOME_INDEX="grch38/genome"
LOGFILE="alignment_log.txt"

# Clear or create logfile
> $LOGFILE

# List of FASTQ files
FILES=(
    "LNCAP_Hypoxia_S2.fastq.gz"
    "LNCAP_Normoxia_S1.fastq.gz"
    "LNCAP_Normoxia_S2.fastq.gz"
    "PC3_Hypoxia_S1.fastq.gz"
    "PC3_Hypoxia_S2.fastq.gz"
    "PC3_Normoxia_S1.fastq.gz"
    "PC3_Normoxia_S2.fastq.gz"
)

# Loop through each file
for f in "${FILES[@]}"; do
    SAMPLE_NAME=$(basename "$f" .fastq.gz)
    echo "Processing $SAMPLE_NAME at $(date)" | tee -a $LOGFILE
    START_TIME=$(date +%s)

    # Run HISAT2 alignment and Samtools sorting/indexing
    hisat2 -q -x $GENOME_INDEX -U $FASTQ_DIR/$f | \
    samtools sort -o ${SAMPLE_NAME}.bam
    samtools index ${SAMPLE_NAME}.bam

    END_TIME=$(date +%s)
    ELAPSED=$((END_TIME - START_TIME))
    echo "Finished $SAMPLE_NAME in $ELAPSED seconds at $(date)" | tee -a $LOGFILE
    echo "--------------------------------------" | tee -a $LOGFILE

done

echo "All files processed successfully at $(date)" | tee -a $LOGFILE



[generating_countsmatrix(2).ipynb](https://github.com/user-attachments/files/22005183/generating_countsmatrix.2.ipynb){
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 1,
   "id": "f8f042ec",
   "metadata": {},
   "outputs": [
    {
     "name": "stderr",
     "output_type": "stream",
     "text": [
      "C:\\Users\\smriti\\Downloads\\jupyter\\Lib\\site-packages\\pandas\\core\\arrays\\masked.py:60: UserWarning: Pandas requires version '1.3.6' or newer of 'bottleneck' (version '1.3.5' currently installed).\n",
      "  from pandas.core import (\n"
     ]
    },
    {
     "name": "stdout",
     "output_type": "stream",
     "text": [
      "Files found: ['C:/Users/smriti/Desktop/quants\\\\LNCAP_Hypoxia_S1_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\LNCAP_Hypoxia_S2_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\LNCAP_Normoxia_S1_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\LNCAP_Normoxia_S2_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\PC3_Hypoxia_S1_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\PC3_Hypoxia_S2_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\PC3_Normoxia_S1_featurecounts.txt', 'C:/Users/smriti/Desktop/quants\\\\PC3_Normoxia_S2_featurecounts.txt']\n",
      "Completed LNCAP_Hypoxia_S1 | Rows: 78894 | Time: 0.02 min\n",
      "Completed LNCAP_Hypoxia_S2 | Rows: 78894 | Time: 0.02 min\n",
      "Completed LNCAP_Normoxia_S1 | Rows: 78894 | Time: 0.02 min\n",
      "Completed LNCAP_Normoxia_S2 | Rows: 78894 | Time: 0.02 min\n",
      "Completed PC3_Hypoxia_S1 | Rows: 78894 | Time: 0.02 min\n",
      "Completed PC3_Hypoxia_S2 | Rows: 78894 | Time: 0.01 min\n",
      "Completed PC3_Normoxia_S1 | Rows: 78894 | Time: 0.02 min\n",
      "Completed PC3_Normoxia_S2 | Rows: 78894 | Time: 0.02 min\n",
      "\\All files processed!\n",
      "Merged matrix shape: (78894, 9)\n",
      "Saved to: C:/Users/smriti/Desktop/quants\\GSE106305_counts_matrix.csv\n"
     ]
    }
   ],
   "source": [
    "import os\n",
    "import glob\n",
    "import pandas as pd\n",
    "import time\n",
    "\n",
    "path = \"C:/Users/smriti/Desktop/quants\"\n",
    "files = glob.glob(os.path.join(path, \"*.txt\"))\n",
    "\n",
    "print(\"Files found:\", files)\n",
    "\n",
    "all_counts = []\n",
    "\n",
    "for file in files:\n",
    "    start_time = time.time()\n",
    "    df = pd.read_csv(file, sep=\"\\t\", comment=\"#\")\n",
    "    \n",
    "    sample_name = os.path.basename(file).replace(\"_featurecounts.txt\", \"\")\n",
    "    df = df[[\"Geneid\", df.columns[-1]]]\n",
    "    df.rename(columns={df.columns[-1]: sample_name}, inplace=True)\n",
    "    \n",
    "    all_counts.append(df)\n",
    "    \n",
    "    elapsed = (time.time() - start_time) / 60  # minutes\n",
    "    print(f\"Completed {sample_name} | Rows: {df.shape[0]} | Time: {elapsed:.2f} min\")\n",
    "\n",
    "counts_matrix = all_counts[0]\n",
    "for df in all_counts[1:]:\n",
    "    counts_matrix = counts_matrix.merge(df, on=\"Geneid\", how=\"outer\")\n",
    "\n",
    "output_file = os.path.join(path, \"GSE106305_counts_matrix.csv\")\n",
    "counts_matrix.to_csv(output_file, index=False)\n",
    "\n",
    "print(\"\\All files processed!\")\n",
    "print(\"Merged matrix shape:\", counts_matrix.shape)\n",
    "print(\"Saved to:\", output_file)"
   ]
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.11.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}

[featurecounts.sh](https://github.com/user-attachments/files/22005191/featurecounts.sh)
#!/bin/bash

# Go to your aligned reads folder
cd /home/smriti_baskworkspace/bulkrnaseq_tutorial/alignedreads

# Loop over all BAM files
for bam in *.bam; do
    start=$(date +%s)  # start time

    echo "Processing $bam ..."
    featureCounts -s 0 -a /home/smriti_baskworkspace/bulkrnaseq_tutorial/Homo_sapiens.GRCh38.114.gtf \
        -o /home/smriti_baskworkspace/bulkrnaseq_tutorial/quants/${bam%.bam}_featurecounts.txt \
        "$bam"

    end=$(date +%s)  # end time
    runtime=$(( (end - start) / 60 ))  # in minutes

    echo "✅ Completed $bam in $runtime minutes."
    echo "------------------------------------"
done


