#!/bin/bash

# Default configuration
CONFIG_FILE="pipeline_config.sh"

# Load configuration
if [[ -f "$CONFIG_FILE" ]]; then
    source "$CONFIG_FILE"
else
    echo "Configuration file '$CONFIG_FILE' not found."
    exit 1
fi

# Get today's date in YYYY_MM_DD format
date=$(date +"%Y_%m_%d")

# Default output directory
outdir="${OUTDIR:-$date}"

# Create output directories
mkdir -p "$SRADIR/$outdir/fastq"
mkdir -p "$SRADIR/$outdir/alignments"
mkdir -p "$SRADIR/$outdir/counts"

# Function to download SRA files
download_sra_files() {
    local sra_file=$1
    echo "Downloading $sra_file..."
    prefetch -O "$SRADIR/$outdir" "$sra_file"
}

# Function to convert SRA file to FASTQ
convert_to_fastq() {
    local srafile="$1"
    local filename=$(basename "$srafile" .sra)
    
    echo "Converting $filename to FASTQ..."
    fastq-dump --outdir "$SRADIR/$outdir/fastq" --split-files "$srafile" --verbose
    echo "Converted $filename.sra to FASTQ files"
}


run_star_alignment() {
    local sra_file=$1
    local bam_file="$SRADIR/$outdir/alignments/${sra_file}Aligned.sortedByCoord.out.bam"
    local log_file="$SRADIR/$outdir/alignments/${sra_file}Log.out"

    # Check if alignment has already been successfully completed
    if [[ -f "$bam_file" && -f "$log_file" && $(tail -n 1 "$log_file") == "ALL DONE!" ]]; then
        echo "Alignment and log check complete for $sra_file, skipping STAR alignment."
        return 0  # Indicate success without rerunning
    else
        echo "Aligning reads for $sra_file..."
        STAR --genomeDir "$GENOME_INDEX_DIR" \
             --runThreadN "$THREADS" \
             --readFilesIn "$SRADIR/$outdir/fastq/${sra_file}_1.fastq" "$SRADIR/$outdir/fastq/${sra_file}_2.fastq" \
             --outFileNamePrefix "$SRADIR/$outdir/alignments/${sra_file}" \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMunmapped Within \
             --outSAMattributes Standard

        # Verify if STAR was successful and create a simple log statement
        if [ $? -eq 0 ]; then
            echo "ALL DONE!" >> "$log_file"  # Append completion text to log file
            return 0  # Success
        else
            echo "STAR alignment for $sra_file failed."
            return 1  # Fail
        fi
    fi
}


# Function to run featureCounts
run_feature_counts() {
    local sra_file=$1
    echo "Running featureCounts for $sra_file..."
    featureCounts -a "$REF_GENOME_GTF" \
                  -o "$SRADIR/$outdir/counts/${sra_file}.counts.txt" \
                  -p -T "$THREADS" \
                  "$SRADIR/$outdir/alignments/${sra_file}Aligned.sortedByCoord.out.bam"
}

# Read the list of SRA accession numbers
mapfile -t SRA_FILES < "$SRA_LIST"

# Download SRA files
echo "Starting download of SRA files..."
for sra_file in "${SRA_FILES[@]}"; do
    download_sra_files "$sra_file" &
done
wait

# Find all SRA files and convert them to FASTQ format
echo "Converting SRA files to FASTQ format..."
sra_files=("$SRADIR/$outdir"/SRR*/*.sra)
for sra_file in "${sra_files[@]}"; do
    convert_to_fastq "$sra_file"
done

# Align and count features
echo "Processing aligned reads..."
for sra_file in "${SRA_FILES[@]}"; do
    run_star_alignment "$sra_file"
    if [ $? -eq 0 ]; then
        run_feature_counts "$sra_file"
    else
        echo "STAR alignment for $sra_file failed. Skipping featureCounts."
    fi
    # Optionally delete fastq files to save space
    rm -f "$SRADIR/$outdir/fastq/${sra_file}"*.fastq
done
