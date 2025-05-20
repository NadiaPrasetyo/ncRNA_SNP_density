# ----------- CONFIG -----------

REFERENCE="GRCh38.primary_assembly.genome.fa"
CSV="data/SNP-densities-and-RNA.csv"
FASTQ_ROOT="data/1000genomes"
OUTPUT_DIR="data/1000genomes/output"
mkdir -p "$OUTPUT_DIR"

# Final combined outputs
COMBINED_BED="$OUTPUT_DIR/all_filtered_reads.bed"
COMBINED_BAM="$OUTPUT_DIR/combined.filtered.bam"
TMP_BAM_DIR="$OUTPUT_DIR/tmp_bams"
mkdir -p "$TMP_BAM_DIR"

> "$COMBINED_BED"  # Clear if exists

# List of FASTQ files to skip
SKIP_FILES=( 
    # ----------- CONFIG -----------

REFERENCE="GRCh38.primary_assembly.genome.fa"
CSV="data/SNP-densities-and-RNA.csv"
FASTQ_ROOT="data/1000genomes"
OUTPUT_DIR="data/1000genomes/output"
mkdir -p "$OUTPUT_DIR"

# Final combined outputs
COMBINED_BED="$OUTPUT_DIR/all_filtered_reads.bed"
COMBINED_BAM="$OUTPUT_DIR/combined.filtered.bam"
TMP_BAM_DIR="$OUTPUT_DIR/tmp_bams"
mkdir -p "$TMP_BAM_DIR"

> "$COMBINED_BED"  # Clear if exists

# List of FASTQ files to skip
SKIP_FILES=( 
    "ERR000047_2.fastq" "ERR000050_2.fastq" "ERR000057_2.fastq" "ERR000058_1.fastq" "ERR000060_1.fastq" "ERR000060_2.fastq" "ERR000064_2.fastq" "ERR000066_2.fastq" "ERR000068_1.fastq" "ERR000069_2.fastq" "ERR000071_1.fastq" "ERR000072_1.fastq" "ERR000074_1.fastq" "ERR000077_2.fastq" "ERR000080_1.fastq" "ERR000080_2.fastq" "ERR000085_1.fastq" "ERR000091_1.fastq" "ERR000094_1.fastq" "ERR000096_1.fastq" "ERR000098_2.fastq" "ERR000102_2.fastq" "ERR000132_1.fastq" "ERR000162_2.fastq" "ERR000166_2.fastq" "ERR000169_1.fastq" "ERR000169_2.fastq" "ERR000216_2.fastq" "ERR000225_1.fastq" "ERR000225_2.fastq" "ERR000228_1.fastq" "ERR000252_1.fastq" "ERR000286_2.fastq" "ERR000288_1.fastq" "ERR000292_1.fastq" "ERR000294_1.fastq" "ERR000300_1.fastq" "ERR000300_2.fastq" "ERR000307_1.fastq" "ERR000316_1.fastq" "ERR000316_2.fastq" "ERR000323_1.fastq" "ERR000333_1.fastq" "ERR000333_2.fastq" "ERR000336_1.fastq" "ERR000341_2.fastq" "ERR000345_1.fastq" "ERR000351_1.fastq" "ERR000368_1.fastq" "ERR000377_2.fastq" "ERR000381_1.fastq" "ERR000381_2.fastq" "ERR000382_1.fastq" "ERR000383_1.fastq" "ERR000389_1.fastq" "ERR000389_2.fastq" "ERR000390_1.fastq" "ERR000391_2.fastq" "ERR000396_1.fastq" "ERR000399_2.fastq" "ERR000483_1.fastq" "ERR000484_2.fastq" "ERR000485_2.fastq" "ERR000486_2.fastq" "ERR000487_2.fastq" "ERR000489_1.fastq" "ERR000489_2.fastq" "ERR000491_1.fastq" "ERR000491_2.fastq" "ERR000495_1.fastq" "ERR000504_1.fastq" "ERR000505_2.fastq" "ERR000508_2.fastq" "ERR000512_2.fastq" "ERR000514_1.fastq" "ERR000514_2.fastq" "ERR000516_1.fastq" "ERR000520_2.fastq" "ERR000523_2.fastq" "ERR000524_2.fastq" "ERR000525_1.fastq" "ERR000528_1.fastq" "ERR000529_2.fastq" "ERR000530_2.fastq" "ERR000531_1.fastq" "ERR000534_1.fastq" "ERR000534_2.fastq" "ERR000537_1.fastq" "ERR000540_1.fastq" "ERR000541_1.fastq" "ERR000545_1.fastq" "ERR000546_1.fastq" "ERR000548_1.fastq" "ERR000549_2.fastq" 
)

# List of FASTQ files to specifically process (empty = process all)
TARGET_FILES=(    
)
# ----------- PREP -----------

BED="$OUTPUT_DIR/regions.bed"
tail -n +2 "$CSV" | awk -F',' 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$BED"

# ----------- MAIN LOOP -----------

FILTERED_BAMS=()

find "$FASTQ_ROOT" -type f -name "*.fastq"| while read -r FASTQ; do
    BASENAME=$(basename "$FASTQ")
    PREFIX="${BASENAME%%.*}"

    # Skip files listed in SKIP_FILES
    if [[ " ${SKIP_FILES[@]} " =~ " ${BASENAME} " ]]; then
        echo "⏭️ Skipping $BASENAME (in skip list)"
        continue
    fi

    # If TARGET_FILES is not empty, only process those in the list
    if [[ "${#TARGET_FILES[@]}" -gt 0 ]]; then
        if [[ ! " ${TARGET_FILES[@]} " =~ " ${BASENAME} " ]]; then
            echo "⏭️ Skipping $BASENAME (not in target list)"
            continue
        fi
    fi

    echo "Processing $FASTQ → $PREFIX"

    SORTED_BAM="$TMP_BAM_DIR/$PREFIX.sorted.bam"
    FILTERED_BAM="$TMP_BAM_DIR/$PREFIX.filtered.bam"

    # 1. Align, convert, sort BAM
    bwa mem -t 2 "$REFERENCE" "$FASTQ" | \
    samtools view -Sb - | \
    samtools sort -o "$SORTED_BAM"

    if [ ! -s "$SORTED_BAM" ]; then
        echo "❌ Alignment failed for $PREFIX, skipping..."
        echo "---"
        continue
    fi

    samtools index "$SORTED_BAM"

    # 2. Filter to regions
    bedtools intersect -abam "$SORTED_BAM" -b "$BED" > "$FILTERED_BAM"
    samtools index "$FILTERED_BAM"

    READ_COUNT=$(samtools view "$FILTERED_BAM" | wc -l)

    if [[ "$READ_COUNT" -gt 0 ]]; then
        echo "✅ $PREFIX has $READ_COUNT reads mapped to regions."
        FILTERED_BAMS+=("$FILTERED_BAM")
        bedtools bamtobed -i "$FILTERED_BAM" | \
        awk -v sample="$PREFIX" 'BEGIN{OFS="\t"} {print $0, sample}' >> "$COMBINED_BED"
    else
        echo "❌ $PREFIX has NO reads mapped to regions."
    fi

    echo "---"

done

# ----------- MERGE FILTERED BAMs -----------

if [ "${#FILTERED_BAMS[@]}" -gt 0 ]; then
    echo "Merging ${#FILTERED_BAMS[@]} filtered BAMs..."
    samtools merge -f "$COMBINED_BAM" "${FILTERED_BAMS[@]}"
    samtools index "$COMBINED_BAM"
    echo "✅ Combined BAM created at: $COMBINED_BAM"
else
    echo "❌ No BAMs to merge. Exiting."
fi
 
)

# List of FASTQ files to specifically process (empty = process all)
TARGET_FILES=(    
)
# ----------- PREP -----------

BED="$OUTPUT_DIR/regions.bed"
tail -n +2 "$CSV" | awk -F',' 'BEGIN{OFS="\t"} {print $1, $2, $3, $5}' > "$BED"

# ----------- MAIN LOOP -----------

FILTERED_BAMS=()

find "$FASTQ_ROOT" -type f -name "*.fastq"| while read -r FASTQ; do
    BASENAME=$(basename "$FASTQ")
    PREFIX="${BASENAME%%.*}"

    # Skip files listed in SKIP_FILES
    if [[ " ${SKIP_FILES[@]} " =~ " ${BASENAME} " ]]; then
        echo "⏭️ Skipping $BASENAME (in skip list)"
        continue
    fi

    # If TARGET_FILES is not empty, only process those in the list
    if [[ "${#TARGET_FILES[@]}" -gt 0 ]]; then
        if [[ ! " ${TARGET_FILES[@]} " =~ " ${BASENAME} " ]]; then
            echo "⏭️ Skipping $BASENAME (not in target list)"
            continue
        fi
    fi

    echo "Processing $FASTQ → $PREFIX"

    SORTED_BAM="$TMP_BAM_DIR/$PREFIX.sorted.bam"
    FILTERED_BAM="$TMP_BAM_DIR/$PREFIX.filtered.bam"

    # 1. Align, convert, sort BAM
    bwa mem -t 2 "$REFERENCE" "$FASTQ" | \
    samtools view -Sb - | \
    samtools sort -o "$SORTED_BAM"

    if [ ! -s "$SORTED_BAM" ]; then
        echo "❌ Alignment failed for $PREFIX, skipping..."
        echo "---"
        continue
    fi

    samtools index "$SORTED_BAM"

    # 2. Filter to regions
    bedtools intersect -abam "$SORTED_BAM" -b "$BED" > "$FILTERED_BAM"
    samtools index "$FILTERED_BAM"

    READ_COUNT=$(samtools view "$FILTERED_BAM" | wc -l)

    if [[ "$READ_COUNT" -gt 0 ]]; then
        echo "✅ $PREFIX has $READ_COUNT reads mapped to regions."
        FILTERED_BAMS+=("$FILTERED_BAM")
        bedtools bamtobed -i "$FILTERED_BAM" | \
        awk -v sample="$PREFIX" 'BEGIN{OFS="\t"} {print $0, sample}' >> "$COMBINED_BED"
    else
        echo "❌ $PREFIX has NO reads mapped to regions."
    fi

    echo "---"

done

# ----------- MERGE FILTERED BAMs -----------

if [ "${#FILTERED_BAMS[@]}" -gt 0 ]; then
    echo "Merging ${#FILTERED_BAMS[@]} filtered BAMs..."
    samtools merge -f "$COMBINED_BAM" "${FILTERED_BAMS[@]}"
    samtools index "$COMBINED_BAM"
    echo "✅ Combined BAM created at: $COMBINED_BAM"
else
    echo "❌ No BAMs to merge. Exiting."
fi
