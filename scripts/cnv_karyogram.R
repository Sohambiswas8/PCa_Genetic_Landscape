# ============================================================
#  CNV Karyogram Script
# ============================================================

# Install required packages (run once)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!require("karyoploteR", quietly = TRUE))
  BiocManager::install("karyoploteR")

# Load libraries
library(karyoploteR)
library(GenomicRanges)

# Read your annotated data
cat("Reading data...\n")
# MODIFIED: Reading the new file with variant type annotations
data <- read.table("annotated_cnvs.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
cat("Found", nrow(data), "variants\n")

# Create GRanges object with variant type metadata
# MODIFIED: Added 'variant_type' to the GRanges object
cnv_ranges <- GRanges(
  seqnames = paste0("chr", data$Chrom),
  ranges = IRanges(start = data$Pos, width = 100000),  # 100kb width for visibility
  cnv_id = data$ID,
  variant_type = data$type
)

# NEW: Define a color map for the variant types
variant_types <- unique(mcols(cnv_ranges)$variant_type)
color_map <- setNames(c("#377EB8", "#E41A1C", "#4DAF4A", "#984EA3"),
                      c("copy number gain", "copy number loss", "copy number variation", "alu deletion"))
# Assign colors to each variant in the GRanges object
mcols(cnv_ranges)$color <- color_map[mcols(cnv_ranges)$variant_type]


# ============================================================
#  Karyogram with Colored Points
# ============================================================
cat("\nCreating  karyogram with annotated points...\n")

# MODIFIED: Swapped width and height for vertical plot
png("cnv_karyogram.png", width=2500, height=2000, res=300)

# MODIFIED: set plot.type=2 for a vertical layout
kp <- plotKaryotype(genome="hg38",
                    plot.type=2,
                    chromosomes=paste0("chr", c(1:22, "X")),
                    main="Chromosomal Distribution of Variants")

# Add CNV markers as points, colored by variant type
# MODIFIED: Using the new color column for plotting
kpPoints(kp, data=cnv_ranges,
         y=0.5,
         col=mcols(cnv_ranges)$color,
         cex=1,
         pch=19)

# NEW: Add a legend to explain the colors
legend("bottomright", legend = names(color_map), fill = color_map, bty = "n", cex = 0.8)

dev.off()
cat("✓ Saved: cnv_karyogram.png\n")


