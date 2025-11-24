#-------------------------------------------------
# WARNING: As this script uses the final gene annotation,
# make sure to run it only after completing the gene annotation steps; 
# i.e., after script 06.3-filter_refine_annotation.sh has been executed.
# Otherwise, it is also possible to visualise only TE densities without genes by
# commenting out the gene-related parts of the code.
#-------------------------------------------------

# Load the circlize package
library(circlize)
library(tidyverse)
library(ComplexHeatmap)

# Create plot directory if it doesn't exist
plot_dir <- "plots"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}

# Load the TE annotation GFF3 file
gff_file <- "Elh-2.asm.bp.p_ctg.fa.mod.EDTA.TEanno.gff3"
gff_data <- read.table(gff_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Load gene annotation file
gene_annotation_file <- "filtered.genes.renamed.gff3"
gene_annotation <- read.table(gene_annotation_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Swap coordinates where end larger than start
gene_annotation <- gene_annotation %>%
  mutate(
    start = pmin(V4, V5),
    end   = pmax(V4, V5)
  )

# Check the superfamilies present in the GFF3 file, and their counts
gff_data$V3 %>% table()

# custom ideogram data
custom_ideogram <- read.table("Elh-2.asm.bp.p_ctg.fa.fai", header = FALSE, stringsAsFactors = FALSE)
custom_ideogram$chr <- custom_ideogram$V1
custom_ideogram$start <- 1
custom_ideogram$end <- custom_ideogram$V2
custom_ideogram <- custom_ideogram[, c("chr", "start", "end")]
custom_ideogram <- custom_ideogram[order(custom_ideogram$end, decreasing = T), ]
sum(custom_ideogram$end[1:20])

# Select only the first 15 longest scaffolds
custom_ideogram <- custom_ideogram[1:15, ]

# Function to filter GFF3 data based on Superfamily
filter_superfamily <- function(gff_data, superfamily, custom_ideogram) {
  filtered_data <- gff_data[gff_data$V3 == superfamily, ] %>%
    as.data.frame() %>%
    mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# Function to filter gene data
filter_genes <- function(gene_data, custom_ideogram) {
  filtered_data <- gene_data %>%
    as.data.frame() %>%
    mutate(chrom = V1, start = start, end = end, strand = V7) %>%
    select(chrom, start, end, strand) %>%
    filter(chrom %in% custom_ideogram$chr)
  return(filtered_data)
}

# Save to plot folder - TE density with genes
pdf(file.path(plot_dir, "02.2-TE_density_with_genes.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

# Initialize the circos plot with the custom ideogram
circos.genomicInitialize(custom_ideogram)

# Plot TE density for Gypsy and Copia
circos.genomicDensity(filter_superfamily(gff_data, "Gypsy_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#2E8B57", track.height = 0.07, window.size = 1e5) # SeaGreen
circos.genomicDensity(filter_superfamily(gff_data, "Copia_LTR_retrotransposon", custom_ideogram), 
                      count_by = "number", col = "#DC143C", track.height = 0.07, window.size = 1e5) # Crimson

# Add gene track
circos.genomicDensity(filter_genes(gene_annotation, custom_ideogram), 
                      count_by = "number", col = "black", track.height = 0.07, window.size = 1e5)

lgd <- Legend(
  title = "Annotations", 
  at = c("Gypsy_LTR_retrotransposon", "Copia_LTR_retrotransposon", "Genes"),
  legend_gp = gpar(fill = c("#2E8B57", "#DC143C", "black"))
)
# Center the legend in the middle of the plot
draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("Plot with genes saved to:", file.path(plot_dir, "02.2-TE_density_with_genes.pdf"), "\n")

# Now plot all your most abundant TE superfamilies in one plot with genes

# Get the most abundant superfamilies
superfamily_counts <- gff_data$V3 %>% table() %>% sort(decreasing = TRUE)
top_superfamilies <- names(superfamily_counts)[1:6]  # Top 6 superfamilies

# Define brighter, more sympathetic colors for each superfamily
superfamily_colors <- c(
  "#2E8B57",  # SeaGreen - Gypsy
  "#DC143C",  # Crimson - Copia
  "#1E90FF",  # DodgerBlue
  "#FF8C00",  # DarkOrange
  "#9370DB",  # MediumPurple
  "#20B2AA",  # LightSeaGreen
  "black"     # For genes
)

# Save all superfamilies plot with genes to plot folder
pdf(file.path(plot_dir, "02.2-TE_density_all_superfamilies_with_genes.pdf"), width = 10, height = 10)
gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)

circos.genomicInitialize(custom_ideogram)

# Plot each superfamily
for(i in seq_along(top_superfamilies)) {
  circos.genomicDensity(filter_superfamily(gff_data, top_superfamilies[i], custom_ideogram), 
                        count_by = "number", col = superfamily_colors[i], 
                        track.height = 0.07, window.size = 1e5)
}

# Add gene track
circos.genomicDensity(filter_genes(gene_annotation, custom_ideogram), 
                      count_by = "number", col = "black", track.height = 0.07, window.size = 1e5)

# Add centered legend using ComplexHeatmap for better positioning
lgd_all <- Legend(
  title = "TE Superfamilies & Genes", 
  at = c(top_superfamilies, "Genes"),
  legend_gp = gpar(fill = superfamily_colors[1:(length(top_superfamilies) + 1)]),
  ncol = 2  # Two columns for better layout
)
draw(lgd_all, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))

circos.clear()
dev.off()

cat("All superfamilies with genes plot saved to:", file.path(plot_dir, "02.2-TE_density_all_superfamilies_with_genes.pdf"), "\n")

# Plot the distribution of Athila and CRM clades (known centromeric TEs in Brassicaceae)

# First, let's properly examine the TSV file structure
tsv_file = "Gypsy_sequences.fa.rexdb-plant.cls.tsv"

# Read and examine the raw file
cat("=== DEBUGGING TSV FILE ===\n")
tsv_raw <- readLines(tsv_file, n = 15)
cat("First 15 lines of TSV file:\n")
cat(tsv_raw, sep = "\n")
cat("===================\n")

# Try different reading methods
tsv_data <- read.table(tsv_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, 
                       comment.char = "", quote = "", fill = TRUE)

cat("TSV file dimensions:", dim(tsv_data), "\n")
cat("Column names:", colnames(tsv_data), "\n")
cat("First few rows:\n")
print(head(tsv_data))

# Check the sTE column specifically
if ("sTE" %in% colnames(tsv_data)) {
  cat("First few sTE values:\n")
  print(head(tsv_data$sTE))
  cat("Class of sTE:", class(tsv_data$sTE), "\n")
  cat("Length of sTE:", length(tsv_data$sTE), "\n")
} else {
  cat("sTE column not found. Available columns:", colnames(tsv_data), "\n")
  # If sTE doesn't exist, use the first column
  tsv_data$sTE <- tsv_data[,1]
}

# Extract IDs from GFF data first
cat("\n=== EXTRACTING GFF IDs ===\n")
extract_id <- function(line) {
  id_match <- regmatches(line, regexpr("TE_[0-9]+", line))
  if (length(id_match) > 0) {
    return(id_match[1])
  } else {
    return(NA)
  }
}

gff_data$id <- sapply(gff_data$V9, extract_id)
cat("GFF IDs extracted. Number of non-NA IDs:", sum(!is.na(gff_data$id)), "\n")
cat("Sample GFF IDs:", head(na.omit(gff_data$id)), "\n")

# Now extract IDs from TSV data - use a safer approach
cat("\n=== EXTRACTING TSV IDs ===\n")
tsv_data$id <- NA_character_

for (i in 1:nrow(tsv_data)) {
  te_string <- tsv_data$sTE[i]
  
  # Method 1: Try to extract TE_ID pattern
  id_match <- regmatches(te_string, regexpr("TE_[0-9]+", te_string))
  if (length(id_match) > 0) {
    tsv_data$id[i] <- id_match[1]
  } else {
    # Method 2: Use the first word before space
    first_part <- strsplit(te_string, " ")[[1]][1]
    tsv_data$id[i] <- first_part
  }
}

cat("TSV IDs extracted. Number of non-NA IDs:", sum(!is.na(tsv_data$id)), "\n")
cat("Sample TSV IDs:", head(na.omit(tsv_data$id)), "\n")

# Filter for Athila and CRM clades
if ("Clade" %in% colnames(tsv_data)) {
  clade_data <- tsv_data %>%
    filter(Clade %in% c("Athila", "CRM")) %>%
    select(id, Clade) %>%
    filter(!is.na(id))
  
  cat("\nClade data after filtering:\n")
  print(clade_data)
  cat("Clade counts:\n")
  print(table(clade_data$Clade))
  
  # Merge with GFF data
  if (nrow(clade_data) > 0) {
    new_data <- gff_data %>%
      inner_join(clade_data, by = "id") %>%
      filter(!is.na(Clade))
    
    cat("\nMerged data dimensions:", dim(new_data), "\n")
    cat("Clade counts in merged data:\n")
    print(table(new_data$Clade))
  } else {
    cat("No matching clades found after filtering.\n")
    new_data <- data.frame()
  }
} else {
  cat("Clade column not found in TSV file.\n")
  new_data <- data.frame()
}

# Only proceed if we have data
if (exists("new_data") && nrow(new_data) > 0) {
  filter_clades <- function(clade_data, clade, custom_ideogram) {
    filtered_data <- clade_data[clade_data$Clade == clade, ] %>%
      as.data.frame() %>%
      mutate(chrom = V1, start = V4, end = V5, strand = V6) %>%
      select(chrom, start, end, strand) %>%
      filter(chrom %in% custom_ideogram$chr)
    return(filtered_data)
  }
  
  # Save clades plot with genes
  pdf(file.path(plot_dir, "02.2-TE_density_clades_with_genes.pdf"), width = 10, height = 10)
  gaps <- c(rep(1, length(custom_ideogram$chr) - 1), 5)
  circos.par(start.degree = 90, gap.after = 1, track.margin = c(0, 0), gap.degree = gaps)
  
  circos.genomicInitialize(custom_ideogram)
  
  # Plot available clades
  available_clades <- unique(new_data$Clade)
  clade_colors <- c("CRM" = "#1E90FF", "Athila" = "#DC143C")
  
  for (clade in available_clades) {
    circos.genomicDensity(filter_clades(new_data, clade, custom_ideogram), 
                          count_by = "number", col = clade_colors[clade], 
                          track.height = 0.07, window.size = 1e5)
  }
  
  # Add gene track
  circos.genomicDensity(filter_genes(gene_annotation, custom_ideogram), 
                        count_by = "number", col = "black", track.height = 0.07, window.size = 1e5)
  
  circos.clear()
  
  # Create legend
  lgd <- Legend(
    title = "Clades & Genes", 
    at = c(available_clades, "Genes"),
    legend_gp = gpar(fill = c(clade_colors[available_clades], "black"))
  )
  draw(lgd, x = unit(0.5, "npc"), y = unit(0.5, "npc"), just = c("center"))
  
  dev.off()
  
  cat("Clades plot with genes saved to:", file.path(plot_dir, "02.2-TE_density_clades_with_genes.pdf"), "\n")
} else {
  cat("Skipping clades plot - no data available after merging.\n")
}