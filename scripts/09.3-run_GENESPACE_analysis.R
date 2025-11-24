library (GENESPACE)

args <- commandArgs(trailingOnly = TRUE)
# get the folder where the genespace workingDirectory is located (passed from the bash script)
wd <- args[1]

# Initialize GENESPACE parameters
# gpar <- init_genespace(wd = wd, path2mcscanx = "/data/courses/assembly-annotation-course/CDS_annotation/softwares/MCScanX")
gpar <- init_genespace(wd = wd, path2mcscanx = "/data/courses/assembly-annotation-course/CDS_annotation/softwares/MCScanX", nCores = 20, verbose = TRUE)

# Run the GENESPACE pipeline
# This performs all-against-all BLAST and synteny analysis (collinearity detection)
out <- run_genespace(gpar, overwrite = TRUE)
pangenome <- query_pangenes(
  out,
  bed = NULL,
  refGenome = "TAIR10", # Use TAIR10 as the reference genome for pangenome calculation
  transform = TRUE,
  showArrayMem = TRUE,
  showNSOrtho = TRUE,
  maxMem2Show = Inf
)

# Save the pangenome object as an R data file (.rds) for future use
saveRDS(pangenome, file = file.path(wd, "pangenome_matrix.rds"))
# in your next script, you can load the pangenome matrix with: pangenome <- readRDS(file.path(wd, "pangenome_matrix.rds")) and then use it for downstream analyses, e.g., calculating core, accessory and specific genes