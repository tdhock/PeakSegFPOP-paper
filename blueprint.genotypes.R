source("packages.R")

geno.txt.vec <- Sys.glob("/mnt/abaproject/Tony_group/wcheung/genotypes_processed/EMC_BluePrint_*.txt")
names(geno.txt.vec) <- sub("_.*", "", sub(".*BluePrint_", "", geno.txt.vec))

genotype.mat <- NULL
for(geno.txt.i in seq_along(geno.txt.vec)){
  geno.txt <- geno.txt.vec[[geno.txt.i]]
  cat(sprintf("%4d / %4d %s\n", geno.txt.i, length(geno.txt.vec), geno.txt))
  person <- names(geno.txt.vec)[[geno.txt.i]]
  geno <- fread(geno.txt)
  if(is.null(genotype.mat)){
    genotype.mat <- matrix(
      NA_character_, nrow(geno), length(geno.txt.vec),
      dimnames=list(
        "Chr:Position"=geno[["Chr:Position"]],
        person=names(geno.txt.vec)))
  }
  mat.has.pos <- geno[["Chr:Position"]] %in% rownames(genotype.mat)
  some.geno <- geno[mat.has.pos]
  genotype.mat[some.geno[["Chr:Position"]], geno.txt.i] <-
    some.geno[, paste0(A1, A2)]
}

sum(genotype.mat=="GA", na.rm=TRUE)

saveRDS(genotype.mat, "blueprint.genotypes.rds")
write.matrix(genotype.mat, "blueprint.genotypes.txt")

sample.info <- fread("blueprint_sample_info.csv")
sample.info[, has.genotype := person %in% names(geno.txt.vec)]
sample.info[has.genotype==TRUE & center != "McGill", table(experiment, cellType)]
