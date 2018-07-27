# https://biobank.ctsu.ox.ac.uk/crystal/docs/ukbgene_instruct.html
for (chr in 1:22) {
  system(glue::glue("./ukbgene cal -c{chr}"))
  system(glue::glue("./ukbgene cal -c{chr} -m "))
}


# list with bed / bim / fam files to be merged
files <- glue::glue(
  "ukb_cal_chr{chr}_v2.bed",
  "ukb_snp_bim/ukb_snp_chr{chr}_v2.bim",
  "ukb{ID1}_cal_chr{chr}_v2_{ID2}.fam",
  .sep = "\t", chr = 1:22, ID1 = 25589, ID2 = "s488363"
)
writeLines(files, "merge-bed-files.txt")

# Related individuals to be removed
system("./ukbgene rel")
rel <- data.table::fread("ukb25589_rel_s488346.dat", data.table = FALSE)

bigsnpr:::write.table2(subset(rel, Kinship > 0.08, c(1, 1)), "rm-rel-indiv.txt")

library(bigsnpr)
plink <- download_plink(".")
system(glue::glue(
  "{plink} --merge-list merge-bed-files.txt",
  " --make-bed --out UKB-merged",
  " --maf 0.02",
  " --mind 0.1",
  " --geno 0.1",
  " --hwe 1e-50",
  " --remove rm-rel-indiv.txt"
))

system.time(
  rds <- snp_readBed("UKB-merged.bed")
) # 17.4 min

system.time(
  ukb <- snp_attach(rds)
)
pryr::object_size(ukb) # 57.2 MB

