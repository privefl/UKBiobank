library(bigsnpr)
library(bigreadr)
library(dplyr)
library(ggplot2)

rel <- fread2("data/ukb25589_rel_s488346.dat")
fam <- fread2("data/ukbb_bed/ukbb_488282.fam")
is_rel <- (fam$V2 %in% rel$ID1 | fam$V2 %in% rel$ID2)
obj.bed <- bed("data/ukbb_bed/ukbb_488282.bed")

## Self-reported ancestry (https://biobank.ctsu.ox.ac.uk/crystal/coding.cgi?id=1001)
code_ancestry <- fread2("data/coding1001.tsv")
csv <- "data/ukb22544.csv"
df0 <- fread2(csv, select = c("eid", "21000-0.0", "22006-0.0"),
              col.names = c("eid", "pop", "is_caucasian")) %>%
  mutate(
    pop  = factor(pop, levels = code_ancestry$coding,
                  labels = code_ancestry$meaning),
    is_caucasian = as.logical(is_caucasian)
  )
pop_UKBB <- df0$pop[match(obj.bed$fam$sample.ID, df0$eid)]
table(pop_UKBB)

ind_used <- c(sample(which(pop_UKBB == "British"), 10e3),
              sample(which(pop_UKBB == "Irish"), 5e3),
              which(!pop_UKBB %in% c("British", "Irish")))
ind.train <- intersect(ind_used, which(!is_rel))
table(pop_UKBB[ind.train])

system.time(
  obj.svd <- bed_autoSVD(obj.bed, ind.row = ind.train,
                         k = 50, ncores = nb_cores())
) # 2.9H
# saveRDS(obj.svd, "tmp-results/SVD_UKBB_restricted.rds")

system.time(
  proj <- bed_projectSelfPCA(obj.svd, obj.bed,
                             ind.row = rows_along(obj.bed)[-ind.train],
                             ncores = nb_cores())
) # 21 min
PC <- matrix(NA_real_, nrow(obj.bed), ncol(obj.svd$u))
PC[ind.train, ]  <- predict(obj.svd)
PC[-ind.train, ] <- proj$OADP_proj


#### UMAP ####

knn <- bigutilsr::knn_parallel(PC, k = 16, ncores = nb_cores())

system.time(
  test2 <- uwot::umap(
    PC, n_epochs = 500, spread = 10, min_dist = 1, bandwidth = 10,
    n_threads = nb_cores(), n_sgd_threads = nb_cores(), scale = FALSE,
    nn_method = list(idx = knn$nn.idx, dist = knn$nn.dists)
  )
)
#     user   system  elapsed
# 3501.513    6.765  361.080

ind <- sample(nrow(PC), 5e3)
qplot(test2[ind, 1], test2[ind, 2], size = I(2), color = pop_UKBB[ind]) +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Population")
plotly::ggplotly()

dist <- bigutilsr::covRob(PC, estim = "pairwiseGK")$dist

qplot(test2[ind, 1], test2[ind, 2], size = I(2), color = dist[ind]) +
  theme_bigstatsr() +
  scale_color_viridis_c(trans = "log") +
  # theme(legend.position = "none") +
  labs(color = "Population")
plotly::ggplotly()

ind <- sample(nrow(PC), 50e3)
qplot(test2[ind, 1], test2[ind, 2], size = I(2), color = dist[ind]) +
  scale_color_viridis_c(trans = "log") +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Distance")
qplot(test2[ind, 1], test2[ind, 2], size = I(2),
      color = as.factor(obj.bed$fam$sex[ind.row[ind]])) +
  scale_color_viridis_d() +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Distance")


knn2 <- bigutilsr::knn_parallel(PC[, 1:15], k = 16, ncores = nb_cores())

system.time(
  test3 <- uwot::umap(
    PC[, 1:15], n_epochs = 200, spread = 10, min_dist = 1, bandwidth = 10,
    n_threads = nb_cores(), n_sgd_threads = nb_cores(), scale = FALSE,
    nn_method = list(idx = knn$nn.idx, dist = knn$nn.dists)
  )
)
#     user   system  elapsed
# 3501.513    6.765  361.080

ind <- sample(nrow(PC), 5e3)
qplot(test3[ind, 1], test3[ind, 2], size = I(2), color = pop[ind]) +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Population")
plotly::ggplotly()

dist <- bigutilsr::covRob(PC[, 1:15], estim = "pairwiseGK")$dist


ind <- sample(nrow(PC), 50e3)
qplot(test3[ind, 1], test3[ind, 2], size = I(2), color = dist[ind]) +
  scale_color_viridis_c(trans = "log") +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Distance")

qplot(test3[ind, 1], test3[ind, 2], size = I(2), color = log(dist)[ind] > 4) +
  scale_color_viridis_d() +
  theme_bigstatsr() +
  theme(legend.position = "none")
qplot(test2[ind, 1], test2[ind, 2], size = I(2),
      color = as.factor(obj.bed$fam$sex[ind.row[ind]])) +
  scale_color_viridis_d() +
  theme_bigstatsr() +
  # theme(legend.position = "none") +
  labs(color = "Distance")
