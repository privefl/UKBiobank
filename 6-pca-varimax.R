
df0 <- readRDS("pheno.rds")

POP <- df0$pop
H <- df0$height
levels(POP) <- c(NA, NA,"White",NA,"Asian","Black","Asian",NA,
                 "White","White","White","White and Black","White and Black",
                 "White and Asian",NA,"Indian","Indian","Indian","Asian (other)",
                 "Black","Black","Black")

pca <- readRDS("UKB_SVD_auto.rds")

system.time(
  vari <- varimax(pca$v)
) # 163 sec

U <- predict(pca)
U2 <- U %*% vari$rotmat

library(ggplot2)
plot2 <- function(pcs = 1:2, ...) {
  ind <- c(which(POP != "White"), sample(which(POP == "White"), 10e3))
  ind <- sample(ind, 10e3)
  p1 <- qplot(U[ind, pcs[1]], U[ind, pcs[2]], color = POP[ind], alpha = I(0.5), ...) +
    theme_bigstatsr()
  p2 <- qplot(U2[ind, pcs[1]], U2[ind, pcs[2]], color = POP[ind], alpha = I(0.5), ...) +
    theme_bigstatsr()
  cowplot::plot_grid(p1, p2, ncol = 1)
}

plot(pca)
plot2()
plot2(3:4)
plot2(5:6)
plot2(7:8)
plot2(9:10)
plot2(11:12)
plot2(13:14)
plot2(15:16)
plot2(17:18)
plot2(19:20)

df <- df0
df$U <- I(U)
df$U2 <- I(U2)
df$absU <- I(abs(U))
df$absU2 <- I(abs(U2))
df$sqU <- I(U^2)
df$sqU2 <- I(U2^2)

ind <- with(df, which(height > 140 & height < 200 & !is.na(deprivation_index)))

summary(mylm0 <- lm(height ~ sex + date, data = df[ind, ]))
# R2 = 0.5155
# Residuals:
#     Min      1Q  Median      3Q     Max
# -33.415  -4.325  -0.035   4.275  34.968



summary(mylm2 <- lm(height ~ sex + date + pop + deprivation_index, data = df[ind, ]))
# 0.5328
# Residuals:
#     Min      1Q  Median      3Q     Max
# -32.569  -4.271  -0.084   4.190  35.058

library(dplyr)
group_by(df, round(deprivation_index)) %>%
  summarise(n = sum(!is.na(height)),
            height = mean(height, na.rm = TRUE))


summary(mylm <- lm(height ~ sex + date + U + absU + absU2, data = df[ind, ]))
# R2 = 0.538

plot2(c(3, 5))

plot3 <- function(pcs = 1:2, ...) {
  ind <- c(which(POP != "White"), sample(which(POP == "White"), 10e3))
  ind <- sample(ind[(H2 > 150 & H2 < 190)[ind]], 10e3)
  p1 <- qplot(U[ind, pcs[1]], U[ind, pcs[2]], color = H2[ind], alpha = I(0.5), ...) +
    theme_bigstatsr() +
    scale_color_viridis_c()
  p2 <- qplot(U2[ind, pcs[1]], U2[ind, pcs[2]], color = H2[ind], alpha = I(0.5), ...) +
    theme_bigstatsr() +
    scale_color_viridis_c()
  cowplot::plot_grid(p1, p2, ncol = 1)
}

plot3()
plot3(3:4)
plot3(5:6)
plot3(7:8)
plot3(9:10)
plot3(11:12)
plot3(13:14)
plot3(15:16)
plot3(17:18)
plot3(19:20)


plot3(c(3, 5))


summary(mylm <- lm(height ~ sex + date + U + absU + sqU, data = df))
# 0.5382
summary(mylm2 <- lm(height ~ sex + date + U2 + absU2 + sqU2, data = df))
# 0.538

as.matrix(df)
plot3(c(1, 7))
plot2(c(1, 7))
