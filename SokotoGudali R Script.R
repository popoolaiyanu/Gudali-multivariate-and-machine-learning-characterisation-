# =====================================================================
# Sexual dimorphism in body conformation of Sokoto Gudali (Bokolo) cattle
# Complete, reproducible analysis script
# Input : Bokolo_clean.xlsx (sheet "Clean_data"; 159 Sokoto Gudali cattle,
#         61 females and 98 males, after quality control)
# R version 4.5.0
# =====================================================================

# setwd("path/to/your/folder")   # <-- edit

pkgs <- c("readxl", "psych", "corrplot", "caret", "randomForest", "kernlab",
          "MASS", "ggplot2", "ggrepel", "patchwork")
# install.packages(pkgs)   # run once if needed
invisible(lapply(pkgs, library, character.only = TRUE))

# =====================================================================
# 1. DATA (cleaned dataset: 159 animals)
# =====================================================================
bok <- as.data.frame(read_excel("Bokolo_clean.xlsx", sheet = "Clean_data"))
names(bok) <- make.names(names(bok))
bok$Sex <- factor(bok$Sex, levels = c("F", "M"))

dim(bok)                 # 159 rows
table(bok$Sex)           # 61 F, 98 M
colSums(is.na(bok))

# =====================================================================
# 2. TRAITS
# =====================================================================
core <- c("Heart.girth", "Body.length", "Height.at.withers", "Hip.height", "Rump.height",
          "Rump.length", "Rump.width", "Thoracic.length", "Shoulder.width",
          "Hump.length", "Neck.length", "Neck.circumference", "Cannon.circumference",
          "Fore.leg.length", "Hind.leg.length", "Tail.length",
          "Head.length", "Head.width", "Distance.btw.eyes...nose", "Ear.length", "Ear.width")
nice <- c("Heart girth", "Body length", "Height at withers", "Hip height", "Rump height",
          "Rump length", "Rump width", "Thoracic length", "Shoulder width",
          "Hump length", "Neck length", "Neck circumference", "Cannon circumference",
          "Foreleg length", "Hind leg length", "Tail length",
          "Head length", "Head width", "Eye-muzzle distance", "Ear length", "Ear width")
tested <- c(core, "Horn.length")                # 22 traits tested for sex differences

# =====================================================================
# 3. DESCRIPTIVES AND SEX DIFFERENCES (Table 2)
# =====================================================================
describe_trait <- function(x) {
  x <- x[!is.na(x)]
  c(n = length(x), Mean = mean(x), SD = sd(x), SE = sd(x) / sqrt(length(x)),
    Min = min(x), Max = max(x), CV = 100 * sd(x) / mean(x))
}
desc_all <- t(sapply(bok[tested], describe_trait))
desc_F   <- t(sapply(bok[bok$Sex == "F", tested], describe_trait))
desc_M   <- t(sapply(bok[bok$Sex == "M", tested], describe_trait))

# Raw difference (Welch t-test), BH-adjusted across the 22 traits
p_raw <- sapply(tested, function(v) t.test(bok[bok$Sex == "F", v], bok[bok$Sex == "M", v])$p.value)
p_raw_adj <- p.adjust(p_raw, method = "BH")

# Size-adjusted difference: trait ~ Heart girth + Sex, BH-adjusted across 21 traits
adj_traits <- setdiff(tested, "Heart.girth")
size_adj <- t(sapply(adj_traits, function(v) {
  d <- bok[!is.na(bok[[v]]), ]
  s <- summary(lm(d[[v]] ~ Heart.girth + Sex, data = d))$coefficients
  c(Diff_MF = s["SexM", "Estimate"], SE = s["SexM", "Std. Error"], p = s["SexM", "Pr(>|t|)"])
}))
size_adj <- cbind(size_adj, p_adj = p.adjust(size_adj[, "p"], method = "BH"))
idx <- match(tested, rownames(size_adj))

fmt <- function(m, se) sprintf("%.2f ± %.2f", m, se)
table2 <- data.frame(
  Trait = tested, n = desc_all[, "n"],
  Overall = fmt(desc_all[, "Mean"], desc_all[, "SE"]), CV = round(desc_all[, "CV"], 1),
  Female = fmt(desc_F[, "Mean"], desc_F[, "SE"]), Male = fmt(desc_M[, "Mean"], desc_M[, "SE"]),
  P_raw_adj = round(p_raw_adj, 3),
  SizeAdj_Diff_MF = round(size_adj[idx, "Diff_MF"], 2), SizeAdj_SE = round(size_adj[idx, "SE"], 2),
  SizeAdj_P_adj = round(size_adj[idx, "p_adj"], 3))
table2
write.csv(table2, "Table2_descriptives_by_sex.csv", row.names = FALSE)

# Sex-specific traits
describe_trait(bok$scrotal.circumference[bok$Sex == "M"])
describe_trait(bok$Udder.diameter[bok$Sex == "F"])

# =====================================================================
# 4. CORRELATIONS WITHIN SEX (Figure 1)
# =====================================================================
cor_F <- cor(bok[bok$Sex == "F", core], use = "pairwise.complete.obs")
cor_M <- cor(bok[bok$Sex == "M", core], use = "pairwise.complete.obs")
dimnames(cor_F) <- dimnames(cor_M) <- list(nice, nice)
cor_FM <- cor_F; cor_FM[upper.tri(cor_FM)] <- cor_M[upper.tri(cor_M)]; diag(cor_FM) <- 1

png("Figure1_correlations_by_sex.png", width = 2400, height = 2200, res = 300)
corrplot(cor_FM, method = "color", type = "full", tl.col = "black", tl.cex = 0.7,
         addCoef.col = "black", number.cex = 0.4, diag = FALSE,
         col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200))
dev.off()
write.csv(round(cor_FM, 2), "TableS_correlations_F_lower_M_upper.csv")

# =====================================================================
# 5. PRINCIPAL COMPONENT ANALYSIS (Table 3, Figure 2, Figure S1)
# =====================================================================
pca_data <- na.omit(bok[, c("Sex", core)])
nrow(pca_data); table(pca_data$Sex)             # 152: 56 F, 96 M
X <- scale(pca_data[, core]); colnames(X) <- nice

KMO(X)$MSA
cortest.bartlett(cor(X), n = nrow(X))

pca0 <- prcomp(X)
eig <- pca0$sdev^2
var_tab <- data.frame(PC = paste0("PC", seq_along(eig)), Eigenvalue = round(eig, 3),
                      Variance_pct = round(100 * eig / sum(eig), 2),
                      Cumulative_pct = round(cumsum(100 * eig / sum(eig)), 2))
var_tab[1:5, ]
fa.parallel(X, fa = "pc", n.iter = 100, show.legend = FALSE, main = "")

png("FigureS1_scree.png", width = 1800, height = 1400, res = 300)
plot(eig[1:10], type = "b", pch = 19, xlab = "Principal component", ylab = "Eigenvalue")
abline(h = 1, lty = 2, col = "red")
dev.off()

load3 <- pca0$rotation[, 1:3] %*% diag(pca0$sdev[1:3]); colnames(load3) <- c("PC1", "PC2", "PC3")
scores <- data.frame(Sex = pca_data$Sex, pca0$x[, 1:3])
# orient PC1 so that larger animals have positive scores
if (mean(load3[, "PC1"]) < 0) { load3[, "PC1"] <- -load3[, "PC1"]; scores$PC1 <- -scores$PC1 }
round(load3[order(-load3[, "PC1"]), ], 2)
write.csv(round(load3, 3), "Table3_PCA_loadings.csv")
write.csv(var_tab, "Table3_PCA_variance.csv", row.names = FALSE)

sapply(c("PC1", "PC2", "PC3"), function(pc)
  round(c(Mean_F = mean(scores[scores$Sex == "F", pc]), Mean_M = mean(scores[scores$Sex == "M", pc]),
          p = t.test(scores[[pc]] ~ scores$Sex)$p.value), 3))

scores$Sex <- factor(scores$Sex, levels = c("F", "M"), labels = c("Female", "Male"))
pA <- ggplot(scores, aes(PC1, PC2, colour = Sex, shape = Sex)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_point(size = 1.8, alpha = 0.8) + stat_ellipse(level = 0.95, linewidth = 0.6) +
  scale_colour_manual(values = c(Female = "#B2182B", Male = "#2166AC")) +
  scale_shape_manual(values = c(Female = 17, Male = 19)) +
  labs(x = sprintf("PC1 - body size (%.1f%%)", var_tab$Variance_pct[1]),
       y = sprintf("PC2 - forequarter vs hindquarter (%.1f%%)", var_tab$Variance_pct[2]), title = "(a)") +
  theme_classic(base_size = 10) + theme(legend.position = c(0.12, 0.9), legend.title = element_blank())
ld <- data.frame(Trait = rownames(load3), PC1 = load3[, "PC1"], PC2 = load3[, "PC2"])
pB <- ggplot(ld, aes(PC1, PC2)) +
  geom_hline(yintercept = 0, linetype = 3, colour = "grey60") +
  geom_vline(xintercept = 0, linetype = 3, colour = "grey60") +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.15, "cm")), colour = "grey40") +
  geom_text_repel(aes(label = Trait), size = 2.6, max.overlaps = Inf, segment.colour = "grey70") +
  coord_cartesian(xlim = c(0, 1.05), ylim = c(-0.65, 0.5)) +
  labs(x = "Loading on PC1", y = "Loading on PC2", title = "(b)") + theme_classic(base_size = 10)
ggsave("Figure2_PCA.png", pA + pB, width = 10, height = 4.8, dpi = 300)

# =====================================================================
# 6. MACHINE-LEARNING CLASSIFICATION OF SEX (Table 4, Figure 3)
# =====================================================================
ml_data <- pca_data
ml_data$Sex <- factor(ml_data$Sex, levels = c("F", "M"), labels = c("Female", "Male"))
ctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10, classProbs = TRUE,
                     summaryFunction = twoClassSummary, savePredictions = "final")

set.seed(2026); m_lda <- train(Sex ~ ., data = ml_data, method = "lda",
                               preProcess = c("center", "scale"), trControl = ctrl, metric = "ROC")
set.seed(2026); m_glm <- train(Sex ~ ., data = ml_data, method = "glm", family = "binomial",
                               preProcess = c("center", "scale"), trControl = ctrl, metric = "ROC")
set.seed(2026); m_rf  <- train(Sex ~ ., data = ml_data, method = "rf", ntree = 500,
                               tuneLength = 5, trControl = ctrl, metric = "ROC")
set.seed(2026); m_svm <- train(Sex ~ ., data = ml_data, method = "svmRadial",
                               preProcess = c("center", "scale"), tuneLength = 5,
                               trControl = ctrl, metric = "ROC")
m_rf$bestTune; m_svm$bestTune

res <- resamples(list(LDA = m_lda, Logistic = m_glm, RF = m_rf, SVM = m_svm))
summary(res)
summary(diff(res, metric = "ROC"))              # Bonferroni-adjusted paired comparisons

ml_tab <- do.call(rbind, lapply(names(res$methods), function(m) {
  r <- res$values[, grep(paste0("^", m, "~"), names(res$values))]
  names(r) <- sub(".*~", "", names(r)); ba <- (r$Sens + r$Spec) / 2
  data.frame(Model = m, AUC = sprintf("%.3f ± %.3f", mean(r$ROC), sd(r$ROC)),
             Sens_F = sprintf("%.3f ± %.3f", mean(r$Sens), sd(r$Sens)),
             Spec_M = sprintf("%.3f ± %.3f", mean(r$Spec), sd(r$Spec)),
             Bal_Acc = sprintf("%.3f ± %.3f", mean(ba), sd(ba)))
}))
ml_tab
write.csv(ml_tab, "Table4_ML_sex_classification.csv", row.names = FALSE)

imp_rf <- varImp(m_rf)$importance
imp_df <- data.frame(Trait = nice[match(rownames(imp_rf), core)], Importance = imp_rf[, 1])
k <- match(core[match(imp_df$Trait, nice)], rownames(size_adj))
imp_df$Larger_in <- ifelse(is.na(k) | size_adj[k, "p_adj"] >= 0.05, "No significant difference",
                           ifelse(size_adj[k, "Diff_MF"] > 0, "Larger in males", "Larger in females"))
p_imp <- ggplot(imp_df, aes(x = reorder(Trait, Importance), y = Importance, fill = Larger_in)) +
  geom_col(width = 0.7) + coord_flip() +
  scale_fill_manual(values = c("Larger in males" = "#2166AC", "Larger in females" = "#B2182B",
                               "No significant difference" = "grey70")) +
  labs(x = NULL, y = "Random Forest importance (scaled 0-100)", fill = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = c(0.72, 0.15))
ggsave("Figure3_RF_importance_sex.png", p_imp, width = 6, height = 5.5, dpi = 300)

sessionInfo()   # package versions for the Methods
