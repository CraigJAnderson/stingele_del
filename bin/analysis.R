library(data.table)
library(see)
library(parameters)

##input is just deletion 5' end locations in a single column
files <- list(
 lineA.con = "AAVS1_no_U.del",
 lineH.con = "HMCES_KO_no_U.del",
 lineA.24  = "AAVS1_24.del",
 lineH.24  = "HMCES_KO_24.del",
 lineA.31  = "AAVS1_31.del",
 lineH.31  = "HMCES_KO_31.del"
)

##classification of sites
abasic_24 <- c(79, 81)
abasic_31 <- c(72, 74)
micro_hom <- c(101, 110)

##use sample names to add clasification to each mutation
df <- rbindlist(lapply(names(files), function(sample_name) {
 strain_construct <- strsplit(sample_name, "\\.")[[1]]
 vec <- fread(files[[sample_name]], header = T)$V5
 data.table(strain = strain_construct[1], construct = strain_construct[2], endpoint = vec)
}))

##function for getting stats
calc_summary <- function(dat) {
 WT <- dat[strain == "lineA"]
 KO <- dat[strain == "lineH"]
 pWT <- WT$abasic / (WT$abasic + WT$microhomology)
 pKO <- KO$abasic / (KO$abasic + KO$microhomology)
 ##format nicely
 results <- matrix(c(WT$abasic, WT$microhomology, KO$abasic, KO$microhomology), nrow = 2, byrow = TRUE, dimnames = list(strain = c("WT", "KO"), endpoint = c("abasic", "microhomology"))  )
 ft <- fisher.test(results)
 data.table(
  WT_proportion   = pWT,
  KO_proportion   = pKO,
  Risk_difference = pKO - pWT,
  Relative_risk   = pKO / pWT,
  Odds_ratio      = (KO$abasic / KO$microhomology) / (WT$abasic / WT$microhomology),
  Fisher_OR       = unname(ft$estimate),
  Fisher_CI_low   = ft$conf.int[1],
  Fisher_CI_high  = ft$conf.int[2],
  Fisher_p        = ft$p.value
 )
}

abasic_analysis <- function(df_all, abasic_range, target_construct) {
 d <- df_all[construct %in% c("con", target_construct)]
 d[, feature := fcase(between(endpoint, micro_hom[1], micro_hom[2]), "microhomology", between(endpoint, abasic_range[1], abasic_range[2]), "abasic")]
 d[, `:=`(feature   = relevel(factor(feature), ref = "microhomology"), construct = relevel(factor(construct), ref = "con"), strain    = relevel(factor(strain), ref = "lineA"))]
 counts <- dcast(d[feature %in% c("abasic", "microhomology"), .N, by = .(strain, construct, feature)], strain + construct ~ feature, value.var = "N", fill = 0L)
 fit <- glm(cbind(abasic, microhomology) ~ strain * construct, family = binomial, data = counts)
 print(summary(fit))
 print(anova(fit, test = "Chisq"))
 print(model_parameters(fit, exponentiate = TRUE, ci_method = "wald"))
 summary_tab <- counts[, calc_summary(.SD), by = construct]
 list(counts = counts, fit = fit, summary = summary_tab)
}

res24 <- abasic_analysis(df, abasic_24, "24")
res31 <- abasic_analysis(df, abasic_31, "31")

counts24 <- res24$counts;  fit24 <- res24$fit;  summary24 <- res24$summary
counts31 <- res31$counts;  fit31 <- res31$fit;  summary31 <- res31$summary
