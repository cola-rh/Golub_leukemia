setwd("/omics/groups/OE0246/internal/guz/cola_hc/examples/Golub_leukemia")
library(cola)

library(golubEsets)
data(Golub_Merge)
m = exprs(Golub_Merge)
colnames(m) = paste0("sample_", colnames(m))
anno = pData(Golub_Merge)

m[m <= 1] = NA
m = log10(m)

m = adjust_matrix(m)

library(preprocessCore)
cn = colnames(m)
rn = rownames(m)
m = normalize.quantiles(m)
colnames(m) = cn
rownames(m) = rn

rh = hierarchical_partition(m, cores = 4, 
    anno = anno[, c("ALL.AML"), drop = FALSE],
    anno_col = c("ALL" = "red", "AML" = "blue"))
saveRDS(rh, file = "Golub_leukemia_cola_rh.rds")

cola_report(rh, output = "Golub_leukemia_cola_rh_report", title = "cola Report for Hierarchical Partitioning - 'Golub_leukemia'")
