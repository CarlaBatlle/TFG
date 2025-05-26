## Counts to TPM

library(IOBR)
setwd("~/Documents/practiques/data")

data <- read.csv("unproc_rnaseq.csv", row.names = 1)
head(data, 10)

data_tpm <- count2tpm(countMat = data, source = "local", idType = "Entrez")
head(data_tpm)

write.csv(data_tpm, "TPM_results.csv", row.names = TRUE)

# BETA-------------
data_BETA <- read.csv("beta_un_merged.csv", row.names = 1)
data_BETA <- as.data.frame(t(data_BETA))
head(data_BETA, 10)

data_BETA_tpm <- count2tpm(countMat = data_BETA, source = "local", idType = "Ensembl")
head(data_BETA_tpm)

write.csv(data_BETA_tpm, "BETA_TPM_results.csv", row.names = TRUE)
