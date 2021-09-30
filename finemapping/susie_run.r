library(reticulate)
library(susieR)

np = import("numpy")

ukb = Sys.getenv("UKB")
args = commandArgs(trailingOnly=TRUE)

phenotype = args[1]
chrom = args[2]
start = args[3]
end = args[4]

dir = paste(
    ukb, '/finemapping/susie_results/', phenotype, '/', chrom, '_', start, '_', end,
sep='')

pheno_residuals = np$load(paste(dir, '/pheno_residuals.npy', sep=''))
gt_residuals = np$load(paste(dir, '/gt_residuals.npy', sep=''))

fitted = susie(gt_residuals, pheno_residuals, L=10)
save(fitted, file = paste(dir, '/susie_fit.RData', sep=''))

