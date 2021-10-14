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

coverage=0.9
fitted = susie(gt_residuals, pheno_residuals, L=10, scaled_prior_variance = 0.005, min_abs_corr = 0, coverage=coverage)
#save(fitted, file = paste(dir, '/susie_fit.RData', sep=''))

write.table(fitted$alpha, paste(dir, '/alpha.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # get pips per variable by 1-apply(1-alpha, 1, prod)
write.table(fitted$lbf, paste(dir, '/lbf.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #log bayes factor of each cs
write.table(fitted$lbf_variable, paste(dir, '/lbf_variable.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #lbf of each variable per cs
write.table(fitted$sigma2, paste(dir, '/sigma2.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated residual posterior variance of y
write.table(fitted$V, paste(dir, '/V.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated prior variance of each CS
write.table(fitted$converged, paste(dir, '/converged.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # TRUE or FALSE, did SuSiE converge within the specified tolerance
write.table(susie_get_lfsr(fitted), paste(dir, '/lfsr.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # local false sign rate for each cs - probability that we've got the sign of the effect wrong

sink(paste(dir, '/requested_coverage.txt', sep=''))
cat(fitted$sets$requested_coverage, '\n')
sink()
count = 0
for (id in fitted$sets$cs_index)  {
  count = count+1
  cs_id = paste('L', id, sep='')
  sink(paste(dir, '/cs', id, '.txt', sep=''))
  cat(fitted$sets$cs[[cs_id]], sep=' ', '\n')
  cat(fitted$sets$coverage[count], '\n')
  cat(fitted$sets$purity[cs_id, 'min.abs.corr'],
      fitted$sets$purity[cs_id, 'mean.abs.corr'],
      fitted$sets$purity[cs_id, 'median.abs.corr'], '\n')
  sink()
}
