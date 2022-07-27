library(susieR)
library(rhdf5)
library(argparse)

parser = ArgumentParser()
parser$add_argument('outdir')
parser$add_argument('pheno_residuals')
parser$add_argument('gt_residuals')
parser$add_argument('L', type='integer')
args = parser$parse_args()

dir = args$outdir
print(args)

pheno_residuals = as.vector(h5read(args$pheno_residuals, 'pheno_residuals'))
gt_residuals = t(h5read(args$gt_residuals, 'gt_residuals'))

coverage=0.9
fitted = susie(
  gt_residuals,
  pheno_residuals,
  L=args$L,
  min_abs_corr = 0.8,
  coverage=coverage,
  verbose=TRUE
)

write.table(susie_get_pip(fitted, prune_by_cs=TRUE), paste(dir, '/pruned_pips', sep=''), sep='\t', row.names=FALSE, col.names=FALSE)
write.table(fitted$pip, paste(dir, '/pip.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE)
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
