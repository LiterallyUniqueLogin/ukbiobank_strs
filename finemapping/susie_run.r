library(susieR)
library(rhdf5)
library(argparse)

parser = ArgumentParser()
parser$add_argument('outdir')
parser$add_argument('pheno_residuals')
parser$add_argument('gt_residuals')
parser$add_argument('L', type='integer')
parser$add_argument('max_iter', type='integer')
parser$add_argument('--snp-p-over-str-p', type='double', default=NULL)
parser$add_argument('--varnames-file')
parser$add_argument('--tol', type='double', default=0.001)
parser$add_argument('--scaled-prior-variance', type='double', default=0.005)
parser$add_argument('--residual-variance', type='double')
args = parser$parse_args()

dir = args$outdir
print(args)

pheno_residuals = as.vector(h5read(args$pheno_residuals, 'pheno_residuals'))
gt_residuals = t(h5read(args$gt_residuals, 'gt_residuals'))

stopifnot(is.null(args$snp_p_over_str_p) == is.null(args$varnames_file))

if (!is.null(args$snp_p_over_str_p)) {
  str_poses = c()
  # from https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line
  con = file(args$varnames_file, 'r')
  while ( TRUE ) {
    line = readLines(con, n=1)
    if (length(line) == 0) {
      break
    }
    str_poses = c(str_poses, startsWith(line, 'STR'))
  }
  prior_weights = double(length(str_poses))
  prior_weights[str_poses] = 1
  prior_weights[!str_poses] = args$snp_p_over_str_p
  prior_weights = prior_weights/sum(prior_weights)
  print('prior weights')
  print(prior_weights)
} else {
  prior_weights = NULL
}

if (!is.null(args$residual_variance)) {
  residual_variance = var(pheno_residuals)*args$residual_variance
  print('residual_variance')
  print(residual_variance)
} else {
  residual_variance = NULL
}

coverage=0.9
fitted = susie(
  gt_residuals,
  pheno_residuals,
  L=args$L,
  max_iter=args$max_iter,
  prior_weights = prior_weights,
  scaled_prior_variance = args$scaled_prior_variance,
  tol = args$tol,
  residual_variance = residual_variance,
  min_abs_corr = 0,
  coverage=coverage,
  verbose=TRUE
)

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
