library(susieR)
set.seed(1)
n    <- 1000
p    <- 500
beta <- rep(0,p)
beta[c(1,200,300,400)] <- 1
X   <- matrix(rnorm(n*p),nrow=n,ncol=p)
X[, 200] = X[, 1]
X[, 300] = -X[, 1]
#X[,1]
#X[,2]
#X[,200]
#X[,300]
y   <- X %*% beta + rnorm(n)
fitted <- susie(X,y,L=10, standardize=FALSE)
dir='susie_test'
write.table(fitted$alpha, paste(dir, '/alpha.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # get pips per va
write.table(fitted$lbf, paste(dir, '/lbf.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #log bayes factor of 
write.table(fitted$lbf_variable, paste(dir, '/lbf_variable.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) #lb
write.table(fitted$sigma2, paste(dir, '/sigma2.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated res
write.table(fitted$V, paste(dir, '/V.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # estimated prior varianc
write.table(fitted$converged, paste(dir, '/converged.txt', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # TRUE or
write.table(susie_get_lfsr(fitted), paste(dir, '/lfsr.tab', sep=''), sep='\t', row.names=FALSE, col.names=FALSE) # local 

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


pdf("susie_test/betas.pdf")
plot(coef(fitted)[-1],pch = 20)
dev.off() 

