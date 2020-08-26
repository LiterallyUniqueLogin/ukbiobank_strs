#!/usr/bin/env Rscript
library(tidyverse)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
plot_dir = args[1]
chrom = args[2]
pos = args[3]
phenotype = args[4]

unit = readLines(sprintf("%s/info.txt", plot_dir))[1]

assoc_info = readLines(sprintf("%s/%s_%s.info", plot_dir, chrom, pos))
ref = assoc_info[1]
period = assoc_info[2]
coeff = as.numeric(assoc_info[3])
intercept = as.numeric(assoc_info[4])
pval = as.numeric(assoc_info[5])

data = read_csv(sprintf("%s/%s_%s.csv", plot_dir, chrom, pos))
data = data %>% drop_na()
res_col = sprintf('%s_residual', phenotype)

# for each category, if it has length < 100, move all its values to extreme
# otherwise move the 25 lowest and 25 highest points to extreme
mark_extreme = function(single_gt_data, ...) {
    nobs = nrow(single_gt_data)
    if (nobs < 200) {
	return(single_gt_data %>% add_column(extreme = TRUE))
    } else {
	return(single_gt_data %>% 
	       arrange(!!as.name(res_col)) %>%
	       mutate(extreme = if_else(row_number() <= 50 | row_number() >= nobs-49, TRUE, FALSE))
        )
    }
}
data = data %>%
	group_by(gt) %>%
	group_modify(mark_extreme) %>% 
	ungroup() %>% 
	mutate(factor_gt = as.factor(gt))

gt_counts = data %>% 
    group_by(gt) %>% 
    summarize(n=n()) %>%
    mutate(xlab = sprintf("%0.1f\nn = %i", gt, n))

#data %>% slice_min(!!as.name(res_col), n=100) %>% print(n=100)

#png(sprintf("%s/%s_%s.png", plot_dir, chrom, pos))

p1 = ggplot() + geom_dotplot(
    data=data %>% filter(extreme == TRUE),
    mapping=aes_string(x='factor_gt', y=res_col),
    binaxis="y",
    stackdir="center",
    dotsize=0.2,
    binwidth=(max(data[[res_col]]) - min(data[[res_col]]))/100
) + geom_violin(
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('factor_gt', res_col),
    scale="count",
    na.rm=TRUE
) + geom_abline(
    slope=coeff,
    intercept=intercept
) + ylab("") + xlab("Repeat unit diff from ref") +
scale_x_discrete(labels=gt_counts$xlab) +
ggtitle("With 50 points at extremes (per gt)")

p2 = ggplot() + geom_violin(
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('factor_gt', res_col),
    scale="count",
    na.rm=TRUE
) + geom_abline(
    slope=coeff,
    intercept=intercept
) + ylab("") + xlab("Repeat unit diff from ref") +
scale_x_discrete(labels=gt_counts$xlab) +
ggtitle("No extrema") + 
annotation_compass(
    sprintf("assoc coeff: %0.2e\nintercept: %0.2e\np-val: %0.2e", coeff, intercept, pval),
    'NE'
)
         	

plot = arrangeGrob(
    p1,
    p2,
    nrow=1,
    top=sprintf("Linear association genotype ~ residual %s at %s:%s\n(residual after accounting for covariates)", phenotype, chrom, pos),
    left=sprintf("residual %s (%s)\n", phenotype, unit),
    bottom=sprintf(
        "Reference allele length: %i, repeat unit size: %s\nreference allele: %s",
	nchar(ref), period, ref
    )
)

ggsave(sprintf("%s/%s_%s.png", plot_dir, chrom, pos), plot=plot, height=7, width=14)

