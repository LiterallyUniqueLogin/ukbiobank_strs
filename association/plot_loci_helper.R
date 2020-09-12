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
data = data %>% mutate(first = ifelse(row_number() <= 10, TRUE, FALSE))
res_col = sprintf('%s_residual', phenotype)

# for each category, if it has length < 100, move all its values to extreme
# otherwise move the 25 lowest and 25 highest points to extreme
mark_extreme = function(single_gt_data, ...) {
    nobs = nrow(single_gt_data)
    if (nobs < 200) {
	return(single_gt_data %>% add_column(extreme = TRUE) %>% add_column(filtered = TRUE))
    } else {
	return(single_gt_data %>% 
	       arrange(!!as.name(res_col)) %>%
	       mutate(extreme = if_else(row_number() <= 50 | row_number() >= nobs-49, TRUE, FALSE)) %>%
	       add_column(filtered = FALSE)
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
    mutate(xlab = sprintf("%0.1f\n%sn=%i", gt, ifelse(row_number() %% 2 == 0, '\n', ''), n))

p1 = ggplot() + geom_dotplot(
    data=data %>% filter(extreme == TRUE),
    mapping=aes_string('gt', res_col, group='factor_gt'),
    binaxis="y",
    stackdir="center",
    dotsize=0.2,
    binwidth=(max(data[[res_col]]) - min(data[[res_col]]))/100
) + geom_violin(
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('gt', res_col, group='factor_gt'),
    scale="count",
    na.rm=TRUE
) + geom_abline(
    slope=coeff,
    intercept=intercept
) + ylab("") + xlab("Avg diff from ref (repeats)") +
scale_x_continuous(breaks=gt_counts$gt, labels=gt_counts$xlab) +
theme(axis.text.x = element_text(size=9)) +
ggtitle("With 50 points at extremes (per gt)") +
labs(caption=" ", alpha=0) #space this graph the same as the one below

p2 = ggplot( 
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('gt', res_col, group='factor_gt')
) + geom_violin(
    scale="count",
    na.rm=TRUE
) + stat_summary(
    data=data %>% filter(filtered == FALSE),
    aes(shape="mean"),
    fun.data='mean_cl_normal',
    color='red',
    na.rm=TRUE,
    size=0.3
) + geom_abline(
    slope=coeff,
    intercept=intercept,
    show.legend=TRUE
) + ylab("") + xlab("Avg diff from ref (repeats)") +
scale_shape_manual("", values=c("mean"=19)) +
scale_x_continuous(breaks=gt_counts$gt, labels=gt_counts$xlab) +
theme(axis.text.x = element_text(size=9)) +
ggtitle("Extrema not plotted") +
labs(caption=sprintf("(%s values are residuals after regressing out covariates)", phenotype))
#+ geom_dotplot(
#    data=data %>% filter(first == TRUE),
#    mapping=aes_string('gt', res_col, group='factor_gt'),
#    binaxis="y",
#    stackdir="center",
#    dotsize=0.5,
#    color='red',
#    binwidth=(max(data[[res_col]]) - min(data[[res_col]]))/100
#)

plot = arrangeGrob(
    p1,
    p2,
    nrow=1,
    top=sprintf("STR length ~ %s\nLinear association at %s:%s\nP-value: %0.2e    Slope: %0.2e    Intercept %0.2f", phenotype, chrom, pos, pval, coeff, intercept),
    left=sprintf("%s (%s)\n", phenotype, unit),
    bottom=sprintf(
        "Reference allele length: %i, repeat unit size: %s\nreference allele: %s",
	nchar(ref), period, ref
    )
)

ggsave(sprintf("%s/%s_%s.png", plot_dir, chrom, pos), plot=plot, height=7, width=14)

