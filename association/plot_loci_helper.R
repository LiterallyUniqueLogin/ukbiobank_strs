#!/usr/bin/env Rscript
library(tidyverse)
library(gridExtra)

args = commandArgs(trailingOnly=TRUE)
plot_dir = args[1]
chrom = args[2]
pos = args[3]
phenotype = args[4]
transform_step = args[5]
unit = args[6]
period = args[7]
ref = args[8]

if (transform_step == 'original') {
    out_fname = sprintf("%s/%s_%s", plot_dir, chrom, pos)
} else {
    out_fname = sprintf("%s/%s_%s_%s", plot_dir, chrom, pos, transform_step)
}

data = read_csv(
    sprintf("%s.csv", out_fname),
    col_names=c('phenotype', 'gt')
)
# data = data %>% drop_na()
data = data %>% mutate(first = ifelse(row_number() <= 10, TRUE, FALSE))

# for each category, if it has length < 100, move all its values to extreme
# otherwise move the 25 lowest and 25 highest points to extreme
mark_extreme = function(single_gt_data, ...) {
    nobs = nrow(single_gt_data)
    if (nobs < 200) {
	return(single_gt_data %>% add_column(extreme = TRUE) %>% add_column(filtered = TRUE))
    } else {
	return(single_gt_data %>% 
	       arrange(!!as.name('phenotype')) %>%
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
    mapping=aes_string('gt', 'phenotype', group='factor_gt'),
    binaxis="y",
    stackdir="center",
    dotsize=0.2,
    binwidth=(max(data[['phenotype']]) - min(data[['phenotype']]))/100
) + geom_violin(
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('gt', 'phenotype', group='factor_gt'),
    scale="count",
    na.rm=TRUE
) + ylab("") + xlab("Avg diff from ref (repeats)") +
  scale_x_continuous(breaks=gt_counts$gt, labels=gt_counts$xlab) +
  theme(axis.text.x = element_text(size=9)) +
  ggtitle("With 50 points at extremes (per gt)") +
  labs(caption=" ", alpha=0) #space this graph the same as the one below

p2 = ggplot( 
    data=data %>% filter(extreme == FALSE),
    mapping=aes_string('gt', 'phenotype', group='factor_gt')
) + geom_violin(
    scale="count",
    na.rm=TRUE
) + stat_summary(
    data=data %>% filter(filtered == FALSE),
    aes(color='mean', shape="mean"),
    fun='mean',
    na.rm=TRUE,
    size=0.3
) + stat_summary(
    data=data %>% filter(filtered == FALSE),
    aes(group="factor_gt",
	color='mean'),
    fun='mean',
    na.rm=TRUE,
    geom="line",
) + ylab("") + xlab("Avg diff from ref (repeats)") +
  scale_shape_manual("", values=c("mean"=19)) + 
  scale_color_manual("", values=c("mean"="red", "best fit linear model"="black")) + 
  guides(shape="none", color=guide_legend("")) +
  scale_x_continuous(breaks=gt_counts$gt, labels=gt_counts$xlab) +
  theme(axis.text.x = element_text(size=9)) +
  ggtitle("Extrema not plotted")

plot = arrangeGrob(
    p1,
    p2,
    nrow=1,
    top=sprintf("STR length ~ %s\nLinear association at %s:%s", phenotype, chrom, pos),
    left=sprintf("%s (%s)\n", phenotype, unit),
    bottom=sprintf(
        "Reference allele length: %i, repeat unit size: %s\nreference allele: %s",
	nchar(ref), period, ref
    )
)

ggsave(sprintf("%s.png", out_fname), plot=plot, height=7, width=14)

