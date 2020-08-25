#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
plot_dir = args[1]
chrom = args[2]
pos = args[3]
phenotype = args[4]

unit = readLines(sprintf("%s/info.txt", plot_dir))[1]

assoc_info = readLines(sprintf("%s/%s_%s.info", plot_dir, chrom, pos))
ref = assoc_info[1]
period = assoc_info[2]
coeff = assoc_info[3]

df = read.csv(sprintf("%s/%s_%s.csv", plot_dir, chrom, pos))
df$factor_gt = as.factor(df$gt)

png(sprintf("%s/%s_%s.png", plot_dir, chrom, pos))

plot = ggplot(df, aes_string('factor_gt', sprintf('%s_residual', phenotype)))
plot + geom_violin(scale="count")

ignore = dev.off()


#
#   fig.suptitle(
#        f"Linear associations at individual loci against trait {args.phenotype} "
#        "(covariate associations not plotted)"
#    )
#    for plot_idx, locus in enumerate(args.locus):
#        ax = axs[plot_idx % 2, plot_idx//2]
#        chrom, pos = locus.split(":")
#
#        with open(f'{assoc_dir}/results/{args.phenotype}.txt') as results:
#            coeff = None
#            for line in results:
#                line_chrom, line_pos = line.split()[:2]
#                if chrom == line_chrom and pos == line_pos:
#                    coeff = line.split()[5]
#                    coeff = float(coeff)
#                    break
#            if coeff is None:
#                raise ValueError(f"Couldn't find {locus} in results file"
#                                 f" {results.name}!")
#                
#
#        vcf = cyvcf2.VCF(
#            f'{ukb}/str_imputed/runs/{args.imputation_run_name}/vcfs/annotated_strs/chr{chrom}.vcf.gz'
#        )
#        record = next(iter(vcf(locus)))
#        ref = record.REF
#        period = record.INFO['PERIOD']
#        period = int(period)
#        bpdiffs = record.INFO['BPDIFFS']
#        bpdiffs = np.array(bpdiffs)
#        bpdiffs_repeat_units = bpdiffs/period
#
#        idx_gts = record.genotype.array()[:, :-1]
#        len_gts = np.zeros(idx_gts.shape)
#        for idx, diff in enumerate(bpdiffs_repeat_units):
#            len_gts[idx_gts == idx] = diff
#        avg_len_gts = np.sum(len_gts, axis=1)/2
#        '''
#        command = f"""
#        source ~/.bashrc ;
#        conda activate ukb ;
#        bcftools query -f "REF %BPDIFFS %PERIOD" \\
#            {ukb}/snpstr/info_field/chr{chrom}.vcf.gz
#        """
#        out = sp.run(command, shell=True, capture_output=True, text=True)
#        if out.stderr:
#            print("Command failed with stderr:\n", out.stderr)
#            sys.exit(1)
#        ref, bpdiffs, period = out.stdout.split()
#        bpdiffs = bpdiffs.split(',')
#        '''
#
#        ax.set_title(f"Locus {locus}")
#        ax.set_ylabel(unit)
#        ax.set_xlabel("Diff from ref in repeat units")
#        ax.annotate(
#            f"Repeat period: {period}\tRef allele: {ref}\tLen: {len(ref)}",
#            (0, -20)
#        )
#
#        phenotypes_by_genotypes = []
#        axis_vals = []
#        unique_gts = np.unique(avg_len_gts)
#        for gt in unique_gts:
#            samples = np.array(vcf.samples)[avg_len_gts == gt]
#            samples = np.char.partition(samples, '_')[:, 0]
#            samples = samples.astype(int)
#            vals = data[np.isin(data[:, 0], samples), 1]
#            if len(vals) == 0:
#                continue
#
#            axis_vals.append(gt)
#            phenotypes_by_genotypes.append(
#                data[np.isin(data[:, 0], samples), 1]
#            )
#        axis_vals = np.array(axis_vals)
#        print(axis_vals)
#        ax.violinplot(
#            phenotypes_by_genotypes,
#            positions=axis_vals,
#        )
#        ax.plot(axis_vals, axis_vals*coeff,
#                label=f"Best fit line\n(slope == {coeff:.2e} {unit}/repeat)")
#        ax.legend()
#    outloc = f'{assoc_dir}/loci_plots_{"_".join(args.locus)}.png'
#    outloc = "_".join(outloc.split(':'))
#    plt.savefig(outloc)
#
#
#if __name__ == "__main__":
#    main()
