version 1.0

import "plot_one_wgs_assoc.wdl"

workflow plot_paper_wgs_assocs {

  call plot_one_wgs_assoc.plot_one_wgs_assoc as fig_4b { input :
    phenotype_name = "platelet_count",
    chrom = 11,
    pos = 119206290
  }

#  call plot_one_wgs_assoc.plot_one_wgs_assoc as fig_4f { input :
#    phenotype_name = "mean_platelet_count",
#    chrom = 11,
#    pos = 119206290
#  }

  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_15ba { input :
    phenotype_name = "total_bilirubin",
    chrom = 3,
    pos = 171009913 
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_16a { input :
    phenotype_name = "platelet_crit",
    chrom = 11,
    pos = 119206290
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_16b { input :
    phenotype_name = "platelet_count",
    chrom = 11,
    pos = 119206290,
    residual = true
  }

  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_17a { input :
    phenotype_name = "mean_sphered_cell_volume",
    chrom = 11,
    pos = 119206290
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_19b { input :
    phenotype_name = "eosinophil_percent",
    chrom = 2,
    pos = 111120967 
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_20a { input :
    phenotype_name = "mean_platelet_volume",
    chrom = 17,
    pos = 29514998
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_21a { input :
    phenotype_name = "haematocrit",
    chrom = 14,
    pos = 64247333
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_22a { input :
    phenotype_name = "red_blood_cell_distribution_width",
    chrom = 17,
    pos = 32142452
  }
  
  call plot_one_wgs_assoc.plot_one_wgs_assoc as supp_fig_23a { input :
    phenotype_name = "mean_platelet_volume",
    chrom = 2,
    pos = 105893985
  }

  output {
    File supp_fig_23a_png = supp_fig_23a.png
    File supp_fig_23a_svg = supp_fig_23a.svg
    File supp_fig_22a_png = supp_fig_22a.png
    File supp_fig_22a_svg = supp_fig_22a.svg
    File supp_fig_21a_png = supp_fig_21a.png
    File supp_fig_21a_svg = supp_fig_21a.svg
    File supp_fig_20a_png = supp_fig_20a.png
    File supp_fig_20a_svg = supp_fig_20a.svg
    File supp_fig_19b_png = supp_fig_19b.png
    File supp_fig_19b_svg = supp_fig_19b.svg
    File supp_fig_17a_png = supp_fig_17a.png
    File supp_fig_17a_svg = supp_fig_17a.svg
    File supp_fig_16a_png = supp_fig_16a.png
    File supp_fig_16a_svg = supp_fig_16a.svg
    File supp_fig_15ba_png = supp_fig_15ba.png
    File supp_fig_15ba_svg = supp_fig_15ba.svg
    File fig_4b_png = fig_4b.png
    File fig_4b_svg = fig_4b.svg
  }
}
