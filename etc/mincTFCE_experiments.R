#/hpf/largeprojects/MICe/matthijs/2016-06-Boughner-embryos/2017-03-TFCE

library(RMINC)

gfE14.0 <- read.csv("/hpf/largeprojects/MICe/matthijs/2016-06-Boughner-embryos/Boughner-to-E14.0_analysis/filename_mapping.csv")
gfE14.0_no_outlier <- gfE14.0[-c(8),]
gfE14.0_no_outlier$abs_jac_fwhm0.2 <- sub(".mnc", "_fwhm0.2.mnc", gfE14.0_no_outlier$abs_jac_no_blur)
gfE14.0_no_outlier$rel_jac_fwhm0.2 <- sub("_with_additional_inverted_absolute", "_inverted_pure_nlin_relative", gfE14.0_no_outlier$abs_jac_fwhm0.2)

vs <- mincLm(rel_jac_fwhm0.2 ~ genotype, gfE14.0_no_outlier
           , mask="/hpf/largeprojects/MICe/matthijs/2016-06-Boughner-embryos/Boughner-to-E14.0/minctracc_nlin/minctracc-nlin-6_rough_mask_fix_zyx.mnc")

tfce <- mincTFCE(vs, R = 500, parallel = c("", 500))
