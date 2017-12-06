library(RMINC)
library(dplyr)

civet_summaries <- 
  read.csv("/hpf/largeprojects/MICe/chammill/POND/volume_summaries/combined_neuroanatomy20170405.csv"
           , stringsAsFactors = FALSE)


civet_files_PND <-
  civet.getAllFilenames(
    filter(civet_summaries, src == "POND") %>% mutate(scan_np = gsub("-", "_", sub("MR160-", "", scan)))
    , idvar = "scan_np", prefix = "MR160_", basedir = "/hpf/largeprojects/MICe/chammill/POND/civetOutputs/pond_2-1_20170227/"
    , civetVersion = "2.1.0"
    , cnf = yaml::yaml.load_file("/hpf/largeprojects/MICe/chammill/POND/civetOutputs/pond_2-1_20170227/088_0002_01_002/CBRAIN_Colosse-384672-1.params.yml")
  )

tvt <- 
  civet_files_PND %>%
  filter(Dx %in% c("CTRL", "ASD")) %>%
  group_by(Dx) %>%
  slice(1:10) %>%
  ungroup

vlm <-  
  vertexLm(tvt$nativeRMS_RSLtlink_left ~ Dx, data = tvt)

left_surf <- read_obj("/hpf/largeprojects/MICe/chammill/packageSources/CIVET_2.0/CIVET_2.0_icbm_avg_mid_sym_mc_left.obj")

plot(left_surf, colour_map = vlm[,"tvalue-DxCTRL"])
