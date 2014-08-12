#
#
#
context("MINC I/O reading")

# silence the output of mincLm, in order to make the test output information is more clear to read
sink("/dev/null"); testfile <- mincGetVolume("/projects/moush/matthijs/2014-04-RMINC-testbed/brain_cut_out.mnc"); sink();
mincextract_output_voxel_0_0_0 <- as.numeric(system("mincextract  -start 0,0,0 -count 1,1,1 /projects/moush/matthijs/2014-04-RMINC-testbed/brain_cut_out.mnc", intern=TRUE))

sink("/dev/null"); mask_10_10_10 <- mincGetVolume("/projects/moush/matthijs/2014-04-RMINC-testbed/mask_at_voxel_10_10_10.mnc"); sink();
sink("/dev/null"); mask_10_31_33 <- mincGetVolume("/projects/moush/matthijs/2014-04-RMINC-testbed/mask_at_voxel_10_31_33.mnc"); sink();
sink("/dev/null"); mask_37_40_28 <- mincGetVolume("/projects/moush/matthijs/2014-04-RMINC-testbed/mask_at_voxel_37_40_28.mnc"); sink();
mincextract_output_voxel_10_10_10 <- as.numeric(system("mincextract  -start 10,10,10 -count 1,1,1 /projects/moush/matthijs/2014-04-RMINC-testbed/brain_cut_out.mnc", intern=TRUE))
mincextract_output_voxel_10_31_33 <- as.numeric(system("mincextract  -start 10,31,33 -count 1,1,1 /projects/moush/matthijs/2014-04-RMINC-testbed/brain_cut_out.mnc", intern=TRUE))
mincextract_output_voxel_37_40_28 <- as.numeric(system("mincextract  -start 37,40,28 -count 1,1,1 /projects/moush/matthijs/2014-04-RMINC-testbed/brain_cut_out.mnc", intern=TRUE))

test_that("mincGetVolume extracts the same value as mincextract", {
    expect_that(testfile[1], equals(mincextract_output_voxel_0_0_0))
    expect_that(testfile[mask_10_10_10 > 0.5], equals(mincextract_output_voxel_10_10_10))
    expect_that(testfile[mask_10_31_33 > 0.5], equals(mincextract_output_voxel_10_31_33))
    expect_that(testfile[mask_37_40_28 > 0.5], equals(mincextract_output_voxel_37_40_28))
})

#
#
#
context("MINC I/O writing")
sink("/dev/null"); mincWriteVolume(testfile, "/tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc", clobber=TRUE); sink();
reread_mincextract_output_voxel_0_0_0 <- as.numeric(system("mincextract  -start 0,0,0 -count 1,1,1 /tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc", intern=TRUE))
reread_mincextract_output_voxel_10_10_10 <- as.numeric(system("mincextract  -start 10,10,10 -count 1,1,1 /tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc", intern=TRUE))
reread_mincextract_output_voxel_10_31_33 <- as.numeric(system("mincextract  -start 10,31,33 -count 1,1,1 /tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc", intern=TRUE))
reread_mincextract_output_voxel_37_40_28 <- as.numeric(system("mincextract  -start 37,40,28 -count 1,1,1 /tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc", intern=TRUE))

test_that("mincGetVolume extracts the same value as mincextract", {
    expect_that(object = reread_mincextract_output_voxel_0_0_0, equals(mincextract_output_voxel_0_0_0, tolerance = 0.0001))
    expect_that(reread_mincextract_output_voxel_10_10_10, equals(mincextract_output_voxel_10_10_10, tolerance = 0.0001))
    expect_that(reread_mincextract_output_voxel_10_31_33, equals(mincextract_output_voxel_10_31_33, tolerance = 0.0001))
    expect_that(reread_mincextract_output_voxel_37_40_28, equals(mincextract_output_voxel_37_40_28, tolerance = 0.0001))
})

#
#
#
context("MINC I/O read from write")
sink("/dev/null"); reread_testfile <- mincGetVolume("/tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc"); sink();

test_that("mincGetVolume extracts the same value as mincextract", {
    expect_that(reread_testfile[1], equals(reread_mincextract_output_voxel_0_0_0))
    expect_that(reread_testfile[mask_10_10_10 > 0.5], equals(reread_mincextract_output_voxel_10_10_10))
    expect_that(reread_testfile[mask_10_31_33 > 0.5], equals(reread_mincextract_output_voxel_10_31_33))
    expect_that(reread_testfile[mask_37_40_28 > 0.5], equals(reread_mincextract_output_voxel_37_40_28))
})

# remove temp file
system("rm -f /tmp/write_out_of_RMINC_test_bed_MINC_IO.mnc")

#
#
#
context("MINC I/O determine volume of segmentation - anatGetAll using jacobians")

filenames <- read.csv("/projects/moush/matthijs/2014-04-RMINC-testbed/filenames.csv")
volumes <- anatGetAll(filenames=filenames$absolute_jacobian, atlas="/projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc", method="jacobians")

# calculate the volume of a structure using minc tools:
# 1) get the exponent of the absolute jacobian file (these are log values)
# 2) sum over a certain area
# 3) multiply that by the voxel size

system("mincmath -clobber -quiet -exp /projects/moush/matthijs/2014-04-RMINC-testbed/absolute_jacobian_file_1.mnc /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_1_exponent.mnc")
system("mincmath -clobber -quiet -exp /projects/moush/matthijs/2014-04-RMINC-testbed/absolute_jacobian_file_2.mnc /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_2_exponent.mnc") 
system("mincmath -clobber -quiet -exp /projects/moush/matthijs/2014-04-RMINC-testbed/absolute_jacobian_file_3.mnc /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_3_exponent.mnc") 

# left (180) and right (181) parieto-temporal lobe
left_parieto_sum_file_1 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 180 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_1_exponent.mnc", intern=TRUE))
left_parieto_sum_file_2 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 180 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_2_exponent.mnc", intern=TRUE))
left_parieto_sum_file_3 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 180 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_3_exponent.mnc", intern=TRUE))
right_parieto_sum_file_1 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 181 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_1_exponent.mnc", intern=TRUE))
right_parieto_sum_file_2 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 181 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_2_exponent.mnc", intern=TRUE))
right_parieto_sum_file_3 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 181 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_3_exponent.mnc", intern=TRUE))

# pons (187) has no separate label for left and right 
pons_sum_file_1 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 187 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_1_exponent.mnc", intern=TRUE))
pons_sum_file_2 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 187 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_2_exponent.mnc", intern=TRUE))
pons_sum_file_3 <- as.numeric(system("mincstats -quiet -sum -mask /projects/moush/matthijs/2014-04-RMINC-testbed/test_segmentation.mnc  -mask_bin 187 /tmp/rminctestdata/RMINC_test_bed_MINC_IO_absolute_jacobian_file_3_exponent.mnc", intern=TRUE))

voxel_volume <- 0.056*0.056*0.056

left_parieto_volume_file_1 <- left_parieto_sum_file_1 * voxel_volume
left_parieto_volume_file_2 <- left_parieto_sum_file_2 * voxel_volume
left_parieto_volume_file_3 <- left_parieto_sum_file_3 * voxel_volume
right_parieto_volume_file_1 <- right_parieto_sum_file_1 * voxel_volume
right_parieto_volume_file_2 <- right_parieto_sum_file_2 * voxel_volume
right_parieto_volume_file_3 <- right_parieto_sum_file_3 * voxel_volume
pons_volume_file_1 <- pons_sum_file_1 * voxel_volume
pons_volume_file_2 <- pons_sum_file_2 * voxel_volume
pons_volume_file_3 <- pons_sum_file_3 * voxel_volume


test_that("anatGetAll extracts the same value as a combination of mincmath and mincstats", {
    expect_that(volumes[, "left cerebral cortex: parieto-temporal lobe"][1], equals(left_parieto_volume_file_1, tolerance = 0.0001))
    expect_that(volumes[, "left cerebral cortex: parieto-temporal lobe"][2], equals(left_parieto_volume_file_2, tolerance = 0.0001))
    expect_that(volumes[, "left cerebral cortex: parieto-temporal lobe"][3], equals(left_parieto_volume_file_3, tolerance = 0.0001))

    expect_that(volumes[, "right cerebral cortex: parieto-temporal lobe"][1], equals(right_parieto_volume_file_1, tolerance = 0.0001))
    expect_that(volumes[, "right cerebral cortex: parieto-temporal lobe"][2], equals(right_parieto_volume_file_2, tolerance = 0.0001))
    expect_that(volumes[, "right cerebral cortex: parieto-temporal lobe"][3], equals(right_parieto_volume_file_3, tolerance = 0.0001))

    expect_that(volumes[, "pons"][1], equals(pons_volume_file_1, tolerance = 0.0001))
    expect_that(volumes[, "pons"][2], equals(pons_volume_file_2, tolerance = 0.0001))
    expect_that(volumes[, "pons"][3], equals(pons_volume_file_3, tolerance = 0.0001))
})


#
#
#
context("MINC I/O determine volume of segmentation - anatCombineStructures using jacobians")

volumes_combined <- anatCombineStructures(vols=volumes, method="jacobians")

test_that("anatGetAll extracts the same value as a combination of mincmath and mincstats", {
    expect_that(volumes_combined[, "cerebral cortex: parieto-temporal lobe"][1], equals(left_parieto_volume_file_1 + right_parieto_volume_file_1, tolerance = 0.0001))
    expect_that(volumes_combined[, "cerebral cortex: parieto-temporal lobe"][2], equals(left_parieto_volume_file_2 + right_parieto_volume_file_2, tolerance = 0.0001))
    expect_that(volumes_combined[, "cerebral cortex: parieto-temporal lobe"][3], equals(left_parieto_volume_file_3 + right_parieto_volume_file_3, tolerance = 0.0001))

    expect_that(volumes_combined[, "pons"][1], equals(pons_volume_file_1, tolerance = 0.0001))
    expect_that(volumes_combined[, "pons"][2], equals(pons_volume_file_2, tolerance = 0.0001))
    expect_that(volumes_combined[, "pons"][3], equals(pons_volume_file_3, tolerance = 0.0001))
})


