# Design grid for the clinician app
#
# Each later position can lose 0-50% of a treatment's effectiveness relative to
# the previous position, i.e. its progression rate is multiplied by (1 + a):
#   step2_attenuation a2: position 2 vs position 1, in {0, 10, ..., 50%}
#   step3_attenuation a3: position 3 vs position 2, in {0, 10, ..., 50%}
#   aml_mds_multiplier:   multiplier of the MDS/AML hazard, in {1, 1.5, ..., 3}
# Attenuation at position 3 is therefore always at least that at position 2.
# The model takes cumulative factors relative to position 1:
#   discount1 = 1 + a2, discount2 = (1 + a2) * (1 + a3)
# 6 x 6 x 5 = 180 design points.
#
# Run from the repository root. Output: design/app_design_grid.rds

library(data.table)

attenuation <- seq(0, 0.5, by = 0.1)
design <- CJ(
  step2_attenuation = attenuation,
  step3_attenuation = attenuation,
  aml_mds_multiplier = seq(1, 3, by = 0.5)
)
design[, `:=`(
  design_id = .I,
  discount1 = 1 + step2_attenuation,
  discount2 = (1 + step2_attenuation) * (1 + step3_attenuation)
)]
setcolorder(design, "design_id")

saveRDS(design, file.path("design", "app_design_grid.rds"))
