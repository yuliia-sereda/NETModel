# maxPro design of the uncertainty set
#
# Three factors:
#   step2 = attenuation of the treatment in position 2 relative to position 1
#   step3 = attenuation of the treatment in position 3 relative to position 2
#   aml_mds_multiplier = multiplier of the MDS/AML hazard
# Each subsequent position can lose up to 50% of effectiveness, i.e. its
# progression rate is up to 1.5 times that of the previous position:
#   1 <= step2, step3 <= 1.5 and 1 <= aml_mds_multiplier <= 3.
# The model takes cumulative factors relative to position 1:
#   discount1 = step2 (1 to 1.5), discount2 = step2 * step3 (1 to 2.25).
# Each factor is padded by 20% of its range on both sides, so the GAM emulator
# is well supported at the edges of the uncertainty set.
#
# Run from the repository root. Output: design/maxPro_design_3000.rds

library(data.table)
library(MaxPro)

set.seed(20260926)
N <- 3000

# Design in the unit cube
Z0 <- MaxProLHD(n = N, p = 3)
Z <- MaxPro(InitialDesign = Z0$Design)

design <- as.data.table(Z$Design)
colnames(design) <- c("unit_step2", "unit_step3", "unit_aml_mds_multiplier")

# Scale a unit-interval value to [lower, upper] padded by `pad` of the range
scale_padded <- function(u, lower, upper, pad = 0.2) {
  width <- upper - lower
  (lower - pad * width) + u * (1 + 2 * pad) * width
}

design[, step2 := scale_padded(unit_step2, 1, 1.5)]                             # 0.9 to 1.6
design[, step3 := scale_padded(unit_step3, 1, 1.5)]                             # 0.9 to 1.6
design[, aml_mds_multiplier := scale_padded(unit_aml_mds_multiplier, 1, 3)]     # 0.6 to 3.4
design[, `:=`(discount1 = step2, discount2 = step2 * step3)]

saveRDS(design, file.path("design", "maxPro_design_3000.rds"))
