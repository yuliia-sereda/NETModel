# NECModel 1.0.0

First public release, accompanying the manuscript on treatment sequencing in
small bowel neuroendocrine tumours.

* Continuous-time microsimulation of 8 sequences of PRRT, everolimus and
  cabozantinib (`Population`, `sim_data()`, `sim_parameter_sets()`).
* Calibration of post-progression mortality and per-dose MDS/AML incidence.
* Choice functions for decisions under ambiguity and Hasse diagrams.
* Analysis scripts, OSCAR (SLURM) pipelines and simulation outputs for both
  parameter sets (`COMPETE`, `NETTER&RADIANT`).
* Shiny clinician app (`app/`), deployed at
  <https://ys-shiny.shinyapps.io/NEC-app/>.
