# Treatment sequencing in small bowel NETs: clinician app

An interactive Shiny app showing simulated outcomes of 8 treatment strategies
for patients with small bowel neuroendocrine tumours who progress on
somatostatin analogues: sequences of peptide receptor radionuclide therapy
(PRRT), everolimus and cabozantinib.

Users can change two assumptions that trials cannot answer — how much less
effective each drug is when given later in the sequence, and how high the risk
of MDS/AML after PRRT is — and see, for every strategy:

- mean and median overall and progression-free survival, survival at 1, 3 and
  5 years, and MDS/AML cumulative incidence at 1, 3 and 5 years, each with a 95%
  range reflecting uncertainty in the trial data;
- the trade-off between survival and MDS/AML;
- differences from a chosen reference strategy;
- which strategies are optimal over a range of assumptions (not uniformly worse,
  best worst case, best best case).

The app only displays precomputed simulation results (`data/`); it does not run
the model. This is a research tool and does not replace clinical judgement for
individual patients.

## Running the app

```r
install.packages(c("shiny", "bslib", "data.table", "ggplot2", "ggrepel", "reactable", "MASS"))
shiny::runApp("app")   # from the repository root
```

## Structure

| Path | Content |
|---|---|
| `app.R` | User interface and server |
| `R/helpers.R` | Strategy definitions, formatting, plots, and choice functions (sourced automatically by Shiny) |
| `data/` | Simulation results prepared for the app |
| `www/model_diagram.png` | Diagram of the model states and transitions (Methods page) |
| `www/prrt_tradeoff.webp` | Illustration of the clinical trade-off: start with PRRT or delay it (About page) |

### Data

| File | Content |
|---|---|
| `design.rds` | The simulated combinations of assumptions: efficacy loss in the 2nd line and additional loss in the 3rd line (0-50% each), MDS/AML risk multiplier (1-3). The app uses those in which the 3rd-line loss is at least the 2nd-line loss |
| `outcomes.rds` | Mean and 95% range of each outcome, by combination and strategy |
| `vs_reference.rds` | Differences between strategies, with 95% ranges and the probability of being better |
| `clouds.rds` | Mean OS, mean PFS and 5-year MDS/AML incidence in each of the 1000 parameter sets, for the uncertainty regions |
| `outcome_info.rds` | Outcome labels and units |
| `inputs.rds` | Model inputs and sources |

The data are produced from the NECModel simulations (1000 parameter sets of
100,000 simulated patients at each combination) by
`analysis/app/04_prepare_app_data.R` in the repository root.

## Deployment

```r
rsconnect::deployApp("app")   # from the repository root
```

## Authors

Thomas A. Trikalinos and Yuliia Sereda.
