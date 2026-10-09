# Treatment sequencing in small bowel NETs: clinician app
#
# Shows simulated outcomes of the NECModel microsimulation for 8 treatment
# strategies. The app only reads precomputed simulation outputs from data/; it
# never runs the model.

library(shiny)
library(bslib)
library(data.table)
library(ggplot2)
library(reactable)

source("R/helpers.R")

# ------------------------------------------------------------------------------
# Data
# ------------------------------------------------------------------------------
design <- readRDS("data/design.rds")
info <- readRDS("data/outcome_info.rds")
outcomes <- readRDS("data/outcomes.rds")
vs_reference <- readRDS("data/vs_reference.rds")
clouds <- readRDS("data/clouds.rds")
inputs <- readRDS("data/inputs.rds")

grouped_outcomes <- outcome_choices(info)

# Sources of the model inputs, with PubMed IDs
references <- data.table(
  citation = c(
    "Walter T, Jann H, Ansquer C, et al. [(177)Lu]Lu-edotreotide versus everolimus for gastroenteropancreatic neuroendocrine tumours (COMPETE): a phase 3, multicentre, randomised, open-label, superiority trial. Lancet 2026; 408(10551): 234-47.",
    "Chan JA, Geyer S, Zemla T, et al. Phase 3 Trial of Cabozantinib to Treat Advanced Neuroendocrine Tumors. N Engl J Med 2025; 392(7): 653-65.",
    "Strosberg JR, Caplin ME, Kunz PL, et al. (177)Lu-Dotatate plus long-acting octreotide versus high-dose long-acting octreotide in patients with midgut neuroendocrine tumours (NETTER-1): final overall survival and long-term safety results from an open-label, randomised, controlled, phase 3 trial. Lancet Oncol 2021; 22(12): 1752-63.",
    "Batista JPT, De La Iglesia M, Kim H, et al. Hematological neoplasms after PRRT: experience from long term follow-up of 173 patients in a tertiary center. Endocrine Oncology 2026; 6(1).",
    "Bernard E, Tuechler H, Greenberg PL, et al. Molecular International Prognostic Scoring System for Myelodysplastic Syndromes. NEJM Evid 2022; 1(7): EVIDoa2200008."
  ),
  pmid = c("42392118", "39282913", "34793718", "42694695", "38319256")
)

# Choices for the scenario dropdowns (values as simulated)
loss_choices <- function(from = 0) {
  values <- seq(from, 50, by = 10)
  setNames(values, paste0(values, "%"))
}
mds_choices <- setNames(seq(1, 3, by = 0.5), c("1× (trial risk)", paste0(seq(1.5, 3, by = 0.5), "×")))

# A range as two dropdowns, "from" and "to" (ids <id>_from and <id>_to)
range_input <- function(id, label, choices, from, to) {
  tagList(
    tags$label(label, class = "control-label"),
    layout_columns(
      col_widths = c(6, 6), gap = "0.5rem",
      selectInput(paste0(id, "_from"), "from", choices = choices, selected = from),
      selectInput(paste0(id, "_to"), "to", choices = choices, selected = to)
    )
  )
}

# Slightly smaller text in the result tables
options(reactable.theme = reactableTheme(style = list(fontSize = "0.88rem")))

# ------------------------------------------------------------------------------
# UI
# ------------------------------------------------------------------------------
scenario_controls <- conditionalPanel(
  "['outcomes', 'compare'].includes(input.nav)",
  h5("Scenario"),
  selectInput("loss2", "How much less effective is a drug when given second rather than first?",
              choices = loss_choices(), selected = 0),
  selectInput("loss3", "...and how much less again when given third (at least as much as when given second)?",
              choices = loss_choices(), selected = 0),
  selectInput("mds", "MDS/AML risk after PRRT, compared with trial data",
              choices = mds_choices, selected = 1),
  actionButton("preset_none", "No efficacy loss, trial MDS/AML risk", class = "btn-sm w-100 mb-2"),
  actionButton("preset_max", "Maximum efficacy loss, 3× MDS/AML risk", class = "btn-sm w-100"),
  conditionalPanel(
    "input.nav == 'compare'",
    hr(),
    selectInput("reference", "Compare all strategies with", choices = strategy_choices, selected = "PCE")
  )
)

optimal_controls <- conditionalPanel(
  "input.nav == 'optimal'",
  h5("Range of assumptions"),
  range_input("opt_loss2", "Efficacy loss when given second", loss_choices(), 0, 50),
  range_input("opt_loss3", "Additional loss when given third", loss_choices(), 0, 50),
  helpText("The additional loss when given third is at least the loss when given second: only such
           combinations are included."),
  range_input("opt_mds", "MDS/AML risk after PRRT", mds_choices, 1, 3),
  selectInput("opt_outcome", "Outcome", choices = grouped_outcomes, selected = "os_mean"),
  uiOutput("opt_margin_ui")
)

ui <- page_navbar(
  title = "Treatment sequencing in small bowel NETs",
  id = "nav",
  theme = bs_theme(version = 5, preset = "shiny"),
  # Pages grow with their content and scroll, rather than squeezing cards into the window
  fillable = FALSE,
  sidebar = sidebar(width = 320, scenario_controls, optimal_controls,
                    conditionalPanel("['about', 'methods'].includes(input.nav)",
                                     p("Choose a page above to explore the results."))),

  nav_panel(
    "About", value = "about",
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(7, 5)),
      card(
        card_header("Background"),
        p("Patients with small bowel neuroendocrine tumours progressing on somatostatin analogues can receive
          peptide receptor radionuclide therapy (PRRT; P), everolimus (E), or cabozantinib (C). The optimal
          treatment sequence is unknown. PRRT provides the longest progression-free survival but carries a risk of
          myelodysplastic syndrome and acute myeloid leukaemia (MDS/AML), which are life-threatening blood malignancies
          Delaying PRRT may reduce this risk, but its effectiveness may decline with later use, as may that of
          everolimus and cabozantinib. How much treatment effects decline, and how large the MDS/AML risk is,
          remain uncertain."),
        p("This creates a clinical trade-off: start with PRRT or delay it?"),
        p("We developed a discrete-event simulation model to compare eight treatment strategies: six using all
          three drugs and two excluding cabozantinib. This app allows users to explore outcomes under different
          assumptions about treatment effect attenuation at later lines and varying MDS/AML risk, and to identify 
          strategies that perform well across a wide range of assumptions.")
      ),
      card(
        tags$img(src = "prrt_tradeoff.webp",
                 alt = "Start with PRRT or delay it? Starting PRRT first gives more effective treatment upfront but a higher risk of treatment-related MDS/AML; delaying it lowers the immediate MDS/AML risk but may reduce the effectiveness of later treatment.",
                 style = "width: 100%; max-width: 520px; height: auto; display: block; margin: auto;")
      )
    ),
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(7, 5)),
      card(
        card_header("How to use this app"),
        tags$dl(
          tags$dt("Explore a scenario"),
          tags$dd("Choose one set of assumptions on the left: how much less effective each drug is when given
                  second or third, and the MDS/AML risk after PRRT. The page shows each strategy's survival and
                  MDS/AML outcomes under these assumptions, with 95% uncertainty bounds that reflect uncertainty in
                  the model parameters, and the trade-off between survival and MDS/AML. Results can be downloaded."),
          tags$dt("Pairwise comparisons"),
          tags$dd("For the same scenario, shows the difference between each strategy and a reference strategy of
                  your choice."),
          tags$dt("Robust strategies"),
          tags$dd("These assumptions cannot be estimated from available trial data. Instead of choosing one scenario, choose
                  the range of assumptions you consider plausible. The page identifies strategies that remain good
                  choices across the whole range: those that are not uniformly worse than another strategy, and
                  those with the best worst case and the best best case."),
          tags$dt("Methods"),
          tags$dd("Model specification, input parameters, and simulation details.")
        )
      ),
      card(
        card_header("Strategies compared"),
        reactableOutput("strategy_table")
      )
    )
  ),

  nav_panel(
    "Explore a scenario", value = "outcomes",
    p(class = "text-muted",
      "Results for the set of assumptions chosen on the left. The 95% uncertainty bounds reflect uncertainty in
      the model parameters."),
    radioButtons("outcome_group", NULL, inline = TRUE, choices = names(grouped_outcomes)),
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(7, 5)),
      card(
        card_header("Outcomes of each strategy (estimate and 95% uncertainty bounds)"),
        reactableOutput("summary_table"),
        card_footer(
          "Survival times in years; survival and MDS/AML in %. Best value in each column in bold.",
          div(class = "mt-2",
              downloadButton("download_scenario", "This scenario (CSV)", class = "btn-sm"),
              downloadButton("download_all", "All scenarios (CSV)", class = "btn-sm"))
        )
      ),
      card(
        full_screen = TRUE,
        card_header("One outcome at a time"),
        uiOutput("interval_outcome_ui"),
        plotOutput("interval_plot", height = "auto")
      )
    ),
    card(
      full_screen = TRUE,
      card_header("Trade-off between survival and MDS/AML"),
      radioButtons("tradeoff_y", NULL, inline = TRUE,
                   choices = c("Mean overall survival" = "os_mean", "Mean progression-free survival" = "pfs_mean")),
      plotOutput("tradeoff_plot", height = "auto"),
      card_footer(
        "Dots show the mean survival and 5-year MDS/AML cumulative incidence of each strategy in the selected
        scenario, averaged over 1,000 simulated populations of 100,000 patients. Outlines show the 95%
        highest-density regions of the two outcomes, reflecting uncertainty propagated through probabilistic
        sampling of model parameters. Better strategies lie towards the top left (longer survival, less MDS/AML)."
      )
    )
  ),

  nav_panel(
    "Pairwise comparisons", value = "compare",
    h5(textOutput("compare_header", inline = TRUE)),
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(5, 7)),
      card(
        card_header("Select outcome"),
        selectInput("compare_outcome", NULL, choices = grouped_outcomes, selected = "os_mean"),
        reactableOutput("difference_table")
      ),
      card(
        full_screen = TRUE,
        card_header("Survival and MDS/AML compared with the reference"),
        radioButtons("compare_y", NULL, inline = TRUE,
                     choices = c("Mean overall survival" = "os_mean", "Mean progression-free survival" = "pfs_mean")),
        plotOutput("difference_plot", height = "auto")
      )
    )
  ),

  nav_panel(
    "Robust strategies", value = "optimal",
    card(
      card_header("Choosing strategies when the assumptions are uncertain"),
      p("The assumptions on the left are not known. Instead of picking one value, choose the range you consider
        plausible. Each strategy's average outcome then lies somewhere within a range, and there are different
        ways to choose among strategies:"),
      tags$ul(
        tags$li(strong("Not uniformly worse (maximality): "),
                "a strategy is excluded only if another strategy is better in every combination of assumptions
                in the range. Usually several strategies remain."),
        tags$li(strong("Best worst case (pessimism): "),
                "the strategy whose outcome in its least favourable combination of assumptions is the best.
                Suits a cautious decision-maker."),
        tags$li(strong("Best best case (optimism): "),
                "the strategy whose outcome in its most favourable combination of assumptions is the best.")
      ),
      p("A difference smaller than the clinical indifference margin is treated as no difference.")
    ),
    layout_columns(
      col_widths = breakpoints(sm = 12, lg = c(6, 6)),
      card(
        card_header("Strategies robust across the selected range of assumptions"),
        reactableOutput("optimal_table")
      ),
      card(
        full_screen = TRUE,
        card_header("Range of each strategy's average outcome across the selected combinations"),
        plotOutput("range_plot", height = "auto")
      )
    )
  ),

  nav_panel(
    "Methods", value = "methods",
    card(
      card_header("How the model works"),
      p(strong("Simulated population. "),
        "We simulated a hypothetical population of patients eligible for all three treatments. The population was
        50% female, with ages following a normal distribution (mean 63 years, standard deviation 9 years) truncated
        to 18–84 years, based on the NETTER-1 trial. Outcomes were simulated over each patient's remaining
        lifetime, up to age 100."),
      p(strong("Model states and transitions. "),
        "All patients entered the simulation with stable disease on the first treatment in their assigned
        sequence. Patients remained on treatment until disease progression, defined using RECIST 1.1, and then
        moved to the next treatment in the sequence. After progression on the final treatment, patients entered
        salvage therapy. Everolimus and cabozantinib were administered daily, while PRRT consisted of four
        infusions given two months apart."),
      p(strong("Competing events. "),
        "The model simulated competing events: disease progression, treatment-related MDS/AML, and death. Each
        PRRT infusion permanently increased the hazard of MDS/AML. Patients who developed MDS/AML stopped their
        neuroendocrine tumour treatment sequence. Death could occur from disease progression after salvage
        therapy, from MDS/AML after its onset, or from other causes at any time. Background mortality depended on
        age, sex, and birth cohort."),
      p(strong("Uncertainty propagation. "),
        "Results are reported as means and 95% uncertainty intervals (2.5th to 97.5th percentiles) across 1,000
        parameter sets drawn from the input distributions, with 100,000 patients simulated per set to minimise
        Monte Carlo error. Each parameter set was simulated at every combination of the assumptions offered in
        the app (105 combinations).")
    ),
    card(
      card_header("Diagram of the model states and transitions"),
      layout_columns(
        col_widths = breakpoints(sm = 12, lg = c(7, 5)),
        tags$img(src = "model_diagram.png",
                 alt = "Diagram of the model: stable disease on treatments I, II and III, salvage treatment, MDS/AML, and death from progression, MDS/AML or other causes, with the transitions between them.",
                 style = "width: 100%; max-width: 720px; height: auto; display: block; margin: auto;"),
        p(class = "text-muted small",
          "λ denotes the rate of each transition. Progression rates (λ12, λ23, λ34) depend on
          the treatment received and change with the efficacy-loss assumptions; the MDS/AML rate (λ∘5)
          changes with the MDS/AML risk assumption. Death from other causes (λ∘8) depends on age, sex,
          and birth cohort.")
      )
    ),
    card(
      card_header("Model inputs"),
      reactableOutput("inputs_table")
    ),
    card(
      card_header("References"),
      tags$ol(
        lapply(seq_len(nrow(references)), function(i) {
          tags$li(
            class = "mb-2",
            references$citation[i], " ",
            tags$a(href = paste0("https://pubmed.ncbi.nlm.nih.gov/", references$pmid[i], "/"),
                   target = "_blank", rel = "noopener", paste0("PubMed ", references$pmid[i]))
          )
        })
      )
    )
  )
)

# ------------------------------------------------------------------------------
# Server
# ------------------------------------------------------------------------------
server <- function(input, output, session) {

  # Plot height from the plot's current width (keeps the shape on any screen),
  # within sensible limits
  height_from_width <- function(output_id, ratio, min = 320, max = 760) {
    function() {
      width <- session$clientData[[paste0("output_", output_id, "_width")]]
      if (is.null(width) || width == 0) return(min)
      max(min(width * ratio, max), min)
    }
  }

  # Scenario as numbers (dropdown values arrive as text)
  scenario <- reactive({
    req(input$loss2, input$loss3, input$mds)
    list(loss2 = as.integer(input$loss2), loss3 = as.integer(input$loss3), mds = as.numeric(input$mds))
  })
  
  # The additional loss when given third is at least the loss when given second:
  # offer only those values, and raise the current one if needed
  observeEvent(input$loss2, {
    loss2 <- as.integer(input$loss2)
    loss3 <- max(as.integer(input$loss3), loss2)
    updateSelectInput(session, "loss3", choices = loss_choices(loss2), selected = loss3)
  })
  
  # Optimal-strategies ranges: "to" is at least "from"; the third-line range is at
  # least the second-line range (from >= second-line from, to >= second-line to)
  opt_value <- function(id) as.numeric(input[[id]])
  update_to <- function(id, choices_from, min_value) {
    updateSelectInput(session, id, choices = choices_from(min_value),
                      selected = max(opt_value(id), min_value))
  }
  observeEvent(input$opt_loss2_from, {
    update_to("opt_loss2_to", loss_choices, opt_value("opt_loss2_from"))
  })
  observeEvent(list(input$opt_loss2_from, input$opt_loss2_to, input$opt_loss3_from), {
    req(input$opt_loss2_from, input$opt_loss2_to, input$opt_loss3_from, input$opt_loss3_to)
    loss3_from <- max(opt_value("opt_loss3_from"), opt_value("opt_loss2_from"))
    updateSelectInput(session, "opt_loss3_from", choices = loss_choices(opt_value("opt_loss2_from")),
                      selected = loss3_from)
    update_to("opt_loss3_to", loss_choices, max(loss3_from, opt_value("opt_loss2_to")))
  })
  observeEvent(input$opt_mds_from, {
    from <- opt_value("opt_mds_from")
    updateSelectInput(session, "opt_mds_to", choices = mds_choices[mds_choices >= from - 1e-9],
                      selected = max(opt_value("opt_mds_to"), from))
  })
  
  set_scenario <- function(loss2, loss3, mds) {
    updateSelectInput(session, "loss2", selected = loss2)
    updateSelectInput(session, "loss3", choices = loss_choices(loss2), selected = loss3)
    updateSelectInput(session, "mds", selected = mds)
  }
  observeEvent(input$preset_none, set_scenario(0, 0, 1))
  observeEvent(input$preset_max, set_scenario(50, 50, 3))

  design_id <- reactive({
    sc <- scenario()
    id <- design[loss2 == sc$loss2 & loss3 == sc$loss3 & abs(mds - sc$mds) < 1e-9, design_id]
    validate(need(length(id) == 1, "No simulated results for this combination."))
    id
  })
  scenario_outcomes <- reactive(outcomes[design_id == design_id()])
  scenario_clouds <- reactive(clouds[J(design_id())])

  # --- About ------------------------------------------------------------------
  output$strategy_table <- renderReactable({
    reactable(
      strategies[, .(Strategy = tx, `Treatment sequence` = name)],
      compact = TRUE, sortable = FALSE
    )
  })

  # --- Outcomes ---------------------------------------------------------------
  output$summary_table <- renderReactable({
    shown <- info[group == input$outcome_group]
    dt <- shown[, .(outcome, unit, higher_is_better)][scenario_outcomes(), on = "outcome", nomatch = NULL]
    dt[, cell := format_estimate(mean, lower95, upper95, unit)]
    dt[, best := fifelse(higher_is_better, mean == max(mean), mean == min(mean)), by = outcome]
    wide <- dcast(dt, tx ~ outcome, value.var = "cell")
    best <- dcast(dt, tx ~ outcome, value.var = "best")
    order <- match(strategy_order, wide$tx)
    wide <- wide[order]
    best <- best[order]
    wide <- add_strategy_info(wide)
    columns <- c(
      list(
        tx = colDef(name = "Strategy", minWidth = 80, sticky = "left",
                    style = function(value, index) list(color = position_palette[[as.character(wide$prrt_position[index])]],
                                                        fontWeight = 700)),
        prrt_position = colDef(show = FALSE), n_drugs = colDef(show = FALSE)
      ),
      setNames(lapply(shown$outcome, function(o) {
        colDef(name = sub("^(OS|PFS) at ", "At ", sub("MDS/AML cumulative incidence at", "At", shown[outcome == o, label])),
               minWidth = 135,
               style = function(value, index) {
                 c(list(whiteSpace = "nowrap"),
                   if (isTRUE(best[[o]][index])) list(fontWeight = 700, background = "#F2F2F2"))
               })
      }), shown$outcome)
    )
    reactable(
      wide[, c("tx", "prrt_position", "n_drugs", shown$outcome), with = FALSE],
      columns = columns,
      columnGroups = list(colGroup(name = input$outcome_group, columns = shown$outcome)),
      compact = TRUE, bordered = TRUE, sortable = FALSE
    )
  })

  output$interval_outcome_ui <- renderUI({
    selectInput("interval_outcome", NULL, choices = grouped_outcomes[[input$outcome_group]])
  })
  
  output$interval_plot <- renderPlot(height = height_from_width("interval_plot", 0.8), expr = {
    req(input$interval_outcome)
    o <- input$interval_outcome
    unit <- info[outcome == o, unit]
    plot_intervals(scenario_outcomes()[outcome == o],
                   paste0(info[outcome == o, label], ", ", unit_suffix(unit)))
  })

  # Downloads: estimate and 95% range of every outcome for every strategy
  results_table <- function(ids) {
    dt <- info[, .(outcome, label, unit)][outcomes[design_id %in% ids], on = "outcome"]
    dt <- design[dt, on = "design_id"]
    dt <- strategies[, .(tx, sequence = name)][dt, on = "tx"]
    dt[, tx := factor(tx, levels = strategy_order)]
    setorder(dt, design_id, tx)
    dt[, .(efficacy_loss_second_pct = loss2, additional_loss_third_pct = loss3, mds_aml_risk_multiplier = mds,
           strategy = as.character(tx), sequence, outcome = label, unit,
           estimate = mean, lower95 = lower95, upper95 = upper95)]
  }
  output$download_scenario <- downloadHandler(
    filename = function() {
      sc <- scenario()
      sprintf("NEC_outcomes_loss2-%d_loss3-%d_mds-%s.csv", sc$loss2, sc$loss3, sc$mds)
    },
    content = function(file) fwrite(results_table(design_id()), file)
  )
  output$download_all <- downloadHandler(
    filename = function() "NEC_outcomes_all_scenarios.csv",
    content = function(file) fwrite(results_table(design[loss3 >= loss2, design_id]), file)
  )
  
  # --- Trade-off --------------------------------------------------------------
  output$tradeoff_plot <- renderPlot(height = height_from_width("tradeoff_plot", 0.6), expr = {
    y_var <- input$tradeoff_y
    cloud <- scenario_clouds()[, .(tx = as.character(tx), x = mds_aml_5y_cir * 100, y = get(y_var))]
    plot_cloud(cloud, "MDS/AML at 5 years, %",
               if (y_var == "os_mean") "Mean overall survival, years" else "Mean progression-free survival, years")
  })

  # --- Compare ----------------------------------------------------------------
  output$compare_header <- renderText(paste0("Difference from ", input$reference, " (",
                                             strategies[tx == input$reference, name], ")"))

  output$difference_plot <- renderPlot(height = height_from_width("difference_plot", 0.8), expr = {
    y_var <- input$compare_y
    diff <- difference_cloud(scenario_clouds(), input$reference, y_var)
    outcome_name <- if (y_var == "os_mean") "mean overall survival" else "mean progression-free survival"
    plot_cloud(diff, paste0("Difference in MDS/AML at 5 years vs ", input$reference, ", percentage points"),
               paste0("Difference in ", outcome_name, " vs ", input$reference, ", months"),
               reference = input$reference) +
      ggplot2::annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.5, colour = "grey40",
                        label = "Longer survival, less MDS/AML") +
      ggplot2::annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.8, colour = "grey40",
                        label = "Shorter survival, more MDS/AML")
  })

  output$difference_table <- renderReactable({
    o <- input$compare_outcome
    dt <- vs_reference[design_id == design_id() & reference == input$reference & outcome == o]
    unit <- dt$diff_unit[1]
    dt[, difference := format_estimate(mean_diff, lower95, upper95, unit)]
    dt <- add_strategy_info(dt)
    dt[, tx := factor(tx, levels = strategy_order)]
    setorder(dt, tx)
    reactable(
      dt[, .(Strategy = as.character(tx), Difference = difference, prrt_position)],
      columns = list(
        Strategy = colDef(style = function(value, index) list(color = position_palette[[as.character(dt$prrt_position[index])]],
                                                              fontWeight = 700)),
        Difference = colDef(name = paste0("Difference vs ", input$reference, ", ", unit, " (95% range)"),
                            minWidth = 200),
        prrt_position = colDef(show = FALSE)
      ),
      compact = TRUE, bordered = TRUE, sortable = TRUE
    )
  })

  # --- Optimal strategies -----------------------------------------------------
  output$opt_margin_ui <- renderUI({
    unit <- info[outcome == input$opt_outcome, unit]
    choices <- if (unit == "years") {
      c("None" = 0, "1 month" = 1, "2 months" = 2, "3 months" = 3)
    } else {
      c("None" = 0, "1 percentage point" = 1, "2 percentage points" = 2, "5 percentage points" = 5)
    }
    radioButtons("opt_margin", "Clinical indifference margin", choices = choices, selected = 0)
  })

  optimal <- reactive({
    req(input$opt_margin, input$opt_loss2_from, input$opt_loss2_to, input$opt_loss3_from,
        input$opt_loss3_to, input$opt_mds_from, input$opt_mds_to)
    o <- input$opt_outcome
    unit <- info[outcome == o, unit]
    higher_better <- info[outcome == o, higher_is_better]
    margin <- as.numeric(input$opt_margin) / if (unit == "years") 12 else 1
    ids <- design[loss3 >= loss2 &
                    loss2 >= opt_value("opt_loss2_from") & loss2 <= opt_value("opt_loss2_to") &
                    loss3 >= opt_value("opt_loss3_from") & loss3 <= opt_value("opt_loss3_to") &
                    mds >= opt_value("opt_mds_from") - 1e-9 & mds <= opt_value("opt_mds_to") + 1e-9, design_id]
    validate(need(length(ids) > 0, "No combinations in this range: widen the ranges on the left."))
    m <- outcomes[outcome == o & design_id %in% ids]
    ranges <- optimal_strategies(m, higher_is_better = higher_better, margin = margin)
    list(ranges = ranges, m = m, unit = unit, label = info[outcome == o, label],
         higher_better = higher_better, n_settings = length(ids))
  })

  output$optimal_table <- renderReactable({
    res <- optimal()
    dt <- add_strategy_info(res$ranges[, .(tx = as.character(tx), worst, best, maximality,
                                           pessimism, optimism)])
    tick <- function(value) if (isTRUE(value)) "✔" else ""
    tick_style <- list(color = "#1B7837", fontWeight = 700, fontSize = "1.2em", textAlign = "center")
    reactable(
      dt[, .(Strategy = tx, worst = format_value(worst, res$unit),
             best = format_value(best, res$unit), maximality, pessimism, optimism, prrt_position)],
      columns = list(
        Strategy = colDef(minWidth = 70, style = function(value, index) list(
          color = position_palette[[as.character(dt$prrt_position[index])]], fontWeight = 700)),
        worst = colDef(name = paste0("Worst case, ", unit_suffix(res$unit))),
        best = colDef(name = paste0("Best case, ", unit_suffix(res$unit))),
        maximality = colDef(name = "Not uniformly worse", cell = tick, style = tick_style),
        pessimism = colDef(name = "Best worst case", cell = tick, style = tick_style),
        optimism = colDef(name = "Best best case", cell = tick, style = tick_style),
        prrt_position = colDef(show = FALSE)
      ),
      columnGroups = list(colGroup(name = paste0(res$label, " across ", res$n_settings, " combinations"),
                                   columns = c("worst", "best"))),
      compact = TRUE, bordered = TRUE, sortable = FALSE
    )
  })

  output$range_plot <- renderPlot(height = height_from_width("range_plot", 0.8), expr = {
    res <- optimal()
    dt <- res$ranges[, .(tx = as.character(tx), lower95 = pmin(worst, best), upper95 = pmax(worst, best))]
    plot_intervals(dt, paste0(res$label, ", ", unit_suffix(res$unit)), show_point = FALSE)
  })

  # --- Methods ----------------------------------------------------------------
  output$inputs_table <- renderReactable({
    reactable(inputs, columns = list(transition = colDef(name = "Transition", minWidth = 200),
                                     estimate = colDef(name = "Estimate", minWidth = 160),
                                     source = colDef(name = "Source", minWidth = 300)),
              compact = TRUE, bordered = TRUE, sortable = FALSE)
  })
}

shinyApp(ui, server)
