# Tests for the choice functions and Hasse diagrams

toy_bounds <- data.table::data.table(
  treatment = c("A", "B", "C"),
  os_mean..min = c(5.0, 3.0, 4.5),
  os_mean..max = c(6.0, 4.0, 5.5),
  # A column whose name contains the outcome must not be picked up
  os_mean_se..min = c(0, 0, 0)
)

test_that("pairwise relations reject bounds that are not of length 2", {
  expect_error(pairwise_relation_optimism(c(1, 2, 3), c(1, 2)))
  expect_error(pairwise_relation_pessimism(c(1, 2), 1))
  expect_error(pairwise_relation_interval_dominance(c(1, 2), c(1, 2, 3)))
  expect_error(pairwise_relation_hurwicz(c(1, 2, 3), c(1, 2), alpha = 0.5))
})

test_that("adjacency_matrix_from_bounds() points edges from the worse to the better treatment", {
  m <- adjacency_matrix_from_bounds(toy_bounds, pairwise_relation_interval_dominance, out_var = "os_mean")
  # A and C both dominate B; A and C overlap
  expect_true(m["B", "A"] == 1 && m["B", "C"] == 1)
  expect_true(m["A", "B"] == 0 && m["C", "B"] == 0)
  expect_true(m["A", "C"] == 0 && m["C", "A"] == 0)
  # For an outcome where lower is better, the direction flips
  r <- adjacency_matrix_from_bounds(toy_bounds, pairwise_relation_interval_dominance,
                                    out_var = "os_mean", reverse = TRUE)
  expect_true(r["A", "B"] == 1 && r["C", "B"] == 1 && r["B", "A"] == 0)
  # A margin larger than the gap removes the edge
  g <- adjacency_matrix_from_bounds(toy_bounds, pairwise_relation_interval_dominance,
                                    out_var = "os_mean", added_margin = 1)
  expect_equal(g["B", "A"], 0)
})

test_that("adjacency_matrix_maximality() uses the worst-case pairwise difference", {
  differences <- data.table::data.table(
    tx1 = "A", tx2 = "B",
    os_mean_delta = c(0.5, 0.2, 0.3),
    mds_aml_5y_cir_delta = c(0.01, 0.02, 0.03)
  )
  m <- adjacency_matrix_maximality(differences, out_var = "os_mean")
  expect_true(m["B", "A"])
  expect_false(m["A", "B"])
  # A 2-month margin (0.167) is smaller than the worst case (0.2); a 0.25 margin is not
  expect_true(adjacency_matrix_maximality(differences, out_var = "os_mean", added_margin = 2 / 12)["B", "A"])
  expect_false(adjacency_matrix_maximality(differences, out_var = "os_mean", added_margin = 0.25)["B", "A"])
  # Lower CIR is better: B is preferred to A
  r <- adjacency_matrix_maximality(differences, out_var = "mds_aml_5y_cir", reverse = TRUE)
  expect_true(r["A", "B"])
  expect_false(r["B", "A"])
})

test_that("get_bounds() defaults to outcome-first names", {
  predictions <- data.table::data.table(
    treatment = rep(c("A", "B"), each = 2),
    os_mean = 1:4, os_median = 1:4, pfs_mean = 1:4, pfs_median = 1:4, mds_aml_5y_cir = 1:4
  )
  b <- get_bounds(predictions)
  expect_true(all(paste0(rep(c("os_mean", "os_median", "pfs_mean", "pfs_median", "mds_aml_5y_cir"), each = 2),
                         c("..min", "..max")) %in% names(b)))
  expect_equal(b[treatment == "B", os_mean..max], 4)
})

test_that("reduce_graph() removes exactly the transitive edges", {
  skip_if_not_installed("igraph")
  # A -> B -> C -> D plus shortcuts A -> C, A -> D, B -> D, and a separate E -> D
  adj <- matrix(0, 5, 5, dimnames = list(LETTERS[1:5], LETTERS[1:5]))
  adj["A", c("B", "C", "D")] <- 1
  adj["B", c("C", "D")] <- 1
  adj["C", "D"] <- 1
  adj["E", "D"] <- 1
  g <- hasse_graph(adj)
  edges <- sort(apply(igraph::as_edgelist(g), 1, paste, collapse = ">"))
  expect_equal(edges, c("A>B", "B>C", "C>D", "E>D"))
})

test_that("Hasse diagrams work without attaching igraph", {
  skip_if_not_installed("igraph")
  # hasse_layout() used to call decompose() unqualified, which resolves to
  # stats::decompose() unless igraph is attached.
  expect_false("package:igraph" %in% search())
  m <- adjacency_matrix_from_bounds(toy_bounds, pairwise_relation_interval_dominance, out_var = "os_mean")
  g <- hasse_graph(m)
  lay <- hasse_layout(g)
  expect_equal(dim(lay), c(3L, 2L))
  expect_equal(rownames(lay), igraph::V(g)$name)
})
