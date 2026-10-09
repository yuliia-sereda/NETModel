# igraph is only needed for the Hasse diagrams, so it is a suggested package.
check_igraph <- function() {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for Hasse diagrams: install.packages(\"igraph\")")
  }
}

#' Reduce a directed graph to its Hasse diagram
#'
#' @description
#' Removes all transitive edges from a directed acyclic graph (DAG) while
#' preserving the underlying partial order. The resulting graph is the Hasse
#' diagram of the original relation, containing only the cover relations
#' necessary to represent the ordering.
#' @param g a directed graph of relations
#' @return a reduced (hasse) graph of the relations)
#' 
#' @export
reduce_graph <- function(g) {
  check_igraph()
  stopifnot(igraph::is_dag(g))
  if (igraph::ecount(g) == 0L) {
    return(g)
  }
  
  A <- igraph::as_adjacency_matrix(g, sparse = FALSE) > 0
  # R[w, v]: v is reachable from w by a path of length >= 1
  R <- is.finite(igraph::distances(g, mode = "out")) & !diag(igraph::vcount(g))
  # An edge u -> v is transitive if u -> w -> ... -> v for some other w.
  transitive <- A & (A %*% R > 0)
  
  edges <- igraph::as_edgelist(g, names = FALSE)
  edges_to_remove <- which(transitive[edges])
  return(igraph::delete_edges(g, edges = edges_to_remove))
}


#' Create a Hasse diagram from an adjacency matrix
#'
#' @description
#' Constructs a directed graph from an adjacency matrix representing a binary
#' relation and removes all transitive edges to obtain the corresponding Hasse
#' diagram.
#' @param adjmat an adjacency matrix of relations
#' @return a hasse graph of the relations
#' 
#' @export
hasse_graph <- function(adjmat) {
  check_igraph()
  g <- reduce_graph(igraph::graph_from_adjacency_matrix(adjmatrix = adjmat, mode = "directed"))
  return(g)
}

#' Create a hierarchical layout for a connected directed graph
#'
#' @description
#' Assigns two-dimensional coordinates to the vertices of a connected directed
#' graph for use in a Hasse-style hierarchical plot. Vertices with zero
#' in-degree are placed in the first level, and subsequent levels are formed
#' iteratively from their outgoing neighbors.
#'
#' Vertices within the same level are spaced evenly along the horizontal axis,
#' while the vertical coordinate indicates the graph level. The function
#' assumes that `g` is connected and that vertex names are available through
#' the `"name"` vertex attribute.
#'
#' @param g An `igraph` object representing a connected directed graph.
#'
#' @return A numeric matrix with one row per vertex and two columns, `"x"` and
#'   `"y"`. Row names correspond to vertex names. The `"x"` coordinate gives
#'   the within-level horizontal position and the `"y"` coordinate gives the
#'   hierarchical level.
#'
#' @param horizontal_spacing Horizontal distance between vertices.
#' @param vertical_spacing Vertical distance between levels.
#'
#' @keywords internal
get_layout_connected_subgraph <- function(
    g,
    horizontal_spacing = 2,
    vertical_spacing = 3
) {
  stopifnot(igraph::is_dag(g))
  
  vertex_names <- igraph::V(g)$name
  
  if (is.null(vertex_names)) {
    stop("The graph must have a vertex 'name' attribute.")
  }
  
  # Handle an isolated vertex
  if (igraph::vcount(g) == 1L) {
    out <- matrix(
      c(0, 0),
      nrow = 1,
      dimnames = list(vertex_names, c("x", "y"))
    )
    return(out)
  }
  
  # Assign levels using the longest path from any source.
  # This guarantees that every child is above all its parents.
  topo_names <- igraph::as_ids(igraph::topo_sort(g, mode = "out"))
  
  level <- stats::setNames(
    rep(0L, igraph::vcount(g)),
    vertex_names
  )
  
  for (v in topo_names) {
    parent_names <- igraph::as_ids(
      igraph::neighbors(g, v, mode = "in")
    )
    
    if (length(parent_names) > 0L) {
      level[v] <- max(level[parent_names]) + 1L
    }
  }
  
  names_per_level <- split(names(level), level)
  
  layout <- matrix(
    NA_real_,
    nrow = igraph::vcount(g),
    ncol = 2,
    dimnames = list(vertex_names, c("x", "y"))
  )
  
  # Initial ordering within each level
  for (i in seq_along(names_per_level)) {
    current_names <- sort(names_per_level[[i]])
    n <- length(current_names)
    
    x <- seq_len(n) * horizontal_spacing
    x <- x - mean(x)
    
    layout[current_names, "x"] <- x
    layout[current_names, "y"] <- as.numeric(names(names_per_level)[i]) *
      vertical_spacing
  }
  
  # Reorder each non-source level using the average x-position
  # of its parents. This generally reduces edge crossings.
  if (length(names_per_level) > 1L) {
    for (i in 2:length(names_per_level)) {
      current_names <- names_per_level[[i]]
      
      parent_centres <- vapply(
        current_names,
        function(v) {
          parent_names <- igraph::as_ids(
            igraph::neighbors(g, v, mode = "in")
          )
          
          if (length(parent_names) == 0L) {
            0
          } else {
            mean(layout[parent_names, "x"])
          }
        },
        numeric(1)
      )
      
      current_names <- current_names[
        order(parent_centres, current_names)
      ]
      
      n <- length(current_names)
      x <- seq_len(n) * horizontal_spacing
      x <- x - mean(x)
      
      layout[current_names, "x"] <- x
    }
  }
  
  layout
}

#' Create a layout for a Hasse diagram
#'
#' @description
#' Generates plotting coordinates for a directed acyclic graph. Each weakly
#' connected component is laid out hierarchically. All terminal vertices,
#' defined as vertices with zero out-degree, are placed on the same vertical
#' level and spaced evenly along the horizontal axis.
#'
#' @param g A directed acyclic `igraph` object with a vertex `"name"` attribute.
#' @param horizontal_spacing Horizontal distance between vertices.
#' @param vertical_spacing Vertical distance between hierarchical levels.
#' @param component_gap Horizontal gap between nontrivial graph components.
#'
#' @return A numeric matrix with columns `"x"` and `"y"`. Rows are ordered to
#'   match `V(g)`.
#'
#' @export
hasse_layout <- function(
    g,
    horizontal_spacing = 2,
    vertical_spacing = 2,
    component_gap = 2
) {
  check_igraph()
  vertex_names <- igraph::V(g)$name
  
  if (is.null(vertex_names)) {
    stop("The graph must have vertex names.")
  }
  
  graph_list <- igraph::decompose(g, mode = "weak")
  
  # Separate isolated vertices from components containing edges
  isolated <- vapply(
    graph_list,
    igraph::vcount,
    numeric(1)
  ) == 1L
  
  connected_graphs <- graph_list[!isolated]
  isolated_graphs  <- graph_list[isolated]
  
  L <- NULL
  right_edge <- 0
  
  # Lay out and combine non-isolated components
  for (subgraph in connected_graphs) {
    tmp <- get_layout_connected_subgraph(subgraph)
    
    # Increase spacing within the component
    tmp[, "x"] <- tmp[, "x"] * horizontal_spacing
    tmp[, "y"] <- tmp[, "y"] * vertical_spacing
    
    # Shift component to the right of previous components
    if (!is.null(L)) {
      shift <- right_edge + component_gap - min(tmp[, "x"])
      tmp[, "x"] <- tmp[, "x"] + shift
    }
    
    L <- rbind(L, tmp)
    right_edge <- max(L[, "x"])
  }
  
  # Highest level reached by the connected graph
  terminal_y <- if (is.null(L)) {
    0
  } else {
    max(L[, "y"])
  }
  
  # Add isolated vertices at the terminal level
  if (length(isolated_graphs) > 0L) {
    isolated_names <- vapply(
      isolated_graphs,
      function(x) igraph::V(x)$name,
      character(1)
    )
    
    isolated_layout <- cbind(
      x = rep(0, length(isolated_names)),
      y = rep(terminal_y, length(isolated_names))
    )
    
    rownames(isolated_layout) <- isolated_names
    
    L <- if (is.null(L)) {
      isolated_layout
    } else {
      rbind(L, isolated_layout)
    }
  }
  
  # Detect all terminal vertices automatically
  out_degree <- igraph::degree(g, mode = "out")
  terminal_vertices <- names(out_degree)[out_degree == 0]
  
  # Keep connected terminals ordered by their current position.
  # Isolated terminals are ordered by name.
  in_degree <- igraph::degree(g, mode = "in")
  
  connected_terminals <- terminal_vertices[
    in_degree[terminal_vertices] > 0
  ]
  
  isolated_terminals <- terminal_vertices[
    in_degree[terminal_vertices] == 0
  ]
  
  connected_terminals <- connected_terminals[
    order(L[connected_terminals, "x"])
  ]
  
  isolated_terminals <- sort(isolated_terminals)
  
  # Put half the isolated nodes on each side
  n_left <- floor(length(isolated_terminals) / 2)
  
  left_isolated <- if (n_left > 0) {
    isolated_terminals[seq_len(n_left)]
  } else {
    character(0)
  }
  
  right_isolated <- if (n_left < length(isolated_terminals)) {
    isolated_terminals[
      seq.int(n_left + 1, length(isolated_terminals))
    ]
  } else {
    character(0)
  }
  
  terminal_order <- c(
    left_isolated,
    connected_terminals,
    right_isolated
  )
  
  # Evenly space all terminal vertices along the top level
  terminal_x <- seq(
    from = 0,
    by = horizontal_spacing,
    length.out = length(terminal_order)
  )
  
  terminal_x <- terminal_x - mean(terminal_x)
  
  L[terminal_order, "x"] <- terminal_x
  L[terminal_order, "y"] <- terminal_y
  
  # Center the entire layout
  L[, "x"] <- L[, "x"] - mean(range(L[, "x"]))
  
  # Return in original vertex order
  L[vertex_names, , drop = FALSE]
}

#' Create an adjacency matrix for maximality
#'
#' @description
#' Constructs an adjacency matrix representing pairwise treatment preferences
#' based on simulated differences in an outcome.
#' 
#' @param bounds a data.table of bounds
#' @param pw_relation a function that takes two bounds (and other arguments) and returns a boolean encoding the relation
#' @param out_var the name of the outcome variable
#' @param digits the number of digits to round the comparison of the bounds to
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param reverse whether to reverse the bounds (use for outcomes where lower is better)
#' @param alpha the alpha value for the hurwicz function
#' @return an adjacency matrix of relations
#' 
#' @export
adjacency_matrix_from_bounds <- function(
    bounds, 
    pw_relation, 
    out_var = "os_mean", 
    digits = 3, 
    added_margin = 0, 
    reverse = FALSE, 
    alpha = NULL # hurwicz only 
) {
  # Column names used in data.table expressions (declared for R CMD check)
  res_min <- res_max <- NULL
  
  
  stopifnot(length(out_var) == 1)

  bound_cols <- paste0(out_var, c("..min", "..max"))
  stopifnot(all(bound_cols %in% names(bounds)))
  dat <- bounds[, c("treatment", bound_cols), with = FALSE]
  
  if(isFALSE(reverse)) {
    setnames(dat, old = bound_cols, new = c("res_min", "res_max"))
  } else {
    setnames(dat, old = bound_cols, new = c("res_max", "res_min"))
    dat[, `:=`(res_min = -res_min, res_max = - res_max)]
  }
  dat[, `:=`(
    res_min = round(res_min, digits = digits), 
    res_max = round(res_max, digits = digits))]
  
  
  mat <- matrix(rep(0, nrow(dat)^2), nrow = nrow(dat))
  colnames(mat) <- dat$treatment
  rownames(mat) <- dat$treatment
  
  for(row in seq_len(nrow(mat))){
    for(col in seq_len(ncol(mat))){
      if (row != col) {
        mat[row, col] <- pw_relation(
          tx1_bounds = dat[col, c(res_min, res_max)], 
          tx2_bounds = dat[row, c(res_min, res_max)], 
          reverse = FALSE, 
          added_margin = added_margin,
          alpha = alpha)  
      }
    }
  }
  
  return(mat)
}


#' Create an adjacency matrix for maximality from differences
#' 
#' @description
#' Constructs an adjacency matrix representing pairwise treatment preferences
#' under the maximality criterion using simulated pairwise outcome differences.
#' For each ordered pair of treatments, the smallest simulated difference is
#' computed across all uncertainty scenarios. A directed edge from treatment A
#' to treatment B is included when the minimum difference (A - B) exceeds the
#' specified decision margin, indicating that A is preferred to B under the
#' maximality criterion.
#' 
#' @param differences a data.table of differences
#' @param out_var the name of the outcome variable
#' @param digits the number of digits to round the comparison of the bounds to
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param reverse whether to reverse the differences (use for outcomes where lower is better)
#' @return an adjacency matrix of relations
#' 
#' @export
adjacency_matrix_maximality <- function(
    differences, 
    out_var = "os_mean", 
    digits = 3, 
    added_margin = 0, 
    reverse = FALSE
) {
  # Column names used in data.table expressions (declared for R CMD check)
  tx1 <- tx2 <- delta <- NULL
  
  stopifnot(length(out_var) == 1)
  delta_col <- paste0(out_var, "_delta")
  stopifnot(delta_col %in% names(differences))
  
  sign <- if (isTRUE(reverse)) -1 else 1
  dat <- differences[, list(tx1, tx2, delta = sign * round(get(delta_col), digits = digits))]
  
  # For each pair, the worst-case difference in each direction
  dt <- dat[, list(min_delta_12 = min(delta), min_delta_21 = min(-delta)), by = c("tx1", "tx2")]
  
  tx_names <- unique(c(dt$tx1, dt$tx2))
  mat <- matrix(FALSE, nrow = length(tx_names), ncol = length(tx_names),
                dimnames = list(tx_names, tx_names))
  
  # Edge from the dominated to the preferred treatment
  i1 <- match(dt$tx1, tx_names)
  i2 <- match(dt$tx2, tx_names)
  mat[cbind(i2, i1)] <- dt$min_delta_12 > added_margin
  mat[cbind(i1, i2)] <- dt$min_delta_21 > added_margin
  return(mat)
}

#' Calculate uncertainty bounds
#'
#' @description
#' Computes lower and upper bounds for one or more simulated outcomes by
#' treatment. For each outcome, the minimum and maximum simulated values are
#' returned and can subsequently be used for interval-based decision analyses.
#' 
#' @param predictions a data.table of predictions
#' @param out_vars the names of the outcome variables
#' @param by_vars the names of the by variables
#' @return a data.table of bounds
#' 
#' @export
get_bounds <- function(
    predictions, 
    out_vars = c("os_mean", "os_median", "pfs_mean", "pfs_median", "mds_aml_5y_cir"), 
    by_vars = "treatment") {
  pred_min <- predictions[, lapply(.SD, function(x) min(x)), 
                          .SDcols = out_vars, 
                          by = by_vars]
  setnames(x = pred_min, old = out_vars, new = c(paste0(out_vars, "..min")))
  
  pred_max <- predictions[, lapply(.SD, function(x) max(x)),
                          .SDcols = out_vars, 
                          by = by_vars]
  setnames(x = pred_max, old = out_vars, new = c(paste0(out_vars, "..max")))
  return(pred_min[pred_max, on = by_vars])
}


#' Optimistic pairwise preference relation
#'
#' @description
#' Determines whether one treatment is preferred over another under the
#' optimism criterion by comparing the upper bounds of their uncertainty
#' intervals.
#' 
#' @param tx1_bounds the bounds of the first treatment
#' @param tx2_bounds the bounds of the second treatment
#' @param reverse whether to reverse the bounds
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param ... other arguments
#' @return a boolean encoding the relation
#' 
#' @export
pairwise_relation_optimism <- function(
    tx1_bounds, 
    tx2_bounds, 
    reverse = FALSE, 
    added_margin = 0,
    ...    
){
  stopifnot(is.vector(tx1_bounds) && length(tx1_bounds) == 2)
  stopifnot(is.vector(tx2_bounds) && length(tx2_bounds) == 2)
  if (isTRUE(reverse)) {
    tx1_bounds <- -tx1_bounds[c(2,1)]
    tx2_bounds <- -tx2_bounds[c(2,1)]
  }
  result <- FALSE 
  if(tx1_bounds[2] > tx2_bounds[2] + added_margin) result <- TRUE
  return(result)
}


#' Pessimistic pairwise preference relation
#'
#' @description
#' Determines whether one treatment is preferred over another under the
#' pessimism criterion by comparing the lower bounds of their uncertainty
#' intervals.
#' 
#' @param tx1_bounds the bounds of the first treatment
#' @param tx2_bounds the bounds of the second treatment
#' @param reverse whether to reverse the bounds
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param ... other arguments 
#' @return a boolean encoding the relation
#' 
#' @export
pairwise_relation_pessimism <- function(
    tx1_bounds, 
    tx2_bounds, 
    reverse = FALSE, 
    added_margin = 0,
    ...
){
  stopifnot(is.vector(tx1_bounds) && length(tx1_bounds) == 2)
  stopifnot(is.vector(tx2_bounds) && length(tx2_bounds) == 2)
  if (isTRUE(reverse)) {
    tx1_bounds <- -tx1_bounds[c(2,1)]
    tx2_bounds <- -tx2_bounds[c(2,1)]
  }
  result <- FALSE 
  if(tx1_bounds[1] > tx2_bounds[1] + added_margin) result <- TRUE
  return(result)
}

#' Hurwicz pairwise preference relation
#'
#' @description
#' Determines whether one treatment is preferred over another using the
#' Hurwicz decision criterion, which evaluates each uncertainty interval as a
#' weighted average of its lower and upper bounds. The weighting parameter
#' `alpha` controls the degree of optimism, with `alpha = 0` corresponding to
#' the pessimism criterion and `alpha = 1` corresponding to the optimism
#' criterion.
#' 
#' @param tx1_bounds the bounds of the first treatment
#' @param tx2_bounds the bounds of the second treatment
#' @param reverse whether to reverse the bounds
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param alpha the alpha value for the hurwicz function. 0 is same as pessimism, 1 is same as optimism. 
#' @return a boolean encoding the relation
#' 
#' @export
pairwise_relation_hurwicz <- function(
    tx1_bounds, 
    tx2_bounds, 
    reverse = FALSE, 
    added_margin = 0,
    alpha = NULL
){
  stopifnot(is.vector(tx1_bounds) && length(tx1_bounds) == 2)
  stopifnot(is.vector(tx2_bounds) && length(tx2_bounds) == 2)
  stopifnot(!is.null(alpha) && alpha >=0 && alpha <=1)
  
  if (isTRUE(reverse)) {
    tx1_bounds <- -tx1_bounds[c(2,1)]
    tx2_bounds <- -tx2_bounds[c(2,1)]
  }
  result <- FALSE 
  if((1-alpha) * tx1_bounds[1] + alpha * tx1_bounds[2] > 
     (1-alpha) * tx2_bounds[1] + alpha * tx2_bounds[2] + added_margin) {
    result <- TRUE
  }
  return(result)
}

#' Interval-dominance pairwise preference relation
#'
#' @description
#' Determines whether one treatment interval strictly dominates another by
#' requiring the lower bound of the first interval to exceed the upper bound
#' of the second interval (after applying any decision margin).
#' 
#' @param tx1_bounds the bounds of the first treatment
#' @param tx2_bounds the bounds of the second treatment
#' @param reverse whether to reverse the bounds
#' @param added_margin the margin to add to the comparison (e.g. 2/12 for 2 months when time is in years)
#' @param ... other arguments 
#' @return a boolean encoding the relation
#' 
#' @export
pairwise_relation_interval_dominance <- function(
    tx1_bounds, 
    tx2_bounds, 
    reverse = FALSE, 
    added_margin = 0,
    ...
){
  stopifnot(is.vector(tx1_bounds) && length(tx1_bounds) == 2)
  stopifnot(is.vector(tx2_bounds) && length(tx2_bounds) == 2)
  if (isTRUE(reverse)) {
    tx1_bounds <- -tx1_bounds[c(2,1)]
    tx2_bounds <- -tx2_bounds[c(2,1)]
  }
  result <- FALSE 
  if(tx1_bounds[1] > tx2_bounds[2] + added_margin) result <- TRUE
  return(result)
}
