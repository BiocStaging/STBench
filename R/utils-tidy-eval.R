# ============================================================
# Tidy evaluation and global variables
# ============================================================

#' @importFrom rlang .data .env
NULL

#' @importFrom magrittr %>%
NULL

# Declare column names used in dplyr and ggplot2
# non-standard evaluation.

utils::globalVariables(
  c(
    "LR_strength",
    "communication_strength",
    "cutoff",
    "interaction_name",
    "mean_strength",
    "method",
    "pathway_name",
    "pathway_strength",
    "prob",
    "removed",
    "target",
    "type",
    "value"
  )
)
