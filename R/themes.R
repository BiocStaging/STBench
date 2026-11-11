
# Plot themes



.get_plot_theme <- function(
    theme_name = "classic"
) {

  switch(
    theme_name,

    classic =
      ggplot2::theme_classic(),

    bw =
      ggplot2::theme_bw(),

    minimal =
      ggplot2::theme_minimal(),

    gray =
      ggplot2::theme_gray(),

    stop(
      "Unknown theme: ",
      theme_name
    )
  )
}
