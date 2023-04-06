# adoptr designs used for testing
designad <- get_example_design()
designgs <- adoptr::GroupSequentialDesign(
  n1 = 29.53980042851903320411,
  c1f =                        0.8563037186428685831885,
  c1e = 2.211178640465977007779,
  n2_pivots = 25.97660166032095929722,
  c2_pivots = c(
    2.10305739106749,
    1.95299501057497,
    1.71048548318654,
    1.41728822790258,
    1.12411036243437,
    0.881629343021562,
    0.731701254448486
  ),
  7
) |>
  cache_design_splines()


