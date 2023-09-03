library(rpact)

designIN <- getDesignInverseNormal(kMax = 2, futilityBounds = c(0))


datasetExample <- getDataset(
  means1 = c(112.3, 113.1),
  means2 = c(98.1, 99.3),
  stDevs1 = c(44.4, 42.9),
  stDevs2 = c(46.7, 41.1),
  n1 = c(34, 31),
  n2 = c(37, 33)
)

results <- getAnalysisResults(design = designIN, dataInput = datasetExample, stage = 2)
results



results <- getAnalysisResults(
  design = designIN, datasetExample,
  stage = 2, nPlanned = 60
)

