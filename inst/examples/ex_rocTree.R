data(simDat)

## Fitting a pruned survival tree
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = FALSE)

## Fitting a unpruned survival tree
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = FALSE,
        control = list(numFold = 0))

\dontrun{
## Fitting the ensemble algorithm (default)
rocTree(Surv(Time, death) ~ z1 + z2, id = id, data = simDat, ensemble = TRUE)
}
