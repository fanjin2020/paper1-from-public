args <- commandArgs(TRUE)
if (length(args) != 1L) stop('Usage: Rscript publication_calibration/gaussian_geometry.R NEW_OUTPUT.csv')
if (file.exists(args[1])) stop('Output already exists; choose a new file')
if (!file.exists('publication_calibration/tail_design.json')) stop('Run from the repository root')
source('analysis/common.R', encoding='UTF-8')
source('paired_count_calibration/engine.R', encoding='UTF-8')
p <- jsonlite::read_json('publication_calibration/tail_design.json', simplifyVector=FALSE)
result <- lapply(p$scenarios, function(s) {
  design <- make_design(s)
  n <- s$n
  a0 <- solve(crossprod(design$x0), t(design$x0))
  a1 <- solve(crossprod(design$x1), t(design$x1))
  contrast <- a1[2, ] - a0[2, ]
  H <- design$x1 %*% a1
  M <- diag(n) - H
  B <- M %*% diag(contrast^2 / (1-diag(H))^2) %*% M
  tr <- sum(diag(B))
  data.frame(scenario=s$id, n=n, rho=s$rho, reference_df=n-ncol(design$x1),
             mean_HC3_variance_over_true=tr/sum(contrast^2),
             moment_equivalent_df=tr^2/sum(B*B),
             numerator_residual_orthogonality=max(abs(M %*% contrast)),
             max_leverage=max(diag(H)))
})
write.csv(do.call(rbind, result), args[1], row.names=FALSE)
