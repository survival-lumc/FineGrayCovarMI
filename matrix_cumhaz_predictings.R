

# Multiply whole cumhaz by other vector

cumhaz
HRs
HRs %*% diag(cumhaz)
matrix(r)

mat <- matrix(
  data = rep(cumhaz, times = length(HRs)),
  ncol = length(HRs)
)
mat %*% diag(HRs)

cumhaz <- seq_len(1000)
HRs <- seq_len(1000)

microbenchmark::microbenchmark(
  "mats" = {
    mat <- matrix(
      data = rep(cumhaz, times = length(HRs)),
      ncol = length(HRs)
    )
    mat %*% diag(HRs)
  },
  "sapply" = {
    sapply(HRs, function(hr) hr * cumhaz)
  },
  times = 10
)

# Matrix way is way slower for larger data (because cost of setting up the matrix)
