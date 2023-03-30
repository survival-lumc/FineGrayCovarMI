# Add parameters for MI imps per kmi imp??
smcfcs.finegray <- function(originaldata,
                            smtype,
                            smformula,
                            method,
                            m = 5,
                            numit = 10,
                            rjlimit = 5000,
                            ...) {

  kmi_imps <- kmi::kmi(
    Surv(time, D != 0) ~ 1, # or > 0, add cens code
    data = data.frame(originaldata),
    etype = D,
    failcode = 1,
    nimp = m
  )

  kmi_imps

  # Do invisible and capture output bits here..

}
