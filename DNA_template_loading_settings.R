library(deSolve)

g0.range <- seq(0, 2E-9, 1E-10)

df <- data.frame(g0 = numeric(),
                 z0 = numeric())

n <- 1

for (i in g0.range) {
  init.values <- vector(mode = 'numeric')
  init.values <- append(init.values, c(i))
  
  state.variables <- c(gp = i, p0 = 1E-9, xs = 0, m0 = 0,
                       r0 = 1E-9, mr = 0, AA = 1E-9, NT = 1E-9,
                       z0 = 0)
  
  
  parm.values <- c(k1 = 6E9, k2 = 600, k3 = 1E10,
                   k4 = 6E9, k5 = 135, k6 = 1E10,
                   k7 = 18)
  
  con.profile <- function(t, state.variables, parm.values){
    with(as.list(c(state.variables, parm.values)),{
      
      dgp <- -k1*gp*p0 + k2*xs + k3*NT*xs
      dp0 <- -k1*gp*p0 + k2*xs + k3*NT*xs
      dxs <-  k1*gp*p0 - k2*xs - k3*NT*xs
      dm0 <-  k3*NT*xs - k4*m0*r0 + k5*mr + k6*AA*mr - k7*m0
      dr0 <- -k4*m0*r0 + k5*mr + k6*AA*mr
      dmr <-  k4*m0*r0 - k5*mr - k6*AA*mr
      dAA <- -k6*AA*mr
      dNT <- -k3*NT*xs
      dz0 <-  k6*AA*mr
      
      return(list(c(dgp, dp0, dxs, dm0, dr0, dmr, dAA,
                    dNT, dz0)))
    })
  }
  
  out <- ode(y = state.variables,
             parms = parm.values,
             method = "bdf",
             times = seq(0, 30, 0.001),
             func = con.profile)
  
  zx.value <- unname(out[nrow(out), ncol(out)])
  
  init.values <- append(init.values, zx.value)
  
  df[n,] <- init.values
  
  n <- n + 1

}
