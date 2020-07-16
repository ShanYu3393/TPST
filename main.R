rm(list = ls())
data(Tr1)
data(V1)

ngrid.x=40
ngrid.y=20
ngrid.t=10

xx=seq(-0.89,3.39,length.out=ngrid.x)
yy=seq(-0.89,0.89,length.out=ngrid.y)
ss=expand.grid(xx,yy)
tt=(0:(ngrid.t-1))/(ngrid.t-1)

Data=data.frame(x=rep(ss[,1],ngrid.t),y=rep(ss[,2],ngrid.t),
                t=rep(tt,each=dim(ss)[1]))

knots=c(0.2,0.4,0.6,0.8)
Boundary.knots=c(0,1)

# library(BPST)
# library(splines2)
# library(MGLM)
# library(Matrix)
library(TPST)
d <- 2
r <- 1
rho <- 3
Basis1 <- basis.tensor(ss = ss, tt = tt, V = V1, Tri = Tr1,
                       d = d, r = r, time.knots = knots, rho = rho,
                       time.bound = Boundary.knots, line.up = TRUE)

Basis2 <- basis.tensor(ss = Data[,1:2], tt = Data[,3],
                       V = V1, Tri = Tr1, d = d, r = r,
                       time.knots = knots, rho = rho,
                       time.bound = Boundary.knots, line.up = FALSE)

which(Basis1$Psi.Q2 != Basis2$Psi.Q2)
