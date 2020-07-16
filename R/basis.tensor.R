#' Tensor-product Splines over a Triangular Prism Partition
#'
#' This function generates the basis for tensor-product splines over a triangular prism partition.
#'
#' @import Matrix
#' @import splines2
#' @import MGLM
#' @import BPST
#' @param ss If \code{line.up=TRUE}, \code{ss} is the spatial cooridinates of dimension \code{nS} by two. If \code{line.up=TRUE}, \code{ss} is the spatial cooridinates of dimension \code{n} by two. Each row is the spatial coordinates of a point.
#' \cr
#' @param tt If \code{line.up=TRUE}, the temporal cooridinates of length \code{nT}. If \code{line.up=FALSE}, the temporal cooridinates of length \code{nT}. Each value is the temporal coordinates of a point.
#' \cr
#' @param V The \code{N} by two matrix of vertices of a triangulation, where \code{N} is the number of vertices. Each row is the coordinates for a vertex.
#' \cr
#' @param Tri The triangulation matrix of dimention \code{nTr} by three, where \code{nTr} is the number of triangles in the triangulation. Each row is the indices of vertices in \code{V}.
#' \cr
#' @param d The degree of piecewise polynomials -- default is 2, and usually \code{d} is greater than one. -1 represents piecewise constant.
#' \cr
#' @param r The smoothness parameter -- default is 1, and 0 \eqn{\le} \code{r} \eqn{<} \code{d}.
#' \cr
#' @param time.knots The vector of interior time.knots for univariate spline.
#' \cr
#' @param time.bound The vector of two. The boundary of univariate spline.
#' \cr
#' @param rho The order of univaraite spline.
#' \cr
#' @param line.up The indicator of whether the observed points are temporally lined up or not -- default is \code{FALSE}.
#' \cr
#' @return A list of vectors and matrices, including:
#' \item{Psi}{The spline basis function of dimension \code{n} by \code{(l)}\code{(nTr)}\code{{(d+1)(d+2)/2}}, where \code{l} is the number of univariate basis, \code{n} is the number of observationed points, \code{nTr} is the number of triangles in the given triangulation, and \code{d} is the degree of the spline. If some points do not fall in the triangulation, the generation of the spline basis will not take those points into consideration.}
#' \item{Psi.Q2}{The spline basis function after QR decomposition}
#' \item{H}{The smoothness matrix for bivariate spline.}
#' \item{Q2}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H}.}
#' \item{H.all}{The smoothness matrix for tensor-product spline.}
#' \item{Q2.all}{The Q2 matrix after QR decomposition of the smoothness matrix \code{H.all} for tensor-product spline.}
#' \item{dimB}{The number of bivariate spline basis functions.}
#' \item{dimU}{The number of univariate spline basis functions.}
#' \item{P1, P2}{The penalty matrices from energy functions.}
#' \item{K1, K2}{The penalty matrices from energy functions QR decomposition.}
#' @examples
#' # load need libraries.
#' # Packages BPST and Triangulation could be downloaded from github.
#' rm(list = ls())
#' library(devtools)
#' install_github("funstatpackages/BPST")
#' install_github("funstatpackages/Triangulation")
#' library(BPST)
#' library(Triangulation)
#' library(splines2)
#' library(MGLM)
#' library(Matrix)
#' library(TPST)
#' data(Tr1)
#' data(V1)
#'
#' ngrid.x=40
#' ngrid.y=20
#' ngrid.t=10
#'
#' xx=seq(-0.89,3.39,length.out=ngrid.x)
#' yy=seq(-0.89,0.89,length.out=ngrid.y)
#' ss=expand.grid(xx,yy)
#' tt=(0:(ngrid.t-1))/(ngrid.t-1)
#'
#' Data=data.frame(x=rep(ss[,1],ngrid.t),y=rep(ss[,2],ngrid.t),
#'               t=rep(tt,each=dim(ss)[1]))
#' knots=c(0.2,0.4,0.6,0.8)
#' Boundary.knots=c(0,1)
#'
#' d <- 2
#' r <- 1
#' rho <- 3
#' Basis1 <- basis.tensor(ss = ss, tt = tt, V = V1, Tri = Tr1,
#'                      d = d, r = r, time.knots = knots, rho = rho,
#'                      time.bound = Boundary.knots, line.up = TRUE)
#'
#' Basis2 <- basis.tensor(ss = Data[,1:2], tt = Data[,3],
#'                       V = V1, Tri = Tr1, d = d, r = r,
#'                     time.knots = knots, rho = rho,
#'                     time.bound = Boundary.knots, line.up = FALSE)
#'
#' which(Basis1$Psi.Q2 != Basis2$Psi.Q2)
#' @export

basis.tensor <- function(ss, tt, V, Tri, d = 2, r = 1, time.knots,
                         time.bound, rho = 3, line.up = FALSE) {
  require(BPST)
  require(splines2)
  require(MGLM)
  require(Matrix)
  # To reduce the computation burden, we first identify unique spatial points and obtain basis functions
  # based on these points.
  ss <- as.matrix(ss)
  ss.uni <- unique(ss)
  index.ss <- match(data.frame(t(ss)), data.frame(t(ss.uni)))
  B.uni <- basis(V, Tri, d, r, ss.uni)
  BB.uni <- B.uni$B
  B0.uni <- matrix(NA, nrow(ss.uni), ncol(BB.uni))
  B0.uni[B.uni$Ind.inside, ] <- as.matrix(BB.uni)
  # Bivariate spline basis for all the points
  B0 <- B0.uni[index.ss, ]
  # Q2
  Q2 <- B.uni$Q2
  BQ2.uni <- B0.uni %*% Q2
  BQ2 <- BQ2.uni[index.ss, ]

  # generate univariate spline basis
  tt.uni <- unique(tt)
  index.tt <- match(tt, tt.uni)
  U0.uni <- bSpline(tt.uni, knots = time.knots, intercept = TRUE, degree = rho,
                    time.bound = time.bound)
  U0 <- U0.uni[index.tt, ]

  U0 <- Matrix(U0, sparse = TRUE)
  B0 <- Matrix(B0, sparse = TRUE)

  # tensor inner product
  if (line.up == TRUE) {
    # need to be checked
    Psi <- Matrix::kronecker(U0, B0)
    Psi.Q2 <- Matrix::kronecker(U0, BQ2)
  } else {
    # cat(class(U0))
    # cat(class(B0))
    Psi <- Matrix::t(Matrix::KhatriRao(Matrix::t(U0), Matrix::t(B0)))
    Psi.Q2 <- Matrix::t(Matrix::KhatriRao(Matrix::t(U0), Matrix::t(BQ2)))

    # Psi <- kr(U0, B0)
    # Psi.Q2 <- kr(U0, BQ2)
  }

  # generate Q2.all and H.all
  dimU <- ncol(U0)
  Q2.all <- kronecker(diag(dimU), Q2)
  H <- B.uni$H
  H.all <- kronecker(diag(dimU), H)

  Energ <- energy.tensor(V, Tri, d, time.knots, degr = rho, time.bound)
  P1 <- Energ$Eng1
  P2 <- Energ$Eng2
  D1 <- crossprod(Q2.all, as.matrix(P1)) %*% Q2.all
  D2 <- crossprod(Q2.all, as.matrix(P2)) %*% Q2.all


  list(
    Psi = Psi, Psi.Q2 = Psi.Q2, Q2 = Q2, H = H, Q2.all = Q2.all, H.all = H.all,
    dimB = ncol(B0), dimU = dimU, P1 = P1, P2 = P2, D1 = D1, D2 = D2
  )
}
