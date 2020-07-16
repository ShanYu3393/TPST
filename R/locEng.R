#' Calculate local energy matrix
#'
#' This function computes the local
#' energy matrix with respect to the triangle \code{(V1,V2,V3)}.
#' @param V1,V2,V3
#' three vectors of length two indicating the vertices locations of a triangle.
#' @param Mat inner coefficients matrix. This arguments could be
#' calculated by function \code{build}.
#' @param d degree of polynomials.
#' @return energy matrix with respect to the triangle \code{(V1,V2,V3)}.
#' @export
#' This function computes the local
#' energy matrix with respect to the triangle \code{(V1,V2,V3)}.
#' @param V1,V2,V3
#' three vectors of length two indicating the vertices locations of a triangle.
#' @param Mat inner coeffients matrix. This arguments could be
#' calculated by function \code{build}.
#' @param d degree of polynomials.
#' @return energy matrix with respect to the triangle \code{(V1,V2,V3)}.
locEng <- function(V1, V2, V3, Mat, d) {
  Mat <- as.matrix(Mat)
  m <- (d + 1) * (d + 2) / 2
  Id <- diag(m)
  vx <- c(1, 0)
  vy <- c(0, 1)
  lamx <- tcord(V1, V2, V3, vx)
  lamy <- tcord(V1, V2, V3, vy)
  Dx <- dirder(Id, lamx[1], lamx[2], lamx[3])
  Dxx <- dirder(Dx, lamx[1], lamx[2], lamx[3])
  Dxy <- dirder(Dx, lamy[1], lamy[2], lamy[3])
  Dy <- dirder(Id, lamy[1], lamy[2], lamy[3])
  Dyy <- dirder(Dy, lamy[1], lamy[2], lamy[3])

  area <- diag(rep(abs(triarea(V1, V2, V3)), m), m)
  K <- area %*% (crossprod(Dxx, Mat) %*% Dxx + crossprod(Dyy, Mat) %*% Dyy)
  return(K)
}

#source('triarea.R')
locEng=function(V1,V2,V3,Mat,d){
	Mat=as.matrix(Mat)
	m=(d+1)*(d+2)/2
	Id=diag(m)
	vx=c(1,0)
	vy=c(0,1)
	lamx=tcord(V1,V2,V3,vx)
	lamy=tcord(V1,V2,V3,vy)
	Dx=dirder(Id,lamx[1],lamx[2],lamx[3])
	Dxx=dirder(Dx,lamx[1],lamx[2],lamx[3])
	Dxy=dirder(Dx,lamy[1],lamy[2],lamy[3])
	Dy=dirder(Id,lamy[1],lamy[2],lamy[3])
	Dyy=dirder(Dy,lamy[1],lamy[2],lamy[3])
	# K=abs(triarea(V1,V2,V3))*(t(Dxx)%*%Mat%*%Dxx+
		# 2*t(Dxy)%*%Mat%*%Dxy+t(Dyy)%*%Mat%*%Dyy)
	area=diag(rep(abs(triarea(V1,V2,V3)),m),m)
	K=area%*%(crossprod(Dxx,Mat)%*%Dxx+crossprod(Dyy,Mat)%*%Dyy)
	return(K)
}
