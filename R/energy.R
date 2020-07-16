energy <- function(V, Tr, d) {
  n <- dim(Tr)[1]
  m <- (d + 1) * (d + 2) / 2
  msqr <- m * m
  Mat <- BPST::build(d - 2)
  Indx1 <- matrix(rep(0, n * msqr), ncol = 1)
  Indx2 <- Indx1
  S <- Indx1
  place <- 1
  for (k in 1:n) {
    LocK <- locEng(V[Tr[k, 1], ], V[Tr[k, 2], ], V[Tr[k, 3], ], Mat, d)
    TmpX <- as(LocK, "dgTMatrix")
    i <- TmpX@i + 1
    j <- TmpX@j + 1
    s <- TmpX@x
    L <- length(i)
    Indx1[place:(place + L - 1)] <- (k - 1) * m + i
    Indx2[place:(place + L - 1)] <- (k - 1) * m + j
    S[place:(place + L - 1)] <- s
    place <- place + L
  }
  K <- sparseMatrix(Indx1[1:(place - 1)], Indx2[1:(place - 1)],
                    x = S[1:(place - 1)],
                    dims = c(n * m, n * m)
  )

  Indx1 <- matrix(rep(0, n * msqr), ncol = 1)
  Indx2 <- Indx1
  S <- Indx1
  place <- 1
  for (k in 1:n) {
    area <- diag(rep(abs(triarea(V[Tr[k, 1], ], V[Tr[k, 2], ], V[Tr[k, 3], ])), m), m)
    LocK <- area %*% as.matrix(BPST::build(d))
    TmpX <- as(LocK, "dgTMatrix")
    i <- TmpX@i + 1
    j <- TmpX@j + 1
    s <- TmpX@x
    L <- length(i)
    Indx1[place:(place + L - 1)] <- (k - 1) * m + i
    Indx2[place:(place + L - 1)] <- (k - 1) * m + j
    S[place:(place + L - 1)] <- s
    place <- place + L
  }
  K1 <- sparseMatrix(Indx1[1:(place - 1)], Indx2[1:(place - 1)],
                     x = S[1:(place - 1)],
                     dims = c(n * m, n * m)
  )
  return(list(K2deri = K, K2 = K1))
}
