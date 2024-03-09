#' @title mapping
#' @description Computes the mappings of several radial axes methods
#' @param algorithm defines the radial method (Posible values: SC, RadViz, Adaptable, Adaptable exact, Adaptable ordered
#' @param X is an N by n matrix wose rows contain the n dimensional data samples.
#' @param V is an n by m matrix whose rows define the method's axis vectors.
#' @param W is an n by n diagonal matrix defining nonnegative weights for each variable
#' @param vector_norm is the vector norm associated with adaptable radial axes plots
#' @param chosen_variable is the selected attribute for constrained adaptable radial axes plots
#' @return The low-dimensional embeddings are stored in the N by m matrix P.
#' @export mapping
#' @examples
#' N = 100
#' n = 5
#' m = 2
#' X <- matrix(rnorm(N, n), nrow = 20, ncol = n)
#' V <- matrix(rnorm(N, m), nrow = n, ncol = m)
#' W <- diag(1, nrow = n)
#' mapping('SC', X, V)
#' mapping('RadViz', X, V)
#' mapping('Adaptable', X, V, W, 1)
#' mapping('Adaptable', X, V, W, 2)
#' mapping('Adaptable', X, V, W, 'Inf')
#' mapping('Adaptable exact', X, V, W, 1, 1)
#' mapping('Adaptable exact', X, V, W, 2, 1)
#' mapping('Adaptable exact', X, V, W, 'Inf', 1)
#' mapping('Adaptable ordered', X, V, W, 1, 1)
#' mapping('Adaptable ordered', X, V, W, 2, 1)
#' mapping('Adaptable ordered', X, V, W, 'Inf', 1)

mapping <- function (algorithm, X, V, W, vector_norm, chosen_variable) {

  N <- nrow(X)
  n <- ncol(X)
  m <- ncol(V)

  if (algorithm == 'SC'){
    P <- X %*% V

  } else if (algorithm == 'RadViz'){
    minimum <- apply(X, 2, min)
    maximum <- apply(X, 2, max)
    X_Radviz <- (X - repmat(minimum, N, 1)) / (repmat(maximum, N, 1) - repmat(minimum, N, 1))
    sum_row <- rowSums(X_Radviz)
    for (i in 1:N){
      if (sum_row[i] == 0){
        X_Radviz[i,] <- ones(1, n) / n
      } else {
        X_Radviz[i,] <- X_Radviz[i,] / sum_row[i]
      }
    }
    P <- X_Radviz %*% V

  } else if (algorithm == 'Adaptable') {
    nrow_A <- nrow(V) * 2
    if (vector_norm == 1) {
      A <- rbind(cbind(-diag(1, n), -W %*% V), cbind(-diag(1, n), W %*% V))
      f <- c(cbind(matrix(1, 1, n), matrix(0, 1, m)))
      signos <- c(rep('<=', nrow_A))
      P <- zeros(N, m)
      bounds <- list(lower = list(ind = 1:length(f), val=c(rep(-Inf, length(f))),
                                  upper = list(ind = 1:length(f), val=c(rep(Inf, length(f))))))
      for (i in 1:N) {
        b <- c(rbind(-W %*% X[i,], W %*% X[i,]))
        p_star <- Rsymphony_solve_LP(f, A, signos, b, bounds)
        P[i,] <- t(tail(p_star$solution, m))
      }

    } else if (vector_norm == 2) {
        P <- X %*% t(pinv(W %*% V) %*% W)

    } else if (vector_norm == 'Inf') {
      A <- rbind(cbind(-ones(n, 1), -W %*% V), cbind(-ones(n, 1), W %*% V))
      f <- c(1, zeros(1, m))
      signos <- c(rep('<=', nrow_A))
      P <- zeros(N, m)
      bounds <- list(lower = list(ind = 1:length(f), val=c(rep(-Inf, length(f))),
                                  upper = list(ind = 1:length(f), val=c(rep(Inf, length(f))))))
      for (i in 1:N){
        b <- rbind(-W %*% X[i,], W %*% X[i,])
        p_star <- Rsymphony_solve_LP(f, A, signos, b, bounds)
        P[i,] <- t(tail(p_star$solution, m))
      }
    }
  } else if (algorithm == 'Adaptable exact') {
    v <- V[chosen_variable,]
    x <- X[,chosen_variable]
    nrow_A <- nrow(V) * 2
    if (vector_norm == 1){
      A <- rbind(cbind(-diag(n), -V), cbind(-diag(n), V))
      f <- c(ones(1, n), 0, 0)
      Aeq <- t(matrix(c(zeros(1, n), v)))
      A <- rbind(A, Aeq)
      signos <- c(rep('<=', nrow_A), '==')
      bounds <- list(lower = list(ind = 1:length(f+2), val=c(rep(-Inf, length(f+2))),
                                  upper = list(ind = 1:length(f+2), val=c(rep(Inf, length(f+2))))))
      P <- zeros(N, m)
      for (i in 1:N){
        b <- c(-X[i,], X[i,], x[i])
        p_star <- Rsymphony_solve_LP(f, A, signos, b, bounds)
        P[i,] <- t(tail(p_star$solution, m))
      }

    } else if (vector_norm ==  2) {
      v <- V[chosen_variable,]
      x <- X[,chosen_variable]
      P <- matrix(0, N, m)
      for (i in 1:N) {
        p <- lsqlincon(V, X[i,], NULL, NULL, v, x[i], NULL, NULL)
        P[i,] = t(p)
      }

    } else if (vector_norm == 'Inf'){
      A <- rbind(cbind(-ones(n, 1), -V), cbind(-ones(n, 1), V))
      f <- c(1, zeros(1, 2))
      Aeq <- t(matrix(c(0, v)))
      A <- rbind(A, Aeq)
      signos <- c(rep('<=', nrow_A), '==')
      bounds <- list(lower = list(ind = 1:length(f+2), val=c(rep(-Inf, length(f+2))),
                                  upper = list(ind = 1:length(f+2), val=c(rep(Inf, length(f+2))))))
      P <- zeros(N, m)
      for (i in 1:N) {
        b <- c(-X[i,], X[i,], x[i])
        p_star <- Rsymphony_solve_LP(f, A, signos, b, bounds)
        P[i,] <- t(tail(p_star$solution, m))
      }
    }
  } else if (algorithm == 'Adaptable ordered') {
    k <- chosen_variable
    I <- order(c(X[,k]))
    constraints <- list()
    if (vector_norm == 1){
      P = Variable(N, m)
      z <- 0
      for (i in 2:N) {
        z = z + cvxr_norm(V %*% t(P[i,]) - X[i,], 1)
      }
      obj <- Minimize(z)
      for (i in 2:N-1) {
        constraints <- list.append(constraints, P[I[i],] %*% V[k,] <= P[I[i+1],] %*% V[k,])
      }
      prob <- Problem(obj, constraints)
      solution <- solve(prob)
      P <- solution$getValue(P)

    } else if (vector_norm == 2) {
      P <- Variable(N, m)
      obj <- Minimize(cvxr_norm(P %*% t(V) - X, 'fro'))
      for (i in 2:N-1) {
        constraints <- list.append(constraints, P[I[i],] %*% V[k,] <= P[I[i+1],] %*% V[k,])
      }
      prob <- Problem(obj, constraints)
      solution <- solve(prob)
      P <- solution$getValue(P)

    } else if (vector_norm == 'Inf') {
      P <- Variable(N, m)
      x <- Variable()
      obj <- Minimize(x)
      for (j in 1:n) {
        z <- cvxr_norm(P %*% V[j,] - X[,j], "inf")
        constraints <- list.append(constraints, z <= x)
      }
      for (i in 2:N-1) {
        constraints <- list.append(constraints, P[I[i],] %*% V[k,] <= P[I[i+1],] %*% V[k,])
      }
      prob <- Problem(obj, constraints)
      solution <- solve(prob)
      P <- solution$getValue(P)
    }
  }
  return(P)
}
