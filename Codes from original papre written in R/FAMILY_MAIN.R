function (X, Z, Y = numeric(), lambdas = c(), alphas = c(), family = c("gaussian", 
    "binomial"), rho = 1, B = NULL, norm = "l2", quad = TRUE, 
    iter = 500, e.abs = 0.001, e.rel = 0.001, maxiter.B = 50, 
    tol.B = 1e-04, verbose = FALSE) 
{
    start.time <- proc.time()[3]
    if (!is.numeric(X) | !is.numeric(Z) | !is.numeric(Y)) {
        stop("Data must be a numeric matrix/vector")
    }
    if (!is.vector(Y)) {
        stop("Response Y must be a vector")
    }
    if (any(lambdas < 0)) {
        stop("Lambda values must be positive")
    }
    if (any(alphas <= 0) | any(alphas >= 1)) {
        stop("Alpha values must be between 0 and 1")
    }
    if (nrow(X) != nrow(Z) | nrow(X) != length(Y) | nrow(Z) != 
        length(Y)) {
        stop("Mismatched dimensions of data. Make sure nrow(X)=nrow(Z)=length(Y)")
    }
    family = family[1]
    if (family == "binomial") {
        if (!all(sort(unique(Y)) == c(0, 1))) {
            stop("Response vector must be a binary vector")
        }
    }
    if (!(family == "gaussian" | family == "binomial")) {
        stop("Unrecognized family name. Family must be 'gaussian' or 'binomial' ")
    }
    if (is.null(B)) {
        B <- matrix(0, ncol = ncol(Z) + 1, nrow = ncol(X) + 1)
        if (quad == FALSE) {
            diag(B)[-1] = NaN
        }
    }
    if (quad == FALSE) {
        diag(B)[-1] <- NaN
    }
    matD = matE = matF = B
    g1 <- rep(0, ncol(B) * nrow(B))
    g2 <- rep(0, ncol(B) * nrow(B))
    g3 <- rep(0, ncol(B) * nrow(B))
    g.list <- list(matrix(g1, ncol = ncol(B)), matrix(g2, ncol = ncol(B)), 
        matrix(g3, ncol = ncol(B)))
    length.alpha <- length(alphas)
    length.lambda <- length(lambdas)
    fin.result <- vector("list", length(alphas))
    cat("Computing w...")
    if (quad) {
        w <- generate.w(X = X, Z = Z, quad)
    }
    else {
        w <- generate.w(X = X, Z = Z, quad)
        w.full <- generate.w(X = X, Z = Z, TRUE)
    }
    cat("done.\n")
    num <- length(Y)
    cat("Starting svd...")
    svd.w <- svd(w)
    svd.w$tu <- t(svd.w$u)
    svd.w$tv <- t(svd.w$v)
    cat("done.\n")
    p1 <- ncol(X)
    p2 <- ncol(Z)
    for (i in 1:length.alpha) {
        alpha <- alphas[i]
        b.array <- vector("list", length.lambda)
        if (verbose) 
            cat("Fitting model for alpha =", round(alpha, 2), 
                "and lambda =", round(lambdas[length.lambda], 
                  2), "\n")
        if (family == "gaussian") {
            b.array[[length.lambda]] <- estimate.SH(X, Z, Y, 
                w, svd.w = svd.w, c(lambdas[length.lambda] * 
                  (1 - alpha) * sqrt(p2), lambdas[length.lambda] * 
                  (1 - alpha) * sqrt(p1), alpha * lambdas[length.lambda]), 
                rho = rho, B, matD, matE, matF, g.list, iter = iter, 
                e.abs = e.abs, e.rel = e.rel, quad = quad, norm = norm)
            b.array[[length.lambda]]$alpha = alpha
            b.array[[length.lambda]]$lambda = lambdas[length.lambda]
        }
        else {
            b.array[[length.lambda]] <- estimate.logistic(X, 
                Z, Y, w, svd.w = svd.w, c(lambdas[length.lambda] * 
                  (1 - alpha) * sqrt(p2), lambdas[length.lambda] * 
                  (1 - alpha) * sqrt(p1), alpha * lambdas[length.lambda]), 
                rho = rho, B, matD, matE, matF, g.list, iter = iter, 
                e.abs = e.abs, e.rel = e.rel, quad = quad, norm = norm)
            b.array[[length.lambda]]$alpha = alpha
            b.array[[length.lambda]]$lambda = lambdas[length.lambda]
        }
        if (!b.array[[length.lambda]]$conv) 
            warning(paste("The algorithm did not converge for alpha= ", 
                alpha, " and lambda= ", lambdas[length.lambda], 
                "did not converge"))
        for (lam in (length.lambda - 1):1) {
            if (verbose) 
                cat("Fitting model for alpha =", round(alpha, 
                  2), "and lambda =", round(lambdas[lam], 2), 
                  "\n")
            if (family == "gaussian") {
                b.array[[lam]] <- estimate.SH(X, Z, Y, w, svd.w = svd.w, 
                  c(lambdas[lam] * (1 - alpha) * sqrt(p2), lambdas[lam] * 
                    (1 - alpha) * sqrt(p1), alpha * lambdas[lam]), 
                  rho = b.array[[lam + 1]]$rho, b.array[[lam + 
                    1]]$B, b.array[[lam + 1]]$D, b.array[[lam + 
                    1]]$E, b.array[[lam + 1]]$F, b.array[[lam + 
                    1]]$glist, iter = iter, e.abs = e.abs, e.rel = e.rel, 
                  quad = quad, norm = norm)
                b.array[[lam]]$alpha = alpha
                b.array[[lam]]$lambda = lambdas[lam]
            }
            else {
                b.array[[lam]] <- estimate.logistic(X, Z, Y, 
                  w, svd.w = svd.w, c(lambdas[lam] * (1 - alpha) * 
                    sqrt(p2), lambdas[lam] * (1 - alpha) * sqrt(p1), 
                    alpha * lambdas[lam]), rho = b.array[[lam + 
                    1]]$rho, b.array[[lam + 1]]$B, b.array[[lam + 
                    1]]$D, b.array[[lam + 1]]$E, b.array[[lam + 
                    1]]$F, b.array[[lam + 1]]$glist, iter = iter, 
                  e.abs = e.abs, e.rel = e.rel, quad = quad, 
                  norm = norm)
                b.array[[lam]]$alpha = alpha
                b.array[[lam]]$lambda = lambdas[lam]
            }
            if (!b.array[[lam]]$conv) 
                warning(paste("The algorithm did not converge for alpha= ", 
                  alpha, " and lambda= ", lambdas[lam], "did not converge"))
        }
        fin.result[[i]] <- b.array
    }
    end.time <- proc.time()[3]
    run.time <- end.time - start.time
    output.obj <- list(Estimate = fin.result, alpha = alphas, 
        lambda = lambdas, Y.train = Y, X.train = X, Z.train = Z, 
        family = family, quad = quad, call = match.call(), time = run.time, 
        w = w)
    class(output.obj) <- "FAMILY"
    return(output.obj)
}

