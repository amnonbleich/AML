# Variable importance for pls, implementation taken from
# http://mevik.net/work/software/VIP.R

## VIPjh returns the VIP of variable j with h components
VIPjh <- function(object, j, h) {
    if (object$method != "oscorespls")
        stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
    if (nrow(object$Yloadings) > 1)
        stop("Only implemented for single-response models")

    b <- c(object$Yloadings)[1:h]
    # coefficients of the model Y=Q1*PC1+Q2*PC2... with h components
    T <- object$scores[,1:h, drop = FALSE]
    # projections of the points (genes) on the directions of the first h pls
    # components
    # T matrix: how much each gene contributed to each component.
    SS <- b^2 * colSums(T^2)
    # sum the scores of the genes for each individual component and multiply by
    # the coefficient (variance of Y)
    W <- object$loading.weights[,1:h, drop = FALSE]
    # contribution (weight) of each gene to each of the 10 components
    Wnorm2 <- colSums(W^2)
    # total contribution for the ten component, normalization factor for the
    # weight
    v <- sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
    # colname(v) = rownames(W)[j]
    # formula for the lecture.
}
