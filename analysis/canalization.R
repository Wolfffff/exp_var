A = c(1, 1, 1, 1)
B = rep(0, 4)
P = c(0.5, 0.5, 0.5, 0.5) + c(0, 0, 1, 1)

Normalize = function(x) 
  return(x / sqrt(sum(x^2)))
}
Norm <- function(x){
    return(sqrt(sum(x^2)))
}
getDist = function(P, A, B = rep(0, length(A))){
    pa = P - A
    ba = B - A
    t = (pa %*% ba) / (ba %*% ba)
    d = Norm(pa - c(t) * ba)
    d
}
getDist(P, A)


