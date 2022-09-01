######################################################
#' Incidence Matrix, Information Matrix, Canonical efficiency factor, Variance between associates and average variance
#'
#' @param design p-rep design in matrix form considering rows as blocks
#'
#'@description This function generates incidence matrix, information matrix, canonical efficiency factor,variance factor between associates and average variance for the input design
#' @return This function generates incidence matrix, information matrix, canonical efficiency factor, variance factor between associates and average variance factor for the input design
#' @export
#'
#' @examples
#'library(pRepDesigns)
#'design=matrix(1:9, nrow = 3, ncol = 3)
#'NCEV(design)
#'@references
#'Williams E, Piepho HP, Whitaker D (2011)<https://doi.org/10.1002/bimj.201000102>
######################################################
NCEV=function(design){
  v = max(design)
  b = nrow(design)
  k = ncol(design)
  N = matrix(0, v, b)
N_matrix=function(design)
{

  for (i in 1:b) {
    for (j in 1:k) {
      N[design[i, j], i] = N[design[i, j], i] + 1
    }
  }
  N
}
N_mat=N_matrix(design)
###########################################################
C_matrix=function(N_mat){
  v1=nrow(N_mat)
  b1=ncol(N_mat)
  r=3
  K=diag(colSums(N_mat), b1, b1)
  R=diag(rowSums(N_mat), v1, v1)
  kvec=colSums(N_mat)
  Kinv = diag(1/kvec, nrow = b1, ncol = b1)
  C = R - N_mat %*% Kinv %*% t(N_mat)
  C
}
C_mat<-C_matrix(N_mat)
##########################################################
C_Efficiency=function(C_mat){
  E=eigen(C_mat, only.values = T)
  r=3
  E1=unlist(E)
  E_positive=E1[E1>=0.000000001]
  n=length(E_positive)
  C_Efficiency=n/(r*sum(c(1/E_positive)))
  C_Efficiency
}
C_E=C_Efficiency(C_mat)
##############################################################
p_matrix=matrix(,nrow=0,ncol=v)
i=1
j=1
while(i<=(choose(v,2))){
  j=i+1
  while(j<=v){
    p1<-matrix(0,nrow=1,ncol=v)
    p1[i]<-1
    p1[j]<--1
    p_matrix<-rbind(p_matrix,p1)
    j=j+1
  }
  i=i+1
}
p_matrix
p_invC_Pprme=(p_matrix)%*%MASS::ginv(C_mat)%*%t(p_matrix)
var<-diag(p_invC_Pprme)
var1<-round(var,digits=4)
var2<-unique(var1)
Average_var<-mean(var)
###########################################################
results=list( "Incidence matrix"=N_mat, "C-matrix"=round(C_mat, 4), "Canonical Efficiency factor"= round(C_E,4), "Variance factor b/w associates"=var2,"Average variance factor"=round(Average_var,4))
print(results)
}
##################################################

