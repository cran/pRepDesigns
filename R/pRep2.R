#p-rep 1 using graphs
#######################3######################
#' Second series of p-rep designs
#'
#' @param v Total number of treatments or breeding lines or entries
#' @param p positive integer (>=2)
#'
#' @return
#' This function generates a new series of p-rep designs with parameters:
#'
#' v = 10p   number of treatments,
#'
#' e = 2     number of environments,
#'
#' b = 10    blocks of size k = 3p and
#'
#' r = 3     number of replications.
#'
#' This function also generates canonical efficiency factor and average variance factor of the generated p-rep design.
#' @export
#' @description
#' For the specified values of v and p, this function generates the second series of p-rep designs. The input should meet the condition that v=10*p where, p>=2.
#'
#' @examples
#' library(pRepDesigns)
#' pRep2(20,2)
#' @references
#'Williams E, Piepho HP, Whitaker D (2011)<https://doi.org/10.1002/bimj.201000102>

#############################################
pRep2=function(v,p){
  if(p>=2 && v==10*p){
    v=10*p
    b=10
    r=3
    k=3*p
    e=2
    vv=1
    z=c()
    while(vv<=v){
      x=c(vv)
      z=c(z, x)
      vv=vv+1
    }
    z=matrix(z, nrow =10, ncol=p, byrow= T)
    block1=c(z[c(1,6,10), 1:p])
    block2=c(z[c(2,6,7), 1:p])
    block3=c(z[c(3,7,8), 1:p])
    block4=c(z[c(4,8,9), 1:p])
    block5=c(z[c(5,9,10), 1:p])
    block6=c(z[c(3,5,6), 1:p])
    block7=c(z[c(1,4,7), 1:p])
    block8=c(z[c(2,5,8), 1:p])
    block9=c(z[c(1,3,9), 1:p])
    block10=c(z[c(2,4,10), 1:p])
    e1=rbind(block1,block2, block3,block4, block5)
    e2=rbind(block6, block7, block8, block9, block10)
    result=list("Environment 1"=e1,"Environment 2"=e2)
    design=rbind(e1,e2)
######################################################
    N_matrix=function(design)
    {
      v = max(design)
      b = nrow(design)
      k = ncol(design)
      N = matrix(0, v, b)
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
      v=nrow(N_mat)
      b=ncol(N_mat)
      r=3
      K=diag(colSums(N_mat), b, b)
      R=diag(rowSums(N_mat), v, v)
      kvec=colSums(N_mat)
      Kinv = diag(1/kvec, nrow = b, ncol = b)
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
    A1=c("Number of treatments","Number of blocks","Number of replications","Block size ","Number of environments")
    A2=c("v", "b", "r", "k", "e")
    A3=c(v, b, r, k, e)
    A=cbind(A1, A2, A3)
    prmatrix(A,rowlab=,collab=rep("",ncol(A)),quote=FALSE,na.print="")
    message("\n")
#######################################################################
    message("p-rep design")
    message("\n","Environment_1 of p-rep design")
    prmatrix(e1,rowlab=,collab=rep("",ncol(e1)),quote=FALSE,na.print="")
    message(c("\n","Environment_2 of p-rep design"))
    prmatrix(e2,rowlab=,collab=rep("",ncol(e2)),quote=FALSE,na.print="")
    message("\n")
########################################################################
    print(c("Canonical Efficiency factor", round(C_E,4)),quote=F)
##########################################################################
    B1=c("variance factor between first associates","variance factor between second associates","variance factor between third associates")
    B2=c(var2[1],var2[2],var2[3])
    B=cbind(B1,B2)
    prmatrix(B,rowlab=,collab=rep("",ncol(B)),quote=FALSE,na.print="")
    message("\n")
#########################################################################
    print(c("Average variance factor", round(Average_var,4)),quote=F)
#########################################################################
  } else {
    message("Please enter v(=10*p, where p>=2)")
  }
}
##############################################################
#pRep2(20,2)


