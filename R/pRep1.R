###################################
#' First series of p-rep designs
#'
#' @param v Total number of treatments or breeding lines or entries
#' @param m positive integer (>=1)
#' @param s positive integer (>=3)
#'@description This function generates first series of p-rep designs for given values of v, m and s. The input should satisfy the condition v=2*m*s*(s-1), m>=1 and s>=3.
#' @return
#' This function generates p-rep designs with parameters:
#'
#' v = 2ms(s-1) : number of treatments,
#'
#' b_1 = 2(s-1) : first set of blocks of size k_1 = ms,
#'
#' b_2 = 2s : second set of blocks of size k_2 = 2m(s-1) and
#'
#' r = 3 : number of replications.
#'
#' e = 2 : number of environments,
#'
#' This function also generates incidence matrix, information matrix, canonical efficiency factor and average variance factor of the generated p-rep design.
#' @export
#'
#' @examples
#' library(pRepDesigns)
#' pRep1(24, 2, 3)
#'
#'@references
#'Williams E, Piepho HP, Whitaker D (2011)<doi:10.1002/bimj.201000102>
########################################################
################################################
pRep1<-function(v,m,s){
  if(m>=1 && s>=3 && v==2*m*s*(s-1)){
    v=2*m*s*(s-1)
    b1=2*(s-1)
    b2=2*s
    r=3
    k1=m*s
    k2=2*m*(s-1)
    e=2
    vv=1
    y=c()
    while(vv<=(v)){
      x=c(vv)
      y=c(y, x)
      vv=vv+1
    }
    y=matrix(y, nrow = s, byrow= F, ncol = 2*m*(s-1))
    y
    ######################################
    z=c()
    while(vv<=v){
      x=c(vv)
      z=c(z, x)
      vv=vv+1
    }
    z=matrix(y, nrow = 2*(s-1), byrow= T, ncol = m*s)
    z
    #####################################
    h<-nrow(z)/2
    cz<-ncol(z)
    xx<-z[1:h,]
    blank<-matrix(,nrow=nrow(xx),ncol=(ncol(y)-ncol(xx)))
    a<-cbind(xx,blank)
    final1<-rbind(a,y)
    row.names(final1)<-c(1:nrow(final1))
    ss=1
    while(ss<=nrow(final1)){
      rownames(final1)[ss] <- paste("block",as.character(ss),sep="")
      ss=ss+1
    }
    ########################################
    xy<-z[(1+nrow(z)/2):nrow(z),]
    blank<-matrix(,nrow=nrow(xy),ncol=(ncol(y)-ncol(xy)))
    b<-cbind(xy,blank)
    final2<-rbind(b,y)
    row.names(final2)<-c(1:nrow(final2))
    ss=nrow(final1)+1
    sp=1
    while(ss<=(2*(nrow(final1)))){
      rownames(final2)[sp] <- paste("block",as.character(ss),sep="")
      sp=sp+1
      ss=ss+1
    }
    #######################################
    design<-rbind(final1,final2)
    design[is.na(design)]=0
    #######################################
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
      CE=n/(r*sum(c(1/E_positive)))
    }
    C_E=C_Efficiency(C_mat)
    ##############################################################
    nc=ncol(C_mat)
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
    var1<-round(var, digits=4)
    var2<-unique(var1)
    Average_var<-mean(var)
    ###########################################################
    A1=c("Number of treatments","First set of blocks","Second set of blocks","Number of replications","Block size of b1","Block size of b2","Number of environments")
    A2=c("v", "b1", "b2","r", "k1", "k2", "e")
    A3=c(v, b1, b2, r, k1, k2, e)
    A=cbind(A1, A2, A3)
    prmatrix(A,rowlab=,collab=rep("",ncol(A)),quote=FALSE,na.print="")
    message("\n")
    ##########################################################
    message("p-rep design")
    message("\n","Environment_1 of p-rep design")
    prmatrix(final1,rowlab=,collab=rep("",ncol(final1)),quote=FALSE,na.print="")
    message(c("\n","Environment_2 of p-rep design"))
    prmatrix(final2,rowlab=,collab=rep("",ncol(final2)),quote=FALSE,na.print="")
    message("\n")
    #############################################################
    print(c("Canonical Efficiency factor", round(C_E,4)),quote=F)
    #############################################################
    B1<-c("variance factor between first associates","variance factor between second associates","variance factor between third associates","variance factor between fourth associates")
    B2<-c(var2[1],var2[2],var2[3],var2[4])
    if(m==1){
      B1<-B1[1:3]
      B2<-B2[1:3]
    }
    B<-cbind(B1,B2)
    prmatrix(B,rowlab=,collab=rep("",ncol(B)),quote=FALSE,na.print="")
    message("\n")
    ###############################################################
    print(c("Average variance factor", round(Average_var,4)),quote=F)
    #################################################################
  } else {
    message("Please enter v(=2*m*s*(s-1), where m>=2 and s>=3)")
  }
}
################################################################################
#pRep1(12, 1, 3)


























