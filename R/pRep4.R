#' p-rep designs with unequal block sizes
#'
#' @param v Total number of treatments
#' @param m Positive integer (>=1)
#' @param s Positive integer (>=3)
#'@description The first step of this function involves generating a higher
#'associate PBIB (partially balanced incomplete block) design consisting of two
#'sets of blocks. Subsequently, by creating various combinations of these two
#'sets of blocks, it generates multiple p-rep designs in different environments
#'with unequal block sizes.
#' @return This function calculates design parameters (v, b1, b2, r, k1, k2),
#' average variance factors, and canonical efficiency factors of generated
#' designs.
#' @export
#'
#' @examples
#' \dontrun{
#' library(pRepDesigns)
#' pRep4(48, 2, 4)}

pRep4=function (v, m, s)
{
  if (m >= 1 && s >= 3 && v == 2 * m * s * (s - 1)) {

    ###################
    v=2*m*s*(s-1)
    b1=2*(s-1)
    b2=2*s
    r=2
    k1=m*s
    k2=2*m*(s-1)

    vv = 1
    y = c()
    while (vv<=(v)){
      x=c(vv)
      y=c(y,x)
      vv=vv+1
    }
    mat1=matrix(y, nrow=s, byrow=F, ncol=2*m*(s-1))
    z1=c()
    while(vv<=v) {
      x = c(vv)
      z = c(z1, x)
      vv = vv + 1
    }
    z1=matrix(y, nrow = 2 * (s - 1), byrow = T, ncol = m * s)
    z2=matrix(0, nrow = 2 * (s - 1), ncol = k2-k1)
    mat2 = cbind(z1,z2)
    mat=rbind(mat1, mat2)

    row.names(mat)=c(1:nrow(mat))
    a1=1
    while(a1<=nrow(mat)){
      rownames(mat)[a1]=paste("Block", as.character(a1), sep = "")
      a1=a1+1
    }
  }
  else {
    message("Please enter v(=2*m*s*(s-1), where m>=1 and s>=3)")
  }



  ####################################Average var and CEF code
  CEAV=function(design){
    v=max(design)
    N_matrix = function(design) {
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
    N_mat = N_matrix(design)
    C_matrix = function(N_mat) {
      v = nrow(N_mat)
      b = ncol(N_mat)
      K = diag(colSums(N_mat), b, b)
      R = diag(rowSums(N_mat), v, v)
      kvec = colSums(N_mat)
      Kinv = diag(1/kvec, nrow = b, ncol = b)
      C = R - N_mat %*% Kinv %*% t(N_mat)
      C
    }
    C_mat <- C_matrix(N_mat)
    C_Efficiency = function(C_mat) {
      E = eigen(C_mat, only.values = T)
      r = length(which(design==max(design)))
      E1 = unlist(E)
      E_positive = E1[E1 >= 1e-09]
      n = length(E_positive)
      CE = n/(r * sum(c(1/E_positive)))
    }
    C_E = C_Efficiency(C_mat)
    nc = ncol(C_mat)
    p_matrix = matrix(, nrow = 0, ncol = v)
    i = 1
    j = 1
    while (i <= (choose(v, 2))) {
      j = i + 1
      while (j <= v) {
        p1 <- matrix(0, nrow = 1, ncol = v)
        p1[i] <- 1
        p1[j] <- -1
        p_matrix <- rbind(p_matrix, p1)
        j = j + 1
      }
      i = i + 1
    }
    p_matrix
    p_invC_Pprme = (p_matrix) %*% MASS::ginv(C_mat) %*% t(p_matrix)
    var <- diag(p_invC_Pprme)
    var1 <- round(var, digits = 4)
    var2 <- unique(var1)
    Average_var <- mean(var)
    listm=list("Average Variance Factor"=Average_var,"Cannonical Efficiency Factor"=C_E)
    return(listm)
  }

  message("__________ PBIB Design __________")
  cat("\n")
  print(mat)
  cat("\n")
  A11=c("Number of treatments (v)","First set of blocks (b1)",
        "Second set of blocks (b2)","Number of replications (r)",
        "Size of b1 blocks (k1)","Size of b2 blocks (k2)")
  A22=c(v, b1, b2, r, k1, k2)
  A1122 = cbind(A11, A22)
  print("Design parameters",quote=F)
  prmatrix(A1122, rowlab = , collab = rep("", ncol(A1122)), quote = FALSE, na.print = "")
  cat("\n")
  print(CEAV(mat))
  cat("\n")



  message("__________ p-Rep designs __________")
  cat("\n")

  #####################################
  if(nrow(mat1)%%2!=0){
    boolean1=F
    max1=1
  }else{
    boolean1=T
    max1=nrow(mat1)/2
  }
  if(nrow(mat2)%%2!=0){
    boolean2=F
    max2=1
  }else{
    boolean2=T
    max2=nrow(mat2)/2
  }
  ############### column gap part
  if(ncol(mat1)>ncol(mat2)){
    gap=ncol(mat1)-ncol(mat2)
    gap_matrix=matrix(0,nrow=nrow(mat2),ncol=gap)
    mat2=cbind(mat2,gap_matrix)
  }else{
    gap=ncol(mat2)-ncol(mat1)
    gap_matrix=matrix(0,nrow=nrow(mat1),ncol=gap)
    mat1=cbind(mat1,gap_matrix)
  }
  ##############for first matrix
  list1=list()
  for(i in 1:max1){
    if(i==1){
      for(j in 1:nrow(mat1)){
        list1=append(list1,list(mat1[j,]))
      }
    }
    if(i>1){
      n=nrow(mat1)/i
      for(j in 1:n){
        list1=append(list1,list(mat1[((i*j)-(i-1)):(i*j),]))

      }

    }
  }

  ###################for second matrix
  list2=list()
  for(i in 1:max2){
    if(i==1){
      for(j in 1:nrow(mat2)){
        list2=append(list2,list(mat2[j,]))
      }
    }
    if(i>1){
      n=nrow(mat2)/i
      for(j in 1:n){
        list2=append(list2,list(mat2[((i*j)-(i-1)):(i*j),]))
      }

    }
  }

  ######################part of first matrix with whole second matrix
  newlist=list()
  for(i in 1:length(list1)){
    final1=rbind(list1[[i]],mat2)
    blk_name=NULL
    for(i in 1:nrow(final1)){
      blk_name=c(blk_name,paste0("Block-",i))
    }
    row.names(final1)<-blk_name
    newlist=append(newlist,list(final1))
  }
  ############ indexing of matrix 1
  vector1=c()
  for(i in 1:max1){
    vector1=c(vector1,nrow(mat1)/i)
  }
  ###############
  finallist1=list()
  total=c(1:sum(vector1))

  for(i in 1:length(vector1)){
    seq=total[1:vector1[i]]
    total=setdiff(total,seq)
    finallist1=append(finallist1,list(newlist[seq]))
  }
  ##############part of second matrix with whole first matrix
  newlist1=list()
  for(i in 1:length(list2)){
    final2=rbind(list2[[i]],mat1)
    blk_name=NULL
    for(i in 1:nrow(final2)){
      blk_name=c(blk_name,paste0("Block-",i))
    }
    row.names(final2)<-blk_name
    newlist1=append(newlist1,list(final2))
  }
  #################
  ############ indexing of matrix 2
  vector2=c()
  for(i in 1:max2){
    vector2=c(vector2,nrow(mat2)/i)
  }
  ###############
  finallist2=list()
  total=c(1:sum(vector2))
  for(i in 1:length(vector2)){
    seq=total[1:vector2[i]]
    total=setdiff(total,seq)
    finallist2=append(finallist2,list(newlist1[seq]))
  }
  ####################
  finallist=append(finallist1,finallist2)
  ####################

  for(i in 1:length(finallist)){
    allenv=NULL
    message(paste("__________ p-Rep Design",i))

    cat("\n")
    for(j in 1:length(finallist[[i]])){
      print(paste0("Environment-",j),quote=F)
      blankmat=finallist[[i]][[j]]
      blankmat[blankmat==0]<-NA
      a=c()
      for(kk in 1:nrow(blankmat)){
        a=c(a,length(sort(blankmat[kk,])))
      }

      print(blankmat,na.print="",quote=F)
      cat("\n")
      allenv=rbind(allenv,finallist[[i]][[j]])
    }
    k1=min(a)
    k2=max(a)
    r=length(which(allenv==1))
    b11=length(which(a==k1))
    b22=length(which(a==k2))
    b1=length(finallist[[i]])*b11
    b2=length(finallist[[i]])*b22
    A1=c("Number of treatments (v)","First set of blocks (b1)",
         "Second set of blocks (b2)","Number of replications (r)",
         "Size of b1 blocks (k1)","Size of b2 blocks (k2)")
    A2=c(v, b1, b2, r, k1, k2)
    A = cbind(A1, A2)
    print("Design parameters",quote=F)
    prmatrix(A, rowlab = , collab = rep("", ncol(A)), quote = FALSE, na.print = "")
    cat("\n")
    print(CEAV(allenv))
  }
}
