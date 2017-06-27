####################################################################

bc = function(node,N,p)
{
  res = rep(NA,N)
  res[seq(node+1,N)] = sapply(X = p[seq(1,N-node)],FUN = rbinom,n=1,size = 1)
  return(res)
}

####################################################################

expected = function(p)
{
  N = length(p)
  PS = ET = ER = EB = matrix(data = NA,nrow = N,ncol = 1)
  PS[1] = ER[1] = p[1]
  ET[1] = EB[1] = 1
  index = seq(1,N)
  for(i in index[index!=1])
  {
    PS[i] = 1 - (1 - p[i])*prod(1 - (p[1:i-1]*rev(PS[1:(i-1)])))
    ET[i] = i + sum(p[1:i-1]*rev(ET[1:(i-1)]))
    ER[i] = sum(p[1:i]) + sum(p[1:i-1]*rev(ER[1:(i-1)]))
    EB[i] = 1 + sum(p[1:i-1]*rev(EB[1:(i-1)]))
  }
  res = data.frame(round(t(cbind(PS,ET,ER,EB)),4),row.names = c("Success Probability","Exp Transmissions","Exp Receptions","Exp Broadcast"))
  colnames(res) = seq(1,N)
  return(res)
}

####################################################################

# rmat = function(N)
# {
#   w <- 1:N
#   n <- length(w)
#   t <- N
#   D <- list()
#   for (j in 0:n) D[[paste(0, j)]] <- list(c())
#   for (i in 1:t) D[[paste(i, 0)]] <- list()
#   for (j in 1:n) {
#     for (i in 1:t){
#       D[[paste(i, j)]] <- do.call(c, lapply(0:floor(i/w[j]), function(r) {
#         lapply(D[[paste(i-r*w[j], j-1)]], function(x) c(x, rep(w[j], r)))
#       }))
#     }
#   }
#   test = D[[paste(t, n)]]
#   return(matrix(data = unlist(lapply(lapply(test, factor,levels=1:N),table)),ncol = N,byrow = TRUE))
# }

####################################################################

MC = function(p,M)
{
  nbv = ntv = nrv = ndv = integer(M)
  pb <- txtProgressBar(min = 0, max = M, style = 3)
  for(i in 1:M)
  {
    model = oportunist(p)
    nbv[i] = model$nb
    ntv[i] = model$nt
    nrv[i] = model$nr
    ndv[i] = model$dv
    setTxtProgressBar(pb, i)
  }
  close(pb)
  mat = round(c(mean(ndv),mean(ntv),mean(nrv),mean(nbv)),5)
  matdf = data.frame(mat,row.names = c("Probability System","Exp Transmissions","Exp Receptions","Exp Broadcast"))
  colnames(matdf) = c("Monte Carlo")
  return(matdf)
}

####################################################################

multi = function(vec){
  return(factorial(sum(vec))/prod(sapply(X = vec,FUN = function(x){factorial(x)})))
}

####################################################################

oportunist = function(p){
  N = length(p)
  m1 = bc(0,N,p)
  M = matrix(m1,nrow = 1,ncol = N)
  pos = seq(1,N-1)[m1[-N]==1]
  while(length(pos[!is.na(pos)]) > 0)
  {
    pos = pos[!is.na(pos)]
    m2 = t(sapply(X = pos,FUN = bc,N=N,p=p))
    M = rbind(M,m2)
    pos = unlist(apply(X = m2,MARGIN = 1,function(x){seq(1,N-1)[x[-N]==1]}))
  }
  nt = length(c(M)[!is.na(c(M))])
  nr = sum(M, na.rm=TRUE)
  nb = dim(M)[1]
  dv = as.integer(any(M[,N]==1))
  return(list(M = M,nt = nt,nr = nr,nb = nb,dv=dv))
}

####################################################################

sols = function(N){
  s = 1:N
  vals = as.matrix(expand.grid(rep(list(0:N),N)))
  ii = 1:dim(vals)[1]
  ii = apply(X = vals,MARGIN = 1,FUN = function(x){ifelse(sum(x) < N,1,0)})
  ab   = vals[N+1,]
  vals = vals[ii==1,]
  vals = rbind(vals[1:N,],ab,vals[-(1:N),])
  jj = apply(X = vals,MARGIN = 1,FUN = function(x){ifelse(t(s)%*%x == N,1,0)})
  return(vals[jj==1,])}

####################################################################

gentext = function(row,pt)
{
  return(gsub("\\^1", "",paste(pt[row!=0],"^",row[row!=0],sep = "",collapse = "*"), perl=TRUE))
}

####################################################################

routes = function(p,delta = 0)
{
  if(length(p)<1) stop("p length must greater or equal than 1")
  if(max(p) > 1 | min(p) < 0) stop("p vector must contain real numbers in [0,1]")
  if(delta ==0){
    N = length(p)
    mat = sols(N)
    freq = cbind(apply(X = mat,MARGIN = 1,FUN = multi))
    NR = sum(freq)
    pt = paste("p",seq_len(N),sep = "")
    probs1 = apply(X = mat,MARGIN = 1,FUN = function(x){prod(p^x)})
    probs2 = apply(X = mat,MARGIN = 1,FUN = gentext,pt=pt)
    res = data.frame(rbind(cbind(freq,probs2,round(probs1,5)),c(NR,"","")),row.names = c(paste("route",seq_len(length(freq)),"  ",sep = " "),"Total"))
    colnames(res) = c("Freq","Probability","Value")
    return(res)
  }else{
    p0 = p1 = rep(0,length(p))
    p0[p!=0] = p[p!=0]-delta
    p1[p!=0] = p[p!=0]+delta
    if(max(p1) > 1 | min(p0) < 0) stop("lower and upper bounds must contain real numbers in [0,1]")
    N = length(p)
    mat = sols(N)
    freq = cbind(apply(X = mat,MARGIN = 1,FUN = multi))
    NR = sum(freq)
    pt = paste("p",seq_len(N),sep = "")
    probs1 = apply(X = mat,MARGIN = 1,FUN = function(x){prod(p^x)})
    probs10 = apply(X = mat,MARGIN = 1,FUN = function(x){prod(p0^x)})
    probs11 = apply(X = mat,MARGIN = 1,FUN = function(x){prod(p1^x)})
    probs2 = apply(X = mat,MARGIN = 1,FUN = gentext,pt=pt)
    res = data.frame(rbind(cbind(freq,probs2,round(probs10,5),round(probs1,5),round(probs11,5)),c(NR,"","","","")),row.names = c(paste("route",seq_len(length(freq)),"  ",sep = " "),"Total"))
    colnames(res) = c("Freq","Probability","p - delta","    p    ","p + delta")
    return(res)
  }
}
