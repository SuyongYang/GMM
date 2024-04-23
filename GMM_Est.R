GMM_Est = function(Y,X,Z,k,l,W,r=0){
  # require tidyverse
  #X : instrument,   Z: regressors,   W: weight, r : restriction on coeff(1, same coeff)
  #k: seq of number of instruments, l:seq of number of regressors
  if(r==1){ # same coeff
    nrowY = nrow(Y) #obs
    ncolY = ncol(Y) #dim
    nrowX = nrow(X)
    ncolX = ncol(X) # sum(K)
    nrowZ = nrowX
    ncolZ = ncol(Z) # sum(L)
    
    for(i in 1:length(k)){
      if(i==1){
        block = X[,1:k[i]]
        block01 = Z[,1:l[i]]
        EsigmaXY  = 1/nrowX*t(block)%*%Y[1:nrowX,]
        ESSigmaXZ = 1/nrowX*t(block)%*%block01
      }else{
        block = X[,(sum(k[1:(i-1)])+1):sum(k[1:i])]
        block01 = Z[,(sum(l[1:(i-1)])+1):sum(l[1:i])]
        EsigmaXY=rbind(EsigmaXY,1/nrowX*t(block)%*%Y[(nrowX*(i-1)+1):(nrowX*i),])
        ESSigmaXZ = rbind(ESSigmaXZ,1/nrowX*t(block)%*%block01)
      }
    }
    GMM.coeff = solve(t(ESSigmaXZ)%*%W%*%ESSigmaXZ)%*%t(ESSigmaXZ)%*%W%*%EsigmaXY
    
    for(i in 1:length(k)){ # i번째 식, # k[i]: i번째 식의 차원
      
        if(i==1){
          Error = Y[1:nrowZ,] - Z[,1:l[i]]%*%GMM.coeff  #첫번째 식의 에러
          g = 1/nrowZ*diag(as.numeric(Error))%*%(X[,1:k[i]])
        }else{
          
          ErrorNext = Y[(nrowZ*(i-1)+1):(nrowZ*i),] - 
                      Z[,(sum(l[1:(i-1)])+1):sum(l[1:i])]%*%GMM.coeff
          gNext = 1/nrowZ*diag(as.numeric(ErrorNext))%*%(X[,(sum(k[1:(i-1)])+1):sum(k[1:i])])
          Error = cbind(Error,ErrorNext)
          g = cbind(g,gNext)
        }
      
    }
    
    W01 = solve(1/nrowZ*t(g)%*%g)
    
    return(list(GMM.coeff,W01))
    #return(GMM.coeff)
  
  }
  else{
    nrowY = nrow(Y) #obs
    ncolY = ncol(Y) #dim
    nrowX = nrow(X)
    ncolX = ncol(X) # K
    #nrowZ = nrow(Z) 
    ncolZ = ncol(Z) # L
    
    for(i in 1:length(k)){
      if(i==1){
        block = X[,1:k[i]]
        block01 = Z[,1:l[i]]
        EsigmaXY  = 1/nrowX*t(block)%*%Y[1:nrowX,]
        #ESSigmaXZ = 1/nrowX*t(block)%*%block01
        ESSigmaXZ = diag(ESSigmaXZ,1/nrowX*t(block)%*%block01)
      }else{
        block = X[,(sum(k[1:(i-1)])+1):sum(k[1:i])]
        block01 = Z[,(sum(l[1:(i-1)])+1):sum(l[1:i])]
        EsigmaXY=rbind(EsigmaXY,1/nrowX*t(block)%*%Y[(nrowX*(i-1)+1):(nrowX*i),])
        #ESSigmaXZ = rbind(ESSigmaXY,1/nrowX*t(block)%*%block01)
        ESSigmaXZ = bdiag(ESSigmaXZ,1/nrowX*t(block)%*%block01)
      }
      
    }
    GMM.coeff = as.matrix(solve(t(ESSigmaXZ)%*%W%*%ESSigmaXZ )%*%t(ESSigmaXZ)%*%W %*%EsigmaXY)
    
    for(i in 1:length(k)){ # i번째 식, # k[i]: i번째 식의 차원
      
      if(i==1){
        Error = Y[1:nrowZ,] - Z[,1:l[i]]%*%GMM.coeff[1:l[i],1]  #첫번째 식의 에러
        g = 1/nrowZ*diag(as.numeric(Error))%*%(X[,1:k[i]])
      }else{
        
        ErrorNext = Y[(nrowZ*(i-1)+1):(nrowZ*i),] - 
          Z[,(sum(l[1:(i-1)])+1):sum(l[1:i])]%*%GMM.coeff[(sum(l[1:(i-1)])+1):sum(l[1:i]),1]
        gNext = 1/nrowZ*diag(as.numeric(ErrorNext))%*%(X[,(sum(k[1:(i-1)])+1):sum(k[1:i])])
        Error = cbind(Error,ErrorNext)
        g = cbind(g,gNext)
      }
      
    }
       
    W01 = solve(1/nrowZ*t(g)%*%g)
    
    return(GMM.coeff,W01)
  }
}
    
