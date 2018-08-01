# parameters and likelihood data
# alphas = optim(c(2, -2, -2, -2, -2), fn, data=lkl.frame, method = "Nelder-Mead", control=list(fnscale=-1))

fn = function(p, data) {
  a0 = p[1]
  a1 = p[2]
  a2 = p[3]
  a3 = p[4]

   #p0 <- 1 - (p1 + p2 + p12 + p1*p2)
  suma = sum(exp(c(a0,a1,a2,a3,a1+a2)))

   lkl.frame.temp <- as.matrix(data)
   lkl.frame.temp[,1] <- lkl.frame.temp[,1] + log(exp(a0)/suma)
   lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(a1)/suma)
   lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(a2)/suma)
   lkl.frame.temp[,4] <- lkl.frame.temp[,4] + log(exp(a1)/suma) + log(exp(a2)/suma)
   lkl.frame.temp[,5] <- lkl.frame.temp[,5] + log(exp(a3)/suma)

   sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
   # penalizing by adding a huge number
   #sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum)) - (10^4)*abs(1-p0-p1-p2-p3-p4)

   return(sumlkl)
}

fn.pw.gwas = function(p, data) {
  a0 = p[1]
  a1 = p[2]
  a2 = p[3]
  a3 = p[4]
  a4 = p[5]
  #print(nrow(data))
  suma = sum(exp(c(a0,a1,a2,a3,a4)))
  #print(paste("Alphas:" , exp(a0)/suma,exp(a1)/suma,exp(a2)/suma,exp(a3)/suma,exp(a4)/suma), sep=" ")
  lkl.frame.temp <- as.matrix(data)
  lkl.frame.temp[,1] <- log(exp(a0)/suma)
  lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(a1)/suma)
  lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(a2)/suma)
  lkl.frame.temp[,4] <- lkl.frame.temp[,4] + log(exp(a3)/suma)
  lkl.frame.temp[,5] <- lkl.frame.temp[,5] + log(exp(a4)/suma)
  #print(lkl.frame.temp[1,])
  #print(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  #print(log(sum(exp(lkl.frame.temp[1,]))))
  sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
  sumlkl = sumlkl + sum(log(c(dnorm(p[1], mean=2, sd=3), (dnorm(p[2:5], mean=-2, sd=3)))))
  #print(sumlkl)
  return(sumlkl)
}

fn.three.hyp = function(p, data) {
  a0 = p[1]
  a1 = p[2]
  a2 = p[3]
  a3 = p[4]

   #p0 <- 1 - (p1 + p2 + p12 + p1*p2)
  suma = sum(exp(c(a0,a1,a2,a3)))

   lkl.frame.temp <- as.matrix(data)
   lkl.frame.temp[,1] <- lkl.frame.temp[,1] + log(exp(a0)/suma)
   lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(a1)/suma)
   lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(a2)/suma)
   lkl.frame.temp[,4] <- logsum(c(lkl.frame.temp[,4], lkl.frame.temp[,5])) + log(exp(a3)/suma)
   lkl.frame.temp[,5] <- 1

   sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum))
   # penalizing by adding a huge number
   #sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum)) - (10^4)*abs(1-p0-p1-p2-p3-p4)

   return(sumlkl)
}



fn.per.SNP = function(p, data) {
  p0 = p[1]
  p1 = p[2]
  p2 = p[3]
  p3 = p[4]

  sump = sum(exp(c(p0,p1,p2,p3)))

  lkl.frame.temp <- as.matrix(data)
  lkl.frame.temp[,1] <- lkl.frame.temp[,1] + log(exp(p0)/sump)
  lkl.frame.temp[,2] <- lkl.frame.temp[,2] + log(exp(p1)/sump)
  lkl.frame.temp[,3] <- lkl.frame.temp[,3] + log(exp(p2)/sump)
  lkl.frame.temp[,5] <- lkl.frame.temp[,5] + log(exp(p3)/sump)

  # penalizing by adding a huge number
  sumlkl = sum(apply(lkl.frame.temp, MAR = 1, FUN = logsum)) 
}


