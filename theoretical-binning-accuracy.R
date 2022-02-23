# R script associated with "A Theoretical Analysis of Taxonomic Binning Accuracy".
# Bianca De Sanctis, Feb 23 2022. bdd28@cam.ac.uk. 

# The main function compute.prob() is defined on line 221. It takes as input N (effective population size), mu (mutation rate for the sequence per generation),  Ttf (true-false species divergence time), Tqt (true-query species divergence time), At, Af and Aq (ages of true, false and query sequences), rT and rF (completeness of true and false reference sequences), prob (the type of assignment probability: correct, incorrect or no), and protocol (either least-mismatch or exact-match), and it outputs the assignment probability.


library(ggplot2)
library(gridExtra)


# 1. Defines functions to compute assignment probabilities.  ################

# Important note: The inside integral is very small, and prone to underflow errors. 
# This means the calculation can be very slow.
# To speed it up, I multiply by a scaling constant of 1e30 on the inside of the integral, and divide it out on the outside.
# This works fine for the wide range of parameters we use in the paper.
# However, with N>30000, the current code will be prone to underflow issues and report inaccurate results, regardless of the scaling constant used.


case1 = function(N,mu,Ttf,
                 Tqt=0,At=0,Af=0,Aq=0,
                 prob="correct",protocol="least-mismatch"){
  ### error catching
  if(! protocol %in% c("least-mismatch","exact-match")){print("error: specify an acceptable protocol")}
  if(! prob %in% c("correct","incorrect","no")){print("error: specify an acceptable prob")}
  
  case1inner = function(t1,t2){
    if(protocol == "least-mismatch"){
      # Theoretically, we want to sum this on l=k to inf, k=0 to infinity. 
      # But very quickly as we increase k, the values will be come negligible. So just go up to some number lmax which is sufficiently large.
      # l is for false species, k is for true
      lmax = max(3,round(Ttf*mu*2*10))
      lks = combn(0:lmax,2)
      if( prob == "correct" & protocol == "least-mismatch"){ # l > k
        ls = lks[2,]; ks = lks[1,]
      }
      if( prob == "no" & protocol == "least-mismatch"){ # l = k
        ls = ks = 0:lmax
      }
      if( prob == "incorrect" & protocol == "least-mismatch" ){ # l < k
        ls = lks[1,]; ks = lks[2,]
      }
      innersumfunc = function(l,k){
        mu^(l+k) * (2*t2 - t1-Af)^l * (t1-At)^k * exp(-mu*(2*t2-Af-At)) / (factorial(l)*factorial(k))
      }
      innersumfunc = Vectorize(innersumfunc)
      inner = rowSums(innersumfunc(ls,ks))
    }
    if(protocol == "exact-match"){
      if (prob=="correct"){
        inner = exp(-mu*(2*t1 - At - Aq)) *(1- exp(-mu*(2*t2-t1 - Af)))
      }
      if(prob =="incorrect"){
        inner = (1-exp(-mu*(t1 - At )))*exp(-mu*(2*t2 - Af - Aq))
      }
      if(prob=="no"){
        inner = 1 - exp(-mu*(2*t1 - At - Aq)) *(1- exp(-mu*(2*t2-t1 - Af))) - (1-exp(-mu*(t1 - At )))*exp(-mu*(2*t2 - Af - Aq))
      }
    }
  front = 1e30* 1/(2*N)^2 * exp(-t1/(2*N)) * exp (-t2/(2*N)) # 1e30 avoids underflow errors
  return(inner*front)
  }
  
  
  F1 = function(t2) {
    fun <- function(t1) case1inner(t1,t2) 
    integrate(fun,Tqt,Ttf)$value # inner integral limits
  }
  F1 = Vectorize(F1)
  answer = integrate(F1,Ttf,Inf)$value/1e30 # outer integral limits
  denominator = exp(-Ttf/(2*N)) * ( exp(-Tqt/(2*N)) - exp(-Ttf/(2*N))   )
  return(answer/denominator)
}

case2.1 = function(N,mu,Ttf,
                   Tqt=0,At=0,Af=0,Aq=0,
                   prob="correct",protocol="least-mismatch"){
  case2.1inner = function(t1,t2){
    if(protocol == "least-mismatch"){
      # Theoretically, we want to sum this on l=k to inf, k=0 to infinity. 
      # But very quickly as we increase k, the values will be come negligible. So just go up to some number lmax which is sufficiently large.
      # l is for false species, k is for true
      lmax = max(3,round(Ttf*mu*2*10))
      lks = combn(0:lmax,2)
      if( prob == "correct" & protocol == "least-mismatch"){ # l > k
        ls = lks[2,]; ks = lks[1,]
      }
      if( prob == "no" & protocol == "least-mismatch"){ # l = k
        ls = ks = 0:lmax
      }
      if( prob == "incorrect" & protocol == "least-mismatch" ){ # l < k
        ls = lks[1,]; ks = lks[2,]
      }
      innersumfunc = function(l,k){
        mu^(l+k) * (2*t2 - t1-Af)^l * (t1-At)^k * exp(-mu*(2*t2-Af-At)) / (factorial(l)*factorial(k))
      }
      innersumfunc = Vectorize(innersumfunc)
      inner = rowSums(innersumfunc(ls,ks))
    }
    if(protocol == "exact-match"){
      if (prob=="correct"){
        inner = exp(-mu*2*t1) *(1- exp(-mu*(2*t2-t1)))
      }
      if(prob =="incorrect"){
        inner = (1-exp(-mu*t1))*exp(-mu*(2*t2))
      }
      if(prob=="no"){
        inner = 1 - exp(-mu*2*t1) *(1- exp(-mu*(2*t2-t1))) - (1-exp(-mu*t1))*exp(-mu*(2*t2))
      }
    }
    front = 1e30* 1/(2*N)^2 * exp(-t1/(2*N)) * exp (-t2/(2*N)) # 1e30 avoids underflow errors
    return(inner*front)
  }
  F1 = function(t1) { # <- outer integral dt
    fun <- function(t2) case2.1inner(t1,t2)  # <- inner integral dt
    integrate(fun,t1,Inf)$value # inner integral limits
  }
  F1 = Vectorize(F1)
  answer = integrate(F1,Ttf,Inf)$value/1e30 # outer integral limits
  denominator = exp(-Ttf/(2*N)) 
  return(answer/denominator)
}

case2.2 = function(N,mu,Ttf,
                   Tqt=0,At=0,Af=0,Aq=0,
                   prob="correct",protocol="least-mismatch"){
  case2.2inner = function(t1,t2){
    if(protocol == "least-mismatch"){
      lmax = max(3,round(Ttf*mu*2*3))
      lks = combn(0:lmax,2)
      if( prob == "correct" & protocol == "least-mismatch"){ # l > k
        ls = lks[2,]; ks = lks[1,]
      }
      if( prob == "no" & protocol == "least-mismatch"){ # l = k
        ls = ks = 0:lmax
      }
      if( prob == "incorrect" & protocol == "least-mismatch" ){ # l < k
        ls = lks[1,]; ks = lks[2,]
      }
      innersumfunc = function(l,k){
        # t1 is vectorized, so doing t1-t1+1 just gets it into the right dimensions
        (t1-t1+1) * mu^(l+k) * (t2-Af)^l * (t2-At)^k * exp(-mu*(2*t2-Af-At)) / (factorial(l)*factorial(k))
      }
      innersumfunc = Vectorize(innersumfunc)
      inner = rowSums(innersumfunc(ls,ks))
    }
    if(protocol == "exact-match"){
      if (prob=="correct"){
        inner = exp(-mu*(2*t1 - At - Aq)) *(1- exp(-mu*(t2 - Af)))
      }
      if(prob =="incorrect"){
        inner = exp(-mu*(2*t1 - Af - Aq)) * (1-exp(-mu*( t2 - At )))
      }
      if(prob=="no"){
        inner = 1 - exp(-mu*(2*t1 - At - Aq)) *(1- exp(-mu*(t2 - Af))) - exp(-mu*(2*t1 - Af - Aq)) * (1-exp(-mu*( t2 - At )))
      }
    }
    
    front = 1e30* 1/(2*N)^2 * exp(-t1/(2*N)) * exp (-t2/(2*N)) # 1e30 avoids underflow errors
    return(inner*front)
  }
  F1 = function(t2) { # <- outer integral dt
    fun <- function(t1) case2.2inner(t1,t2)  # <- inner integral dt
    integrate(fun,t2,Inf)$value # inner integral limits
  }
  F1 = Vectorize(F1)
  answer = integrate(F1,Ttf,Inf)$value/1e30 # outer integral limits
  denominator = exp(-Ttf/(2*N)) 
  return(answer/denominator)
}

case2.3 = function(N,mu,Ttf,
                   Tqt=0,At=0,Af=0,Aq=0,
                   prob="correct",protocol="least-mismatch"){
  case2.3inner = function(t1,t3){
    if(protocol == "least-mismatch"){
      lmax = max(3,round(Ttf*mu*2*3))
      lks = combn(0:lmax,2)
      if( prob == "correct" & protocol == "least-mismatch"){ # l > k
        ls = lks[2,]; ks = lks[1,]
      }
      if( prob == "no" & protocol == "least-mismatch"){ # l = k
        ls = ks = 0:lmax
      }
      if( prob == "incorrect" & protocol == "least-mismatch" ){ # l < k
        ls = lks[1,]; ks = lks[2,]
      }
      innersumfunc = function(l,k){
        # the t1-t1+1 at the beginning is because t1 is vectorized, and this just gets it into the right dimensions
        mu^(l+k) * (t3-Af)^l * (2*t1-t3-At)^k * exp(-mu*(2*t1-Af-At)) / (factorial(l)*factorial(k))
      }
      innersumfunc = Vectorize(innersumfunc)
      inner = rowSums(innersumfunc(ls,ks))
    }
    
    if(protocol == "exact-match"){
      if (prob=="correct"){
        inner = exp(-mu*(2*t3 - At - Aq)) *(1- exp(-mu*(t2 - Af)))
      }
      if(prob =="incorrect"){
        inner = exp(-mu*(2*t2 - Af - Aq)) * (1-exp(-mu*(2*t3-t2 - At)))
      }
      if(prob=="no"){
        inner = 1 - exp(-mu*(2*t3 - At - Aq)) *(1- exp(-mu*(t2 - Af))) - exp(-mu*(2*t2 - Af - Aq)) * (1-exp(-mu*(2*t3-t2 - At)))
      }
    }
      
    front = 1e30* 1/(2*N)^2 * exp(-t1/(2*N)) * exp (-t3/(2*N)) # 1e30 avoids underflow errors
    return(inner*front)
  }
  F1 = function(t3) { # <- outer integral dt
    fun <- function(t1) case2.3inner(t1,t3)  # <- inner integral dt
    integrate(fun,t3,Inf)$value # inner integral limits
  }
  F1 = Vectorize(F1)
  answer = integrate(F1,Ttf,Inf)$value/1e30 # outer integral limits
  denominator = exp(-Ttf/(2*N)) 
  return(answer/denominator)
}

compute.prob <- function(N,mu,Ttf,
                         Tqt=0,At=0,Af=0,Aq=0,rT=1,rF=1,
                         prob="correct",protocol="least-mismatch"){
  prob.case1 = (1-exp(-min(Ttf - Aq , Ttf - At)/(2*N))) 
  if(prob.case1 > 0.9999){
    # no need to compute case 2
    out = case1(N,mu,Ttf,Tqt,At,Af,Aq,prob,protocol)
  }
  else{
    # need to compute case 2
    print("running case 2") 
    prob.case2 = 1-prob.case1
    out = prob.case1 * case1(N,mu,Ttf,Tqt,At,Af,Aq,prob,protocol) + 
      (1/3) * prob.case2 * ( case2.1(N,mu,Ttf,Tqt,At,Af,Aq,prob,protocol) + 
                               case2.2(N,mu,Ttf,Tqt,At,Af,Aq,prob,protocol) +
                               case2.3(N,mu,Ttf,Tqt,At,Af,Aq,prob,protocol) )
  }
  # account for representation here
  if(prob=="correct"){final = (rT * rF * out + rT*(1-rF)) / (1 - (1-rT)*(1-rF))}
  if(prob=="incorrect"){final = (rT * rF * out + rF*(1-rT)) / (1 - (1-rT)*(1-rF))}
  if(prob=="no"){final = (rT * rF * out) / (1 - (1-rT)*(1-rF))}
  return(final)
}



# minimum example: 
# compute.prob(N = 10000,mu = 1e-8,Ttf = 4e5)



# 2. Creates the main paper figure. #####

N.list = seq(from=4000,to=20000,by=2000)
Ttf.list = seq(from=2e5,to=10e5,by = 2e5)
k.list = c(32,64,96,128,160)
rT.list=c(0.2,0.4,0.6,0.8,1)
Tqt.list = seq(from = 0, to = 3e5, by=5e4) 

# base parameter choices
mu.base = 1e-8 
N.base = N.list[4] 
Ttf.base = Ttf.list[2]
k.base = k.list[3]
rT.base = rT.list[5]
Tqt.base = Tqt.list[1]

protocol.list = c("least-mismatch","exact-match")
assignment.list = c("correct","incorrect","no")
tempN = cbind(expand.grid(assignment.list,protocol.list,N.list,stringsAsFactors=FALSE),
              Ttf=Ttf.base,k=k.base,rT=rT.base,Tqt=Tqt.base,
              stringsAsFactors=FALSE)
tempTtf = cbind(expand.grid(assignment.list,protocol.list,Ttf.list,stringsAsFactors=FALSE),
                N=N.base,k=k.base,rT=rT.base,Tqt=Tqt.base,stringsAsFactors=FALSE)[,c(1,2,4,3,5,6,7)]
tempk = cbind(expand.grid(assignment.list,protocol.list,k.list,stringsAsFactors=FALSE),
              N=N.base,Ttf=Ttf.base,rT=rT.base,Tqt=Tqt.base,stringsAsFactors=FALSE)[,c(1,2,4,5,3,6,7)]
temprT = cbind(expand.grid(assignment.list,protocol.list,rT=rT.list,stringsAsFactors=FALSE),
               N=N.base,Ttf=Ttf.base,k=k.base,Tqt=Tqt.base,stringsAsFactors=FALSE)[,c(1,2,4,5,6,3,7)]
tempTqt = cbind(expand.grid(assignment.list,protocol.list,Tqt=Tqt.list,stringsAsFactors=FALSE),
                N=N.base,Ttf=Ttf.base,k=k.base,rT=rT.base,stringsAsFactors=FALSE)[,c(1,2,4,5,6,7,3)]

colnames(tempN) = colnames(tempTtf) = colnames(tempk) = 
  colnames(temprT) = colnames(tempTqt) =  c("assignment","protocol","N","Ttf","k","rT","Tqt")
temp = rbind(tempN,tempTtf,tempk,temprT,tempTqt,stringsAsFactors=FALSE)
temp = cbind(temp,rep(-1,nrow(temp)),stringsAsFactors=FALSE)
out.table = data.frame(temp,stringsAsFactors=FALSE)
colnames(out.table) = c("assignment","protocol","N","Ttf","k","rT","Tqt","probability")
for(row in 1:nrow(out.table)){
  tryCatch({out.table[row,8] = compute.prob(N = out.table$N[row],
                                            mu = mu.base * out.table$k[row],
                                            Ttf = out.table$Ttf[row] ,
                                            prob = out.table$assignment[row],
                                            rT = out.table$rT[row],
                                            Tqt = out.table$Tqt[row],
                                            protocol = out.table$protocol[row])})
} 

col1="#009E73"   # green
col2="#D55E00"  # red
col3="#F0E442"  # yellow
col4="dodgerblue1"

### fig 1a. varying N. unconditional.
out.table$N = out.table$N / 1000 
out.table$Ttf = out.table$Ttf / 1e5
out.table$Tqt = out.table$Tqt / 1e5

relN = out.table[1:27,] 
relN$N = factor(relN$N,levels=unique(relN$N))
relN$protocol = ifelse(relN$protocol=="exact-match","EM","LM")
G1a = ggplot(relN, aes(x = N, y = probability, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') +    
  scale_fill_manual(values = c(col1, col2, 
                               col3)) +  
  labs(x="Effective population size (thousands)",y="")+
  theme(legend.position="none")


### fig 1b. varying N. conditional.
cor = relN[relN$assignment == "correct",]
incor = relN[relN$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.relN = relN
cond.relN = cond.relN[cond.relN$assignment %in% c("correct","incorrect"),]
cond.relN$probability = c(rbind(cond.cor,cond.incor))
G1b_old = ggplot(cond.relN, aes(x = protocol, y = probability, fill = assignment)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ N) +   ggtitle(paste0(
    "Conditional probabilities, varying N, 
    kmer length = ",k.base,", Ttf = ",Ttf.base,", mu = ",mu.base)) 
ratio.cor = cor$probability / incor$probability
ratio.relN = relN
ratio.relN = ratio.relN[ratio.relN$assignment == "correct",]
ratio.relN$probability = ratio.cor
colnames(ratio.relN)[8] = "count"
G1b = ggplot(ratio.relN, aes(x = N, y = count, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = c(col4, col2, 
                               col3)) + 
  labs(x="Effective population size (thousands)",y="")+
  theme(legend.position="none")


### fig 2a. varying Ttf. unconditional.
relTtf = out.table[28:42,]# out.table[(out.table$k == 32) & (out.table$N == 10000),] 
# relTtf[,4] = sapply(relTtf[,4], function(x) paste0("Ttf=",x) )
relTtf$Ttf = factor(round(relTtf$Ttf),levels=unique(round(relTtf$Ttf)))
relTtf$protocol = ifelse(relTtf$protocol=="exact-match","U","LM")
G2a = ggplot(relTtf, aes(x = Ttf, y = probability, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') +    
  scale_fill_manual(values = c(col1, col2, 
                               col3)) + 
  labs(x="True-false species divergence (1e5 generations)",y="")+
  theme(legend.position="none")

### fig 2b. varying Ttf. conditional.
cor = relTtf[relTtf$assignment == "correct",]
incor = relTtf[relTtf$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.relTtf = relTtf
cond.relTtf = cond.relTtf[cond.relTtf$assignment %in% c("correct","incorrect"),]
cond.relTtf$probability = c(rbind(cond.cor,cond.incor))
G2b_old = ggplot(cond.relTtf, aes(x = protocol, y = probability, fill = assignment)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ Ttf) +   ggtitle(paste0(
    "Conditional probabilities, varying divergence times, 
    read length = ",k.base,", N = ",N.base,", mu = ",mu.base)) 
ratio.cor = cor$probability / incor$probability
ratio.relTtf= relTtf
ratio.relTtf = ratio.relTtf[ratio.relTtf$assignment == "correct",]
ratio.relTtf$probability = ratio.cor
colnames(ratio.relTtf)[8] = "count"
G2b = ggplot(ratio.relTtf, aes(x = Ttf, y = count, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = c(col4, col2, 
                               col3)) +  
  labs(x="True-false species divergence (1e5 generations)",y="")+
  theme(legend.position="none")


### fig 3a. varying k. unconditional.
relk = out.table[43:57,] 
relk$k = factor(round(relk$k),levels=unique(round(relk$k)))
relk$protocol = ifelse(relk$protocol=="exact-match","U","LM")
G3a = ggplot(relk, aes(x = k, y = probability, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') +    
  scale_fill_manual(values = c(col1, col2, 
                               col3)) +  
  labs(x="Query sequence length",y="")+
  theme(legend.position="none")

### fig 3b. varying k. conditional.
cor = relk[relk$assignment == "correct",]
incor = relk[relk$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.relk = relk
cond.relk = cond.relk[cond.relk$assignment %in% c("correct","incorrect"),]
cond.relk$probability = c(rbind(cond.cor,cond.incor))
G3b_old = ggplot(cond.relk, aes(x = protocol, y = probability, fill = assignment)) + 
  geom_bar(stat = 'identity', position = 'stack') + facet_grid(~ k) +   ggtitle(paste0(
    "Conditional probabilities, varying k, 
    Ttf = ",Ttf.base,"N =" ,N.base,", mu = ",mu.base))  + theme(text = element_text(size=10))
ratio.cor = cor$probability / incor$probability
ratio.relk= relk
ratio.relk = ratio.relk[ratio.relk$assignment == "correct",]
ratio.relk$probability = ratio.cor
colnames(ratio.relk)[8] = "count"
G3b = ggplot(ratio.relk, aes(x = k, y = count, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = c(col4, col2, 
                               col3))+ 
  labs(x="Query sequence length",y="")+
  theme(legend.position="none")





### fig 4a. varying completeness. unconditional.
relrT = out.table[58:72,] 
relrT$rT = factor(relrT$rT,levels=unique(relrT$rT))
relrT$protocol = ifelse(relrT$protocol=="exact-match","U","LM")
G4a = ggplot(relrT, aes(x = rT, y = probability, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') +    
  scale_fill_manual(values = c(col1, col2, 
                               col3))+
  labs(x="True sequence completeness",y="")+
  theme(legend.position="none")

### fig 3b. varying completeness. conditional.
cor = relrT[relrT$assignment == "correct",]
incor = relrT[relrT$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.relrT = relrT
cond.relrT = cond.relrT[cond.relrT$assignment %in% c("correct","incorrect"),]
cond.relrT$probability = c(rbind(cond.cor,cond.incor))
ratio.cor = cor$probability / incor$probability
ratio.relrT= relrT
ratio.relrT = ratio.relrT[ratio.relrT$assignment == "correct",]
ratio.relrT$probability = ratio.cor
colnames(ratio.relrT)[8] = "count"
G4b = ggplot(ratio.relrT, aes(x = rT, y = count, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = c(col4, col2, 
                               col3))+ 
  labs(x="True sequence completeness",y="")+
  theme(legend.position="none")



### fig 5a. varying Tqt. unconditional.

relTqt = out.table[73:93,] 
relTqt$Tqt = factor(relTqt$Tqt,levels=unique(relTqt$Tqt))
relTqt$protocol = ifelse(relTqt$protocol=="exact-match","U","LM")
G5a = ggplot(relTqt, aes(x = Tqt, y = probability, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') +    
  scale_fill_manual(values = c(col1, col2, 
                               col3))+ 
  labs(x="True-query species divergence (1e5 generations)",y="")+
  theme(legend.position="none")

### fig 5b. varying Tqt. conditional.
cor = relTqt[relTqt$assignment == "correct",]
incor = relTqt[relTqt$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.relTqt = relTqt
cond.relTqt = cond.relTqt[cond.relTqt$assignment %in% c("correct","incorrect"),]
cond.relTqt$probability = c(rbind(cond.cor,cond.incor))
ratio.cor = cor$probability / incor$probability
ratio.relTqt= relTqt
ratio.relTqt = ratio.relTqt[ratio.relTqt$assignment == "correct",]
ratio.relTqt$probability = ratio.cor
colnames(ratio.relTqt)[8] = "count"
G5b = ggplot(ratio.relTqt, aes(x = Tqt, y = count, fill = assignment)) + 
  theme_set(theme_bw() + theme(legend.position="none")) + 
  geom_bar(stat = 'identity', position = 'stack') + 
  scale_fill_manual(values = c(col4, col2, 
                               col3))+
  labs(x="True-query species divergence (1e5 generations)",y="") 


## Plot

# Plot to device
grid.arrange(G1a,G1b,G2a,G2b,G5a,G5b,G3a,G3b,G4a,G4b,
             nrow= 5,ncol=2)


# 3. Creates a figure for the example when the two protocols differ. ####

mu.base = 1e-8
N.base = 10000
Ttf.base =  1e6
k.base = 160
rT.base = 1
Tqt.base =  3e5

protocol.list = c("least-mismatch","exact-match")
assignment.list = c("correct","incorrect","no")
temp = expand.grid(assignment.list,protocol.list,N.base,Ttf=Ttf.base,
                   k=k.base,rT=rT.base,Tqt=Tqt.base,stringsAsFactors=FALSE)
colnames(temp) =  c("assignment","protocol","N","Ttf","k","rT","Tqt")
temp = cbind(temp,rep(-1,nrow(temp)))
out.table = data.frame(temp,stringsAsFactors=FALSE)
colnames(out.table) = c("assignment","protocol","N","Ttf","k","rT","Tqt","probability")
for(row in 1:nrow(out.table)){
  tryCatch({out.table[row,8] = compute.prob(N = out.table$N[row],
                                            mu = mu.base * out.table$k[row],
                                            Ttf = out.table$Ttf[row] ,
                                            prob = out.table$assignment[row],
                                            rT = out.table$rT[row],
                                            Tqt = out.table$Tqt[row],
                                            protocol = out.table$protocol[row])})
}

reld = out.table
reld$protocol = ifelse(reld$protocol=="exact-match","Exact match","Least mismatch")
dplot = ggplot(reld, aes(x = protocol, y = probability, fill = assignment)) + 
  theme_set(theme_bw() +theme(legend.position="none") ) +
  geom_bar(stat = 'identity', position = 'stack')  +  
  ggtitle(paste0("Probability of assignment")) +
  scale_fill_manual(values = c(col1, col2, col3)) + labs(x="",y="")

cor = reld[reld$assignment == "correct",]
incor = reld[reld$assignment == "incorrect",]
cond.cor = cor$probability / (cor$probability + incor$probability)
cond.incor = 1 - cor$probability / (cor$probability + incor$probability)
cond.reld = reld
cond.reld = cond.reld[cond.reld$assignment %in% c("correct","incorrect"),]
cond.reld$probability = c(rbind(cond.cor,cond.incor))
ratio.cor = cor$probability / incor$probability
ratio.reld = reld
ratio.reld = ratio.reld[ratio.reld$assignment == "correct",]
ratio.reld$probability = ratio.cor
colnames(ratio.reld)[8] = "count"
dplot2 = ggplot(ratio.reld, aes(x = protocol, y = count, fill = assignment)) + 
  theme_set( theme_bw() +theme(legend.position="none") ) +
  geom_bar(stat = 'identity', position = 'stack')  +   ggtitle(paste0(
    "Correct assignments per incorrect assignment")) +
  scale_fill_manual(values = c(col4)) + labs(x="",y="")

grid.arrange(dplot,dplot2,nrow= 1) 







