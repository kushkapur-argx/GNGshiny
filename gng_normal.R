rm(list=ls())

gng_normal<- function(delta, n_pl, n_trt, sd_pl,sd_trt,   
                      prior_mu, prior_sd, 
                      PCT_up, PCT_low, 
                      LRV, TV, 
                      seed = seed,  iters = iters) {

df <- n_pl + n_trt - 2
pop_var_pool <- ((n_pl - 1) * sd_pl**2 + (n_trt - 1) * sd_trt**2) / df
sd_pool <- sqrt((pop_var_pool / n_pl) + (pop_var_pool / n_trt))
var_post <- 1 / (sd_pool**(-2) + prior_sd**(-2))
sd_post <- sqrt(var_post)

delta_pr <- rnorm(iters, delta, sd_pool)
del_post <- var_post * (prior_mu / prior_sd**2 + delta_pr / sd_pool**2)

del_post<-del_post[order(del_post)]
    
    PCT_lowthresh <- qnorm(PCT_low, del_post, sd_post)
    PCT_upthresh <- qnorm(PCT_up, del_post, sd_post)
    
    govec <- { PCT_lowthresh > LRV } & {PCT_upthresh > TV }
    considervec <- { PCT_lowthresh <= LRV } & { PCT_upthresh > TV}
    stopvec <- { PCT_upthresh <= TV }
    
    
    govalue<- ifelse((!((iters-sum(govec))+1)>length(del_post)),as.numeric(del_post[(iters-sum(govec))+1]),c(NA))
    stopvalue<- ifelse((!sum(stopvec)==0),as.numeric(del_post[sum(stopvec)]),c(NA))
  
    result <- list(df=df, pop_sd = sqrt(pop_var_pool), sd_pool=sd_pool, sd_post=sd_post, 
                   delta=delta, d_pr=mean(delta_pr), d_pos=mean(del_post),
                   stop=mean(stopvec), consider=mean(considervec), go=mean(govec), 
                   stopval=stopvalue, goval=govalue)
  return(result)
   
}


#default values
n_p <- 15
n_t <- 15
p_mu <- 0
p_sd <- 10000
LowRV = 10
TarV = 25
PCT_u = 0.90
PCT_l = 0.20

seed<-12345
iters<-10000

pop_sdvec_pl<-seq(15,25,by=5)
pop_sdvec_trt<-seq(15,25,by=5)
deltaseq<-seq(0,35,by=5)

resultdf <- data.frame(df=as.numeric(),pop_sd=as.numeric(),sd_pool=as.numeric(),sd_post=as.numeric(),
                       delta=as.numeric(),d_pr=as.numeric(),d_pos=as.numeric(),
                       stop=as.numeric(),consider=as.numeric(),go=as.numeric(),
                       stopval=as.numeric(), goval=as.numeric())

for(j in 1:length(pop_sdvec_pl)){
  for(i in 1:length(deltaseq)){
  delta_i<-deltaseq[i]
  popsd_pl<-pop_sdvec_pl[j]
  popsd_trt<-pop_sdvec_trt[j]
  res<-gng_normal(delta = delta_i, n_pl = n_p, n_trt =  n_t, sd_pl = popsd_pl, sd_trt = popsd_trt,
                prior_mu = p_mu, prior_sd = p_sd, 
                PCT_up = PCT_u, PCT_low =  PCT_l,
                LRV = LowRV, TV = TarV,
                seed=seed, iters = iters)
  resultdf<-rbind(resultdf,res)
}
}




