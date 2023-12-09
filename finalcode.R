rm(list=ls(all=T))

# Backbone code. This model involves the removal of fish from a pond.

#max1 fun, makes it so prob is not greater than 1
max1<-function(a,b){
  c=sum(a,b)
  d=min(c,1)
  return(d)
}

#function to get beta parameters 
BetaParams <- function(mu, var) {
  alpha <- mu*((mu*(1-mu)/var-1))
  beta <- (1-mu)*((mu*(1-mu)/var-1))
  return(c(alpha = alpha, beta = beta))
}

Replicates = function(sd_beta_ind_v, alpha_v,vr_V) {  
  
  replicates = expand.grid(sd_beta_ind_v, alpha_v,vr_V)
  colnames(replicates) = c('sd_beta_ind_v', 'alpha_v','vr_V')
  
  return(replicates)
} 

#growth parameters 
k=0.251 #growth parameter
sd_k=0.002 #var in K
sd_k_ind=0.0005
t0=-0.207 #time at age 0
linf=650 #asymptotic length

#RPS-recruits per spawner
lambda=1
RPS=rpois(1,lambda) 
periods=50 #number of weeks prior to spawning

sd_beta_indiv_v=c(0,0.005,0.01)
vr_v=c(0.5,0.75,1)
alpha_v=c(0,0.1,0.2)
reps=Replicates(sd_beta_indiv_v,alpha_v,vr_v)

#more containers
nrep=10
final_tot_beta_y=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_n=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_mean_length=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_mean_age=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_caught=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_har=array(NA,dim=c(nrow(reps),nrep,periods))
final_tot_trophy=array(NA,dim=c(nrow(reps),nrep,periods))

tot_beta_y=matrix(NA,nrep,periods)
tot_n=matrix(NA,nrep,periods)
tot_mean_length=matrix(NA,nrep,periods)
tot_mean_age=matrix(NA,nrep,periods)
tot_caught=matrix(NA,nrep,periods)
tot_har=matrix(NA,nrep,periods)
tot_trophy=matrix(NA,nrep,periods)


for(q in 1:nrow(reps)){
start=Sys.time()
for(j in 1:nrep){
a_NS=0.1 #near shore capture prob intercept yearly
a_OS=0.01 #offshore capture prob intercept yearly
a_M=0.3 #natural M yearly
M=a_M^1/periods
vr=reps$vr_V[q] #voluntary release
alpha=reps$alpha_v[q] #post release M
PR=alpha+M #post release M + natural M
nstates=6
a=c(a_NS,a_OS)^1/periods #just allows me. to index alpha, puts on smaller timestep
HP=1 #harvest rate under MLM
MLM=300 #harvest every fish under MLM
#states
# 1.  Not Caught alive
# 2. C+R survivor
# 3. Natural M
# 4. C+R Death
# 5. Harvested
# 6. term death

#initialize
n=100
cc=250 #carrying capacity
age=rep(2,n)# stock all fish at 2 years old, say 100 survive rpois(n,4) #initial age
age[age<2]<-2
age[age>12]<-12
sex=rbinom(n,1,0.5) #0=M 1=F
k_beta=BetaParams(mu=k,var=sd_k)
k_ind=rbeta(n, k_beta[1],k_beta[2]) #var in growth
lgth=linf*(1-exp(-k_ind*(age-t0))) #length from von. berf
mean_beta=0.3
sd_beta=0.05
sd_beta_indiv=reps$sd_beta_ind_v[q]
temp_par=BetaParams(mean_beta,sd_beta)
beta_yearly=rbeta(n,temp_par[1],temp_par[2]) #yearly indiv prob of capture
beta=beta_yearly^1/periods
hab=rbinom(n,1,0.5) #habitat 0=NS 1=OS
dat=data.frame(id=1:n,age=age,lgth=lgth,sex=sex+1,k_ind=k_ind,beta=beta,beta_yearly=beta_yearly,hab=hab+1)
dat


spawners=list(NULL,NULL) #lag spawners by 2 years
pop_size=NULL
mean_beta_y=NULL
mean_age=NULL
mean_length=NULL
fish_caught=NULL
fish_har=NULL
trophy_caught=NULL




for(rep in 3:52){ #50 years, starting at. 3 works out a bug. need to think about that
  year=matrix(NA,nrow(dat),periods)
  vr2=year #voluntary release, we want to harvest every fish <300 mm
  lgth_year=year
  year[,1]=dat$age
  lgth_year[,1]=dat$lgth
  for(i in 1:nrow(dat)){
    for(t in 2:periods){
      year[,t]=year[,1]+t/periods
      lgth_year[i,t]=linf*(1-exp(-dat$k_ind[i]*(year[i,t]-t0)))
    }
  }
  
  vr2[lgth_year<MLM]<-0
  vr2[lgth_year>=MLM]=vr
BPM=array(NA,dim=c(nrow(dat),periods,nstates,nstates)) #transition prob
for(i in 1:nrow(dat)){
  for(t in 1:periods){ #don't need t now, give. me. the option to add time varying. params later
    BPM[i,t,,]<-matrix(c(
      (1-M)*(1-max1(a[dat$hab[i]],dat$beta[i])),max1(a[dat$hab[i]],dat$beta[i])*vr2[i,t]*(1-PR),
      M*(1-max1(a[dat$hab[i]],dat$beta[i])),max1(a[dat$hab[i]],dat$beta[i])*vr2[i,t]*PR,max1(a[dat$hab[i]],dat$beta[i])*(1-vr2[i,t]),0,
      (1-M)*(1-max1(a[dat$hab[i]],dat$beta[i])),max1(a[dat$hab[i]],dat$beta[i])*vr2[i,t]*(1-PR),
      M*(1-max1(a[dat$hab[i]],dat$beta[i])),max1(a[dat$hab[i]],dat$beta[i])*vr2[i,t]*PR,max1(a[dat$hab[i]],dat$beta[i])*(1-vr2[i,t]),0,
      0,0,0,0,0,1,
      0,0,0,0,0,1,
      0,0,0,0,0,1,
      0,0,0,0,0,1
    ),nrow=nstates,byrow=T)
  }
}
z=matrix(NA,nrow(dat),periods)
z[,1]=1
# propagate process via pre defined transitions
for(i in 1:nrow(dat)){
  for(t in 2:periods){
    departure.state=z[i,t-1]
    arrival.state=which(rmultinom(1,1,BPM[i,t-1,departure.state,])==1)
    z[i,t]<-arrival.state
  }
}

fish_caught=c(fish_caught,length(z[z==2|z==4|z==5]))
fish_har=c(fish_har,length(z[z==5]))
trophy_caught=c(trophy_caught,length(z[(z==2|z==4|z==5)&(lgth_year>=550)]))

survivors=dat[z[,periods]<3,]
survivorsf=survivors[survivors$sex==1,]
survivorsm=survivors[survivors$sex==2,]

#have survivors reproduce
pairs=min(nrow(survivorsf),nrow(survivorsm))
mates_f=sample(survivorsf$id,pairs)
mates_m=sample(survivorsm$id,pairs,replace = T) #males can mate with more than 1 female


tdat=dat[z[,periods]<3,] #remove dead fish
tdat$id=1:nrow(tdat)
tdat$age=tdat$age+1
tdat<-tdat[!tdat$age>12,]
tdat$lgth=linf*(1-exp(-tdat$k_ind*(tdat$age-t0)))




#containers
recruits=data.frame(id=0,age=0,lgth=0,sex=0,k_ind=0,beta=0,hab=0)
recruits
id=NULL
age=NULL
k_ind=NULL
lgth=NULL
sex=NULL
hab=NULL
beta=NULL
beta_yearly=NULL


for(p in 1:pairs){
  Rec=rpois(1,lambda)
  for(r in 1:Rec){
    age=c(age,2)
    temp=BetaParams(mean(dat$k_ind[dat$id==mates_f[p]],dat$k_ind[dat$id==mates_m[p]]),sd_k_ind)
    k_rec=rbeta(1,temp[1],temp[2])
    k_rec=max(0.1,k_rec) #constrain so k doesn't get too extreme
    k_rec=min(0.4,k_rec) #constrain so k doesn't get too extreme
    k_ind=c(k_ind,k_rec)
    lgth=c(lgth,linf*(1-exp(-k*(2-t0))))
    sex=c(sex,rbinom(1,1,0.5)+1)
    hab=c(hab,rbinom(1,1,0.5)+1)
    mean_p_beta_yearly=mean(dat$beta_yearly[dat$id==mates_f[p]],dat$beta_yearly[dat$id==mates_m[p]]) #yearly beta average
    beta_yearly=c(beta_yearly,max(min(rnorm(1,mean_p_beta_yearly,sd_beta_indiv),1),0.01))
    beta=beta_yearly^1/periods
    recruits=data.frame(id=(1:length(age)),age=age,lgth=lgth,sex=sex,k_ind=k_ind,beta=beta,beta_yearly=beta_yearly,hab=hab)
  }
}
spawners=c(spawners,list(recruits))
#r=c(r,list(recruits))
#recruits$id=recruits$id+nrow(tdat)
#caught up on way to lag recruits, need to think that through....
maxrecruits=cc-nrow(tdat) #maximum number of recruits
thisyearspawners=as.data.frame(spawners[[rep-2]]) # I think this will work
need=min(max(nrow(thisyearspawners),0),maxrecruits)
addspawners=thisyearspawners[sample(1:need,need),]
dat=rbind(tdat,addspawners)
dat<-as.data.frame(dat)
dat$id=1:length(dat$id)
#need to save catch,  mean vuln etc. but right now process is working,  sorta.  I think there is a bug somewhere

#summary
pop_size=c(pop_size,nrow(dat))
mean_beta_y=c(mean_beta_y,mean(dat$beta_yearly))
mean_length=c(mean_length,mean(dat$lgth))
mean_age=c(mean_age,mean(dat$age))
}

tot_n[j,]=pop_size
tot_beta_y[j,]=mean_beta_y
tot_mean_length[j,]=mean_length
tot_mean_age[j,]=mean_age
tot_caught[j,]=fish_caught
tot_har[j,]=fish_har
tot_trophy[j,]=trophy_caught

}
  final_tot_n[q,,]=tot_n
  final_tot_beta_y[q,,]=tot_beta_y
  final_tot_mean_length[q,,]=tot_mean_length
  final_tot_mean_age[q,,]=tot_mean_age
  final_tot_caught[q,,]=tot_caught
  final_tot_har[q,,]=tot_har
  final_tot_trophy[q,,]=tot_trophy
  
  print(c(q,'rep out of',nrow(reps)))
  end=Sys.time()
  time.elapsed=end-start
  print(c(as.numeric(time.elapsed)*(nrow(reps)-q),'min remaining'))
}

plot(tot_n[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_n[i,],type='l',col=i)
}

plot(tot_beta_y[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_beta_y[i,],type='l',col=i)
}

plot(tot_mean_length[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_mean_length[i,],type='l',col=i)
}

plot(tot_mean_age[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_mean_age[i,],type='l',col=i)
}

plot(tot_caught[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_caught[i,],type='l',col=i)
}

plot(tot_har[1,],type='l')
for(i in 2:nrep){
  lines(1:periods,tot_har[i,],type='l',col=i)
}

quartz()
par(mfrow=c(3,3))

for(r in 19:27){
plot(final_tot_trophy[r,1,],type='l',ylim=c(0,40),ylab = "trophy fish",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
axis(1,lwd=2)
axis(2,lwd=2)  
box(lwd=2)
for(i in 2:nrep){
  lines(1:periods,final_tot_trophy[r,i,],type='l',col=i)
  axis(1,lwd=2)
axis(2,lwd=2)  
box(lwd=2)
}
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_beta_y[r,1,],type='l',ylim=c(0,0.5),ylab = "vuln",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_beta_y[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_har[r,1,],type='l',ylim=c(0,20),ylab = "harvest",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_har[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_caught[r,1,],type='l',ylim=c(0,130),ylab = "fish caught",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_caught[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_mean_length[r,1,],type='l',ylim=c(325,475),ylab = "length",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_mean_length[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3),no.readonly = F)
for(r in 1:9){
  plot(colMeans(final_tot_mean_age[r,,]),type='l',ylim=c(3,5),ylab = "",xlab="",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r],"vr=",reps$vr_V[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_mean_age[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}
mtext("age",line=0,side=3,outer=T)

quartz()
par(mfrow=c(3,3))
for(r in 10:18){
  plot(final_tot_mean_age[r,1,],type='l',ylim=c(3,5),ylab = "age",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r],"vr=",reps$vr_V[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_mean_age[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_mean_age[r,1,],type='l',ylim=c(3,5),ylab = "age",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r],"vr=",reps$vr_V[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_mean_age[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

quartz()
par(mfrow=c(3,3))
for(r in 19:27){
  plot(final_tot_n[r,1,],type='l',ylim=c(0,250),ylab = "number of fish",xlab="year",main=paste("sd_b=",reps$sd_beta_ind_v[r],"a=",reps$alpha_v[r]),cex.lab=1.5,cex.main=1.5)
  axis(1,lwd=2)
  axis(2,lwd=2)  
  box(lwd=2)
  for(i in 2:nrep){
    lines(1:periods,final_tot_n[r,i,],type='l',col=i)
    axis(1,lwd=2)
    axis(2,lwd=2)  
    box(lwd=2)
  }
}

