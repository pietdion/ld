#test "marginal correlation" for various copulas

#setwd("C:\Users\sgarcwh\Documents\PhD\marginal dependence")

#normal copula
norm <- normalCopula(0.617,dim=2)
u <- rCopula(100000, norm)

#gumbel copula
gum <- gumbelCopula(1.765,dim=2)
u <- rCopula(100000, gum)

#frank copula
fra <- frankCopula(4.45,dim=2)
u <- rCopula(100000, fra)

#joe copula
joe <- joeCopula(param=2.4, dim = 2)
u <- rCopula(100000, joe )

#t copula
t <- tCopula(0.68,df=1,dim=2)
u <- rCopula(100000, t)

#clayton copula
clay <- claytonCopula(1.51,dim=2)
u <- rCopula(100000, clay)

#structural copulas 1
n=100000; s=3.025;
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f=rgamma(n,1/2,1);
x=s*f+e1; y=s*f+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);

#structural copulas 2
n=100000; s1=1.05; s2=1.05;
e1=rnorm(n,0,1); e2=rnorm(n,0,1);
f1=rgamma(n,1,1);
f2=rgamma(n,1,1);
x=s1*f1-s2*f2+e1; y=s1*f1-s2*f2+e2;
u1=rank(x)/n; u2=rank(y)/n;
u=cbind(u1,u2);

#structural copulas 3
n=50000; 
u1=0; u2=0;
c=rbeta(n,2,2);

for (i in 1:n)
{
z=runif(1,0,1);
u1[i]=(z<=c[i])*runif(1,0,c[i]) + (z>c[i])*runif(1,c[i],1);
u2[i]=(z<=c[i])*runif(1,0,c[i]) + (z>c[i])*runif(1,c[i],1);
}
u=cbind(u1,u2);
rho=cor(u1,u2);

#structural copulas 4
n=30000; 
u1=0; u2=0;

for (i in 1:n)
{
z=runif(1,0,1); s=runif(1,0,1); d=0.57;
c=(s<=d)*runif(1,0,1)+(s>d)*0.5;
u1[i]=(z<=c)*runif(1,0,c) + (z>c)*runif(1,c,1);
u2[i]=(z<=c)*runif(1,0,c) + (z>c)*runif(1,c,1);
}
u=cbind(u1,u2);
rho=cor(u1,u2);
plot(u,cex=0.01);

####################




#compute marginal correlation
data=u; u1=data[,1]; u2=data[,2]; rho=cor(u1,u2);
n=100
cutoff=matrix((1:n)/n);
marcor=0; #marginal correlation


for (i in 1:n)
{
k=cutoff[i]; ind1=(u1<=k); 
marcor[i]=cov(u2,ind1)/cov(u1,ind1);
}

setEPS();
postscript("structural4.eps")
par(mar=c(5,5,1,1));
plot(u1[1:1000],u2[1:1000],type="p",cex=1,
xlab=expression(paste(list(u))),ylab=expression(paste(list(v))),cex.lab=2);
lines(cutoff,marcor,lty=1,lwd=5,col="red");
#lines(cutoff,rho*rep(1,n),lty=2,lwd=5);
#lines(cutoff,rep(0.65,n),lty=2,lwd=5,col="green");
#lines(rep(0.65,n),cutoff,lty=2,lwd=5,col="green");
dev.off();



setEPS();
postscript("structural4vs.eps")
par(mar=c(5,5.5,1,1));
plot((cutoff*gmarcor+1)/2,(cutoff*marcor+1)/2,type="l",lwd=4, xlim=c(0.5,1), ylim=c(0.5,1),
xlab="E(v|u>a) for Gaussian", ylab="E(v|u>a)",cex.lab=2,col="red");
lines(cutoff,cutoff,lty=2,lwd=5);
dev.off();











#gmarcor=marcor;

setEPS();
postscript("structural4abs.eps")
par(mar=c(5,5.5,1,1));
plot(cutoff,marcor-gmarcor,type="l",lwd=4, ylim=c(-0.4,0.4),
xlab=expression(alpha), ylab=expression(rho[alpha]-rho[alpha]^G),cex.lab=2,col="red");
lines(cutoff,rep(0,n),lty=2,lwd=5);
dev.off();

