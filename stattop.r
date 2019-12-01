

library(TDA)
library(RandomFields)
library(mvtnorm)
library(mgcv)
library(plot3D)



#####################################################
setrealdat=function(usesquare=FALSE,trend=3){
#This function reads the GASS data and sets up three R objects, one per region.
#Each object contains the data, empirical correlations, parameter estimates for fitted Matern correlations, persistence diagrams, and locations of local minima and maxima.  The correlations can be short or long range.
#
#The arguments are:
#  usesquare (TRUE/FALSE)  - to select square (TRUE) or cross (FALSE) neighbourhoods to define connections between components;
# trend (1-5) - to select one of fivetrend removal options.
#
#First read the GASS data, split into the three regions, detrend and marginally transform to N(0,1):
obs=matrix(scan("gass.txt"),256,1846,byrow=T)
use1=obs[,1:256]
use2=obs[,796:1051]
use3=obs[,1591:1846]
use1=detrend(use1,trend=trend)
use2=detrend(use2,trend=trend)
use3=detrend(use3,trend=trend)
use1=margtrans(use1)
use2=margtrans(use2)
use3=margtrans(use3)
#Rearrange so matrices match the images
ii=1:256
use1=t(use1[rev(ii),])
use2=t(use2[rev(ii),])
use3=t(use3[rev(ii),])
d=256 # all regions are 256*256
#Loop through the 3 regions:
for(it in 1:3){
fmat=use1
if(it==2) fmat=use2
if(it==3) fmat=use3
#Standardise:
fmat=(fmat-mean(as.vector(fmat)))/sd(as.vector(fmat))
#Find local minima and (by change of sign) local maxima:
fmin=getlocmin(fmat,usesquare=usesquare)
nmin=sum(as.vector(fmin))
ftem=-fmat
fmax=getlocmin(ftem,usesquare=usesquare)
nmax=sum(as.vector(fmax))
#Get empirical correlations using own function:
emp=getempcor1(fmat)
centres=emp$centres
empcor=emp$empcor
nbin=emp$nbin
#Fit a Matern model to these correlations:
matpar=fitmatcor(centres,empcor,nbin)
#Now get the empirical correlations over short distances only, and again fit Matern:
tem=getempcorloc(fmat,delta0=2)
matparloc=fitmatcor(tem$centres, tem$empcor,tem$nbin)
#Get the persistence diagram using a function from the TDA package:
pers=gridDiag(FUNvalues=fmat)
#Make permanent copies:
tem=list(fmat=fmat,centres=centres,empcor=empcor,pers=pers$diagram,nbin=nbin,fmin=fmin,nmin=nmin,fmax=fmax,nmax=nmax,matpar=matpar,matparloc=matparloc)
if(it==1) obs1<<-tem 
if(it==2) obs2<<-tem
if(it==3) obs3<<-tem
}
}
#############################################################
detrend=function(fmat,trend=1){
#Detrend an image using one of five methods (see below).
#
#The arguments are:
# fmat - an arbitrary square matrix;
# trend (1-5) - choice of detrending method. 
d=dim(fmat)[1]
#Vectorise
ff=as.vector(fmat)
#The x's are for polynomial fits:
x1=rep(1:d,each=d)
x2=rep(1:d,length=(d^2))
x3=x1*x2
x4=x1^2
x5=x2^2
x6=x1^3
x7=x2^3
x8=x1*x2^2
x9=x2*x1^2
x10=x1^4
x11=x2^4
x12=x1^3*x2
x13=x1^2*x2^2
x14=x1*x2^3
xmat1=cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14)
x15=x1^5
x16=x2^5
x17=x1^4*x2
x18=x1^3*x2^2
x19=x1^2*x2^3
x20=x1*x2^4
x21=x1^6
x22=x2^6
x23=x1^5*x2
x24=x1^4*x2^2
x25=x1^3*x2^3
x26=x1^2*x2^4
x27=x1*x2^5
x28=x1^7
x29=x2^7
x30=x1^6*x2
x31=x1^5*x2^2
x32=x1^4*x2^3
x33=x1^3*x2^4
x34=x1^2*x2^5
x35=x1*x2^6
xmat2=cbind(x15,x16,x17,x18,x19,x20,x21,x22,x23,x24,x25,x26,x27)
xmat2=cbind(xmat2,x28,x29,x30,x31,x32,x33,x34,x35)
if(trend==1){
#Polynomial trend removal, low order (15 terms, up to order 4):
fit=lm(ff~xmat1)
res=matrix(fit$residuals,d,d)
}
if(trend==2){
#Polynomial trend removal, higher order (36 terms, up to order 7):
fit=lm(ff~xmat1+xmat2)
res=matrix(fit$residuals,d,d)
}
if(trend==3){
#Thin plate smoothing spline using gam function from mgcv package:
fit=gam(ff~s(x1,x2))
res=matrix(fit$residuals,d,d)
}
if(trend==4){
#Thin plate smoothing splne with extra terms:
fit=gam(ff~s(x1,x2,m=10))
res=matrix(fit$residuals,d,d)
}
if(trend==5){res=fmat} #No trend removal
res
}

#####################################################
margtrans=function(fmat){
#Vectorise an image, empirically marginally transform to N(0,1),then re-format.
#
#The argument fmat is an arbitrary square matrix.
d=dim(fmat)[1]
ff=as.vector(fmat)
n=length(ff)+1
ii=rank(ff)
pp=ii/(n+1)
z=qnorm(pp)
matrix(z,d,d)
}
#################################################################
getlocmin=function(fmat, usesquare=FALSE){
#This function locates the local minima in a matrix.
#
#The arguments are:
# fmat - an arbitrary square matrix;
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE).
d=dim(fmat)[1] 
fvec=as.vector(fmat)
fmax=max(as.vector(fmat))+10
#Embed the matrix in a larger matrix with boundaries larger than any original value.  This means we can find the local minima by fmat with shifted versions of itself, with shifts in the horizontal, vertical and (for usesquare==TRUE) diagonal directions.
big=rbind(fmax,fmat,fmax)
big=cbind(fmax,big,fmax)
ii=2:(d+1)
n=d^2
fmin=matrix(TRUE,d,d)
xa=c(1,-1,0,0)
ya=c(0,0,1,-1)
#Compare the elements of fmat with their shared-edge neighbours, one edge at a time
for(it in 1:4){
ftem=big[ii+xa[it],ii+ya[it]]
fmin=fmin&(fvec<ftem)
}
if(usesquare==TRUE){
xa=c(1,-1,1,-1)
ya=c(1,1,-1,-1)
#Repeat for shared vertices if square neighbourhood selected
for(it in 1:4){
ftem=big[ii+xa[it],ii+ya[it]]
fmin=fmin&(fvec<ftem)
}
}
#This produces the locations of local minima
fmin
}

#########################################################
getempcor2=function(fmat,delta0=3){
#Get isotropic empirical correlations.
#It is convenient to set this up as a matrix in the same form as the anisotropic version rather than as a vector.
#
#The arguments are:
# fmat - an arbitrary square matrix;
# delta0 - a positive integer defining the maximum local distance to use.
d=dim(fmat)[1]
R1=matrix(0,2*delta0+1,delta0+1)
#This looks at all pairs of points separated by lags up to delta0: 
for(it in 0:delta0){
for(jt in 0:delta0){
x1=1:(d-it)
y1=1:(d-jt)
x2=(it+1):d
y2=(jt+1):d
x3=(it+1):d
y3=1:(d-jt)
x4=1:(d-it)
y4=(jt+1):d
v1=as.vector(fmat[x1,y1])
v2=as.vector(fmat[x2,y2])
v3=as.vector(fmat[x3,y3])
v4=as.vector(fmat[x4,y4])
R1[it+delta0+1,jt+1]=cor(v1,v2)
if(it>0) R1[delta0+1-it,jt+1]=cor(v3,v4)
}}
#This makes the estimates isotropic:
R=R1[-(1:(delta0)),]
R2=R1[(delta0:1),]
R[2:(delta0+1),]=0.5*(R[2:(delta0+1),]+R2)
R1=0.5*(R+t(R))
R=rbind(R1[(delta0+1):2,],R1)
R
}
##########################################################
getempcorloc=function(fmat,delta0=2,nu=0.5,eta=10){
#Gets all empirical correlations in local neighbourhood assuming isotropy
#Also - for reference - gets Matern values for the same distances.
#
#The arguments are:
# fmat - an arbitrary square matrix;
# delta0 - a positive integer defining the maximum local distance to use;
# nu,eta - Matern shape and scale parameters.
#
#Get local correlations to distance delta0, using own function:
R1=getempcor2(fmat=fmat,delta0=2)
R1=R1[-(1:delta0),]
d=r=NULL
#Find distances and correlations:
for(it in 1:(delta0+1)){
for(jt in it:(delta0+1)){
d=c(d,sqrt((it-1)^2+(jt-1)^2))
r=c(r,R1[it,jt])
}}
ii=order(d)
d=d[ii][-1]
r=r[ii][-1]
#Give equal weight
nbin=d*0+1
#Obtain Matern values
z=sqrt(2*nu)*d/eta
cc=2^(1-nu)*z^nu*besselK(z,nu)/gamma(nu)
comp=cbind(d,r,cc)
dimnames(comp)=list(NULL,c("Dist","Empirical","Matern"))
list(centres=d,empcor=r,nbin=nbin,comp=comp)
}
#########################################################
getempcor1=function(fmat,lags=c(1:10,12,14,16,18,10*(2:7))){
#Gets empirical correlations in horizontal and vertical directions (combined).
#
#The arguments are:
# fmat - an arbitrary square matrix;
# lags - a vector of positive integers giving the distances at which to estimate correlations.
d=dim(fmat)[1]
#Check that the distances do not exceed the matrix dimension:
lagsuse=lags[lags<(d-1)]
m=length(lagsuse)
cc=nn=NULL
#For each distance find pairs of points separated by that amount either vertically or horizontally:
for(it in 1:m){
dd=lagsuse[it]
x=1:(d-dd)
y=1:(d-dd)
f1=as.vector(fmat[x,y])
f3=as.vector(fmat[x+dd,y])
f4=as.vector(fmat[x,y+dd])
fa=c(f1,f1)
fb=c(f3,f4)
nn=c(nn,length(c(f3,f4)))
cc=c(cc,cor(fa,fb))
}
#For convenience us same names as RandomFields package
list(centres=lagsuse,empcor=cc,nbin=nn)
}

#####################################################
fitmatcor=function(centres,empcor,nbin){
#Fits a Matern correlation function to empirical values.
#
#The arguments are:
# centres - a non-negative vector of distances; 
# empcor  - the correlations at these distances;
# nbin    - a weighting to be applied.
#
#Initial values:  
nu=0.5
eta=20
p0=c(nu,eta)
#Estimate using the R function nlm to minimise the weighted sum of squares between observed and fitted correlations:
fit=nlm(getcorss,p=p0,centres=centres,empcor=empcor,nbin=nbin,iterlim=1000)
fit$estimate
}
###########################################################
getcorss=function(p,centres=centres,empcor=empcor,nbin=nbin){
#  Function to obtain the weighted sum of squares between observed and fitted correlations.
#Arguments are:
# p - a 2-vector of parameter estimates;
# centres - a non-negative vector of distances; 
# empcor  - the correlations at these distances;
# nbin    - a weighting to be applied.
nu=p[1]
eta=p[2]
#Return a very large value if negative estimates attempted:
if(min(p)<=0) ss=10000000 else{
ii=nbin>0
w=nbin/sum(nbin)
#Get Matern values:
z=sqrt(2*nu)*centres/eta
cc=2^(1-nu)*z^nu*besselK(z,nu)/gamma(nu)
cc[centres==0]=1
fit=cc
#Get weighted sum of squares:
ss=sum(w[ii]*(fit[ii]-empcor[ii])^2)
}
ss
}
#############################################################
getbetplot=function(levs=c(-2,-1,0,1),useplot=FALSE){
#This produces the simple example for Section 2.1 of the paper.
#
#The arguments are
# levs - the levels for four level sets;
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.
dd=10 #Size of grid
ss=117 #Seed
lbs=c("(a) ", "(b) ","(c) ","(d) ")# Plot labels
#Setting a seed so can reproduce a good illustration
set.seed(ss)
#The simulation is a Gaussian random field with Matern correlation.
#It is generated using the RandomField package.
model=RMmatern(nu=0.5,var=1,scale=2)
f=RFsimulate(model,x=1:dd,y=1:dd)$variable1
#Some processing:
dat=as.vector(round(f,2))
dat=(dat-mean(dat))/sd(dat)
n=length(dat)
xleft=rep(0:(dd-1),length=n)
ybottom=rep((dd-1):0,each=dd)
xright=xleft+1
ytop=ybottom+1
z=sort(as.vector(dat))
#useplot is used to create a postscript file.
if(useplot==TRUE) postscript("gridfig.ps",horizontal=FALSE)
#The four level sets:
par(mfrow=c(2,2))
for(it in 1:4){
lev=levs[it]
ll=paste(lbs[it],"t=",as.character(lev),sep="")
plot(0:dd,0:dd,pch="",xaxt='n',yaxt='n',xlab="",ylab="",bty='n',main=ll,cex.main=2)
for(jt in 0:dd){
lines(c(0,dd),c(jt,jt))
lines(c(jt,jt),c(0,dd))
}
ii=(1:n)[dat<=lev]
rect(xleft[ii],ybottom[ii],xright[ii],ytop[ii],col=1)
}
par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}
############################################################
getexampleplot=function(useplot=FALSE, use3d=TRUE){
#This is a routine to visualise the simple example data.
#
#The arguments are: 
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not;
# use3d - whether to draw a 3d histogram (TRUE) or a 2D image (FALSE).
if(use3d==TRUE) fnm="examplefig1.ps" else fnm="examplefig2.ps"
#Re-generate the data (instead of passing it):
dd=10
ss=117
set.seed(ss)
model=RMmatern(nu=0.5,var=1,scale=2)
f=RFsimulate(model,x=1:dd,y=1:dd)$variable1
dat=as.vector(round(f,2))
dat=(dat-mean(dat))/sd(dat)
dat=t(matrix(dat,dd,dd,byrow=FALSE))

if(useplot==TRUE) postscript(fnm,horizontal=FALSE)
if(use3d==TRUE){
#hist3D is from package plot3D:
hist3D(x=1:10,y=1:10,z=dat,scale = FALSE, expand = 1, bty = "g", phi = 20,
        col = "#0072B2", border = "black", shade = 0.2, ltheta = 90,
        space = 0.3, ticktype = "detailed", d = 2,ylab="x",xlab="y")}else{
dat=dat[dd:1,]
image(t(dat),col=grey(seq(0,1,length=100)),xlab="",ylab="",xaxt='n',yaxt='n')}
if(useplot==TRUE) dev.off()
}
#############################################################
getpersfig=function(useplot=FALSE){
#This routine plots the persistence diagram and convex peels for the region 1 GASS data.
#
#The argument is: 
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.
if(useplot==TRUE) postscript("persfig.ps", horizontal=TRUE)
par(mfrow=c(1,2))
#Plot the persistence diagram:
persplot1(obs1$pers)
#Plot the convex peels, for components, 90% thinning.
tem=convpeelplot1(obs1$pers,usetype=1,prop=0.9,draw=TRUE)
par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}
###################################################################
getallpers=function(useplot=FALSE){
#This routine plots the persistence diagrams for all three GASS regions for both components and holes.
#
#The argument is: 
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.
if(useplot==TRUE) postscript("persall.ps", horizontal=FALSE)
par(mfrow=c(3,2))
#mn is the subplot title.
#usetype chooses components (1,default) or holes (2):
persplot1(obs1$pers,mn="Region 1, components")
persplot1(obs1$pers,usetype=2, mn="Region 1, holes")
persplot1(obs2$pers,mn="Region 2, components")
persplot1(obs2$pers,usetype=2, mn="Region 2, holes")
persplot1(obs3$pers,mn="Region 3, components")
persplot1(obs3$pers,usetype=2, mn="Region 3, holes")

par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}
###################################################################
getallconvpeels=function(useplot=FALSE){
#This routine plots the convex peels for all three GASS regions for both components and holes.
#
#The argument is: 
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.
if(useplot==TRUE) postscript("allpeels.ps", horizontal=FALSE)
par(mfrow=c(3,2))
#mn is the subplot title.
#usetype chooses components (1,default) or holes (2).
#prop is the proportion of thinning.
#draw indicates that the peels should be drawn.
tem=convpeelplot1(obs1$pers,mn="Region 1, components",prop=0.9,draw=TRUE)
tem=convpeelplot1(obs1$pers,usetype=2, mn="Region 1, holes",prop=0.9,draw=TRUE)
tem=convpeelplot1(obs2$pers,mn="Region 2, components",prop=0.9,draw=TRUE)
tem=convpeelplot1(obs2$pers,usetype=2, mn="Region 2, holes",prop=0.9,draw=TRUE)
tem=convpeelplot1(obs3$pers,mn="Region 3, components",prop=0.9,draw=TRUE)
tem=convpeelplot1(obs3$pers,usetype=2, mn="Region 3, holes",prop=0.9,draw=TRUE)

par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}

###############################################################
persplot1=function(pers,usetype=1,mn=""){
#Routine to plot a persistence diagram.
#
#The arguments are:
# pers - a n*3 matrix giving feature type and x and y coordinates;
# usetype (1/2) - plot components (1) or holes (2);
# mn - a plot header.
#The feature types  produced by the TDA package are 0 and 1.
#In other routines they are coded 1 and 2.
x=pers[,2]
y=pers[,3]
typ=pers[,1]+1#Using 1,2 rather than 0,1
yl=range(c(x,y))
x1=x[typ==usetype]
y1=y[typ==usetype]
cl=c(2,4)
plot(x1,y1,pch=20,col=cl[usetype],xlim=yl,ylim=yl,xlab="Birth level",ylab="Death level",main=mn)
lines(yl,yl)
}
###########################################################################
convpeelplot1=function(pers,usetype=1,prop=0.95,draw=TRUE,mn=""){
#Routine to plot convex peels of a persistence diagram.
#
#The arguments are:
# pers - a n*3 matrix giving feature type and x and y coordinates;
# usetype (1/2) - plot components (1) or holes (2);
# prop - the proportion of points to peel away;
# draw (TRUE/FALSE) - an indicator of whether plots should be produced, or the peels simply returned.
#The feature types  produced by the TDA package are 0 and 1.
#In other routines they are coded 1 and 2.

x=pers[,2]
y=pers[,3]
typ=pers[,1]+1#Using 1,2 rather than 0,1
yl=range(c(x,y))
x1=x[typ==usetype]
y1=y[typ==usetype]
if(draw==TRUE) {plot(x1,y1,xlim=yl,ylim=yl,xlab="Birth level",ylab="Death level",pch=" ")
lines(yl,yl)}
n=length(x1)
nleft=n
#The current data:
dat1=cbind(x1,y1)
#Colors for components and holes:
cl=c(2,4)
#Successively peel until only proportion prop left:
while(nleft>(prop*n)){
nleft=length(dat1[,1])
#Get the convex hull of the current data:
ch=getconvhull(dat1)
if(draw==TRUE) lines(ch$xch,ch$ych,col=grey(0.7))
#Peel:
dat1=dat1[-ch$ch[-1],]
}
if(draw==TRUE) points(ch$xch,ch$ych,col=cl[usetype],pch=20,type="b")
#Return the final convex hull:
ch
}
###########################################
getconvhull=function(dat1,draw=FALSE){
#Get the convex hull of a set of points.
#
#The arguments are:
# dat1 - an arbitrary n*2 set of locations;
# draw - an indicator of whether to add the convex hull to an open plot.
#
#Get the locations:
x=dat1[,1]
y=dat1[,2]
n=length(x)
if(draw==TRUE) plot(x,y)
#Find the smallest x and set initial vertex:
c1=(1:n)[x==min(x)]
if(length(c1>1)) c1=c1[1]
chlow=c1  
if(draw==TRUE) points(x[c1],y[c1],pch=20,col=2)
#Iterate round the low side of the convex hull, working from left to right:
for(it in 1:n){
#Find the angles to the other points to right of current vertex, and the smallest angle:
lst=(1:n)[x>x[c1]]
dd=sqrt((x-x[c1])^2+(y-y[c1])^2)
th=asin((y[lst]-y[c1])/dd[lst])
c2=lst[th==min(th)]
if(length(c2>1)) c2=c2[1]
#Add this point to convex hull and re-set the vertex:
chlow=c(chlow,c2)
c1=c2
#Stop when there are no further points to the right:
if(x[c1]==max(x)) break()
}
#Repeat for the high side of the convex hull, working from right to left:
chhigh=c1
for(it in 1:n){
lst=(1:n)[x<x[c1]]
dd=sqrt((x-x[c1])^2+(y-y[c1])^2)
th=asin((y[lst]-y[c1])/dd[lst])
c2=lst[th==max(th)]
if(length(c2>1)) c2=c2[1]
chhigh=c(chhigh,c2)
c1=c2
if(x[c1]==min(x)) break()
}
if(draw==TRUE) {lines(x[chlow],y[chlow],type="b",pch=20,col=2)
lines(x[chhigh],y[chhigh],type="b",pch=20,col=2)}
#The convex hull as indices in the original data:
ch=c(chlow,chhigh)
#The x,y coordinates of the convex hull:
xch=c(x[chlow],x[chhigh])
ych=c(y[chlow],y[chhigh])
#Return these, together with the low and high sides separately:
list(xlow=x[chlow],ylow=y[chlow],xhigh=x[chhigh],yhigh=y[chhigh],ch=ch,xch=xch,ych=ych)
}

##########################################################################
getenmax=function(d=256,corrtype='matern',eta=20,nbr='cross',useedge=TRUE,nu=0.5,empirical=FALSE,anis=FALSE,fmat=obs1$fmat){
#This routine calculates the expected number of local maxima for a Gaussian random field and an approximate standard deviation.
#
#The arguments are:
# d - the size of the GRF, assumed to be a d*d matrix;
# corrtype - either 'matern' or 'exp' for Matern or exponential;
# nu, eta - Matern shape and scale parameters;
# nbr ('cross' or 'square') - neighbourhood in defining connections as edge-only ('cross') or edge or vertex ('square'); 
# useedge - whether or not to adjust for edge effects near boundary of matrix;
# empirical (TRUE/FALSE) - whether to use the empirical correlations or Matern;
# anis (TRUE/FALSE) - whether to assume isotropy (anis=FALSE) or anisotropy (anis=TRUE) if empirical correlations selected;
# fmat - an arbitrary square matrix from which to take empirical correlations if that option is selected.
#
# R will be the correlation matrix for the conditional distributions. We only need local values:
if(empirical==TRUE) {R=getcorns(d=8,fmat=fmat, anis=anis)} else{
 R=getcor(d=8,eta=eta,nu=nu,corrtype=corrtype)}
#Pick out covariances in different directions:
c01=getcovi1i2(x11=c(2,2),x12=c(2,3),type=nbr,R=R)
c11=getcovi1i2(x11=c(2,2),x12=c(3,3),type=nbr,R=R)
c02=getcovi1i2(x11=c(2,2),x12=c(2,4),type=nbr,R=R)
c12=getcovi1i2(x11=c(2,2),x12=c(3,4),type=nbr,R=R)
c22=getcovi1i2(x11=c(2,2),x12=c(4,4),type=nbr,R=R)
#How many points of each type (n1=corners, n2=other edges, n3=other): 
n1=4
n2=4*(d-2)
n3=(d-2)^2
if(useedge==FALSE) n3=d^2
#The different possible configurations:
if(nbr=='cross'){
l1='cornercross'
l2='edgecross'
l3='cross'}else{
l1='cornersquare'
l2='edgesquare'
l3='square'
}
#Get the probability of a local maximum for each configuration, given the correlations:
e3=getei1edge(R=R,nbr=l3)
e1=getei1edge(R=R,nbr=l1)
e2=getei1edge(R=R,nbr=l2)
#Multiply by the appropriate number of points and sum, with variance approximation:
e=n3*e3
v=n3*e3*(1-e3)
if(useedge==TRUE){
e=e+n1*e1
v=v+n1*e1*(1-e1)
e=e+n2*e2
v=v+n2*e2*(1-e2)}
#Edge effects negligible for covariance correction.
n33=(d-4)^2
c12=c22=0
v=v+4*n33*(c01+c11+c02+2*c12+c22)
c(e,sqrt(v))
}
#########################################################################
getcorns=function(d=8,fmat=obs1$fmat, anis=FALSE){
#Routine to get empirical correlations.
#
#Arguments are
# d - positive integer giving maximum separation of interest;
# fmat - arbitrary square matrix from which to take correlations;
# anis - assume isotropy (anis=FALSE) or not (anis=TRUE).
if(anis==FALSE)  R1=getempcor2(fmat,delta0=(d-1)) else R1=getempcor3(fmat,delta0=(d-1))
#Look at all possible distances:
ii=1:d
i1=j1=ipos=NULL
kt=0
for(jt in ii){
for(it in ii){
kt=kt+1
i1=c(i1,it)
j1=c(j1,jt)
ipos=c(ipos,kt)
}
}
n=d^2
n2=n^2
dvec=rep(0,length=n2)
kt=0
for(it in 1:n){
for(jt in 1:n){
kt=kt+1
x1=i1[ipos[it]]
x2=i1[ipos[jt]]
y1=j1[ipos[it]]
y2=j1[ipos[jt]]
if((x2>=x1)&(y2>=y1)){del1=d+x2-x1
del2=y2-y1+1
}
if((x2>=x1)&(y2<y1)){del1=d+x1-x2
del2=y1-y2+1
}
if((x2<x1)&(y2>=y1)){del1=d+x2-x1
del2=y2-y1+1
}
if((x2<x1)&(y2<y1)){del1=d+x1-x2
del2=y1-y2+1
}
dvec[kt]=R1[del1,del2]
}
}
R=matrix(dvec,n,n)
R
}
###########################################################################
getcovi1i2=function(x11=c(2,2),x12=c(2,5),type='cross',R=getcor(d=d,corrtype="exp",eta=eta,nu=0.5)){
#Gets the covariance between the local max indicators at any two locations. 
#
#Arguments are:
# x11, x12 - the two locations;  The first location should not be on a boundary.
# type ('cross' or 'square') - neighbourhood in defining connections as edge-only ('cross') or edge or vertex ('square'); 
# R - the local correlation matrix to use.
#
#The probability of a maximum at each location, and then simultaneously:
e1=getei1(x1=x11,type=type,R=R)
e2=getei1(x1=x12,type=type,R=R)
e12=getei1i2(x11=x11,x12=x12,type=type,R=R)
round((e12-e1*e2),5)
}

#########################################################################
getei1=function(x1=c(2,2),type='cross',R=getcor(d=8,corrtype="exp",eta=10,nu=0.5)){
#This routine finds the probability of a local max at any location given its local correlations.
#
#Arguments are:
# xl - the location;
# R  - the local correlation in prescribed form.  
d=sqrt(dim(R)[1])
# Get the position of xl relative to R:                
loc1=getpos(x1[1],x1[2],d=d)
#Set up conditional distributions:
m2=epoint(i1=x1[1],j1=x1[2],type=type,R=R)
loc2=m2[,3]
R22=R[loc2,loc2]
r=R[loc2,loc1]
zero=loc2*0
one=zero+1
V=R22+one%*%t(one)-one%*%t(r)-r%*%t(one)
#Find the probability of a local maximum using the multivariate normal df function in mvtnorm package:
ploc=pmvnorm(upper=zero,mean=zero,sigma=V,algorithm=Miwa(steps=200))[1]
ploc
}
#############################################################
getpos=function(i1,j1,d){
#Original data are a grid. These get converted to column form to calculate R.
#This is to find where in the column a point is, or a set of points.  It is needed
#to find the relevant elements of R when calculating the mean and variance of the
#number of local max.
n=d^2
ii=jj=ipos=rep(0,length=n)
kt=0
for(jt in 1:d){
for(it in 1:d){
kt=kt+1
ii[kt]=it
jj[kt]=jt
ipos[kt]=kt
}
}
iout=NULL
for(it in 1:length(i1)){
iout=c(iout,ipos[(ii==i1[it])&(jj==j1[it])]  )
}
iout
}
###########################################################
epoint=function(i1=1,j1=1,type='cross',R=getcor(d=8,corrtype="exp",eta=10,nu=0.5)){
#This finds the cross or square neighbours of i1,j1 and their locations in the vector form used for R.
#Arguments are:
# i1, j1 - the two locations in the vector form;  The first location should not be on a boundary.
# type ('cross' or 'square') - neighbourhood in defining connections as edge-only ('cross') or edge or vertex ('square'); 
# R - the local correlation matrix to use.
#
#Obs are numbered 1:n where n=d^2
d=sqrt(dim(R)[1])
ii=1:d
n=d^2
x=rep(ii,length=n)
y=rep(ii,length=n,each=d)
ipos1=(1:n)[(x==i1)&(y==j1)]
i2=c(i1-1,i1,i1,i1+1)
j2=c(j1,j1-1,j1+1,j1)
if(type=='square'){
i2=c(i2,i1-1,i1+1,i1-1,i1+1)
j2=c(j2,j1-1,j1-1,j1+1,j1+1)
}
iuse=(i2>0)&(i2<=d)&(j2>0)&(j2<=d)
i2=i2[iuse]
j2=j2[iuse]
ipos2=NULL
for(it in 1:length(i2)){
ipos2=c(ipos2,(1:n)[(x==i2[it])&(y==j2[it])])
}
cbind(i2,j2,ipos2)
}
####################################################
getcor=function(d=32,corrtype="exp",eta=10,nu=0.5){
#This routine creates a correlation matrix.
#
#The arguments are
# d - the vector length of interest
# corrtype ('exp', 'matern') - exponential or Matern correlations;
# eta - scale parameter;
# nu - shape parameter if Matern selected.
if(eta==0){
n=d^2
R=matrix(0,n,n)
diag(R)=1
}else{
ii=1:d
i1=j1=ipos=NULL
kt=0
for(jt in ii){
for(it in ii){
kt=kt+1
i1=c(i1,it)
j1=c(j1,jt)
ipos=c(ipos,kt)
}
}
n=d^2
n2=n^2
#Sets up all pairwise distances:
dvec=dpos1=dpos2=rep(0,length=n2)
kt=0
for(it in 1:n){
for(jt in 1:n){
kt=kt+1
dpos1[kt]=ipos[it]
dpos2[kt]=ipos[jt]
x1=i1[ipos[it]]
x2=i1[ipos[jt]]
y1=j1[ipos[it]]
y2=j1[ipos[jt]]
dvec[kt]=sqrt((x1-x2)^2+(y1-y2)^2)
}
}
if(corrtype=='exp'){
dmat=matrix(dvec,n,n)
R=exp(-dmat/eta)
}
if(corrtype=='matern'){
z=sqrt(2*nu)*dvec/eta
cc=2^(1-nu)*z^nu*besselK(z,nu)/gamma(nu)
R=matrix(cc,n,n)
diag(R)=1
}
}
R
}
################################################################################
getei1i2=function(x11=c(2,2),x12=c(2,5),type='cross',R=getcor(d=8,corrtype="exp",eta=10,nu=0.5)){
#Finds E[I1*I2] for two locations x11 and x12 with no overlapping neighbours.
#
#Arguments are:
# x11, x12 - the two locations;  The first location should not be on a boundary.
# type ('cross' or 'square') - neighbourhood in defining connections as edge-only ('cross') or edge or vertex ('square'); 
# R - the local correlation matrix to use.
#

d=sqrt(dim(R)[1])
#Get R11
i11=getpos(x11[1],x11[2],d=d)
i12=getpos(x12[1],x12[2],d=d)
loc1=c(i11,i12)
R11=R[loc1,loc1]
#Get R22
m21=epoint(i1=x11[1],j1=x11[2],type=type,R=R)
m22=epoint(i1=x12[1],j1=x12[2],type=type,R=R)
k1=dim(m21)[1]
k2=dim(m22)[1]
k=k1+k2
loc2=c(m21[,3],m22[,3])
R22=R[loc2,loc2]
#Get R12
R12=R[loc2,loc1]
#Get E[i1I2]
J=matrix(0,k,2)
J[1:k1,1]=1
J[k1+(1:k2),2]=1
R11inv=solve(R11)
Rtem=R12%*%R11inv
D=J-Rtem
V=R22-Rtem%*%t(R12)+D%*%R11%*%t(D)
zero=rep(0,length=k)
#The Miwa algorithm is repeatable and supposedly works up to dim 20 only, so
#use GenzBretz for products.
ploc=pmvnorm(upper=zero,mean=zero,sigma=V,algorithm=GenzBretz(abseps=0.0000001))[1]
ploc
}
###############################################################
getei1edge=function(x=c(2,2),nbr='cross',R=getcor(d=8,corrtype="exp",eta=10,nu=0.5)){
#Allows edge effects
d=sqrt(dim(R)[1])
ii=c(0,0,0,-1,-1,-1,1,1,1)
jj=c(0,-1,1,0,-1,1,0,-1,1)
if(nbr=='cross') use=c(1,1,1,1,0,0,1,0,0)
if(nbr=='square') use=c(1,1,1,1,1,1,1,1,1)
if(nbr=='edgecross') use=c(1,1,0,1,0,0,1,0,0)
if(nbr=='edgesquare') use=c(1,1,0,1,1,0,1,1,0)
if(nbr=='cornercross') use=c(1,1,0,1,0,0,0,0,0)
if(nbr=='cornersquare') use=c(1,1,0,1,1,0,0,0,0)
i=ii[use==1]
j=jj[use==1]
p=length(i)
locs=i*0
for(it in 1:p){
locs[it]=getpos(x[1]+i[it],x[2]+j[it],d=d)
}
R1=R[locs,locs]
R22=R1[2:p,2:p]
r=R1[2:p,1]
one=matrix(1,p-1,1)
V=R22+one%*%t(one)-one%*%t(r)-r%*%t(one)
zero=as.vector(one*0)
ploc=pmvnorm(upper=zero,mean=zero,sigma=V,algorithm=Miwa(steps=200))[1]
ploc
}
######################################################
ensims=function(nsim=10,dist='Gauss',d=256,usesquare=FALSE){
#Routine to find distribution of number of local maxima and minima via simulations.
#
#Arguments are:
# nsim - the number of simulations;
# dist - the distribution to use, selected from 'Gauss', 'chisq1', 'chisq3', 'T', 'F';
# d - we assume a d*d field
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE).


if(usesquare==TRUE) nbr='Square' else nbr='Cross'

nmin=nmax=rep(0,length=nsim)
for(jt in 1:nsim){
cat("d=",d," Nbr=",nbr," Sim=",jt,"\n")
dat=simdat(d=d,dist=dist,usesquare=usesquare)
nmin[jt]=dat$nmin
nmax[jt]=dat$nmax
}
#Report mean and standard deviation of numbers of local max or min:
out1=c(mean(nmin),sd(nmin),sd(nmin)/sqrt(nsim))
out2=c(mean(nmax),sd(nmax),sd(nmax)/sqrt(nsim))
list(nmax=out2,nmin=out1)
}


#####################################################
plotsims=function(useplot=FALSE,useseed=TRUE){
#This routine illustrates fields simulated under different models.
#
#The arguments are:
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not;
# useseed (TRUE/FALSE) - an indicator of whether or not to use a seed for repeatability.

if(useseed==TRUE) set.seed(50)
if(useplot==TRUE) postscript("simexamples.ps",horizontal=FALSE)
par(mfrow=c(3,2))
#Plot an image of one simulation from each of five distributions with default parameter values:
image(simdat(dist='Gauss',onlydat=TRUE)$fmat,main="Gaussian")
image(simdat(dist='chisq1',onlydat=TRUE)$fmat,main="Transformed Chisq(1)")
image(simdat(dist='chisq3',onlydat=TRUE)$fmat,main="Transformed Chisq(3)")
image(simdat(dist='T',onlydat=TRUE)$fmat,main="Transformed T(3)")
image(simdat(dist='F',onlydat=TRUE)$fmat,main="Transformed F(3,3)")
par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}

#####################################################
simdat=function(d=256,corrtype='matern',eta1=20,nu1=0.5,dist='Gauss',usesquare=FALSE,usedefault=TRUE,onlydat=FALSE){
#This is the main simulation routine.  It will generate a square matrix/field from one of five different distributions. The root (generating) correlation matrix is exponential or Matern, with either user-specified parameters or with default values chosen to give final (target) correlations close to exponential with scale 20. The final field is marginally transformed to Gauss.
#
#The arguments are:
# d - the size of the field, ie a d*d matrix;
# corrtype ('exp', 'matern') - exponential or Matern correlations for root correlation;
# eta1, nu1 - non-negative scale and shape parameters for root correlation;
# dist - the distribution to use, selected from 'Gauss', 'chisq1', 'chisq3', 'T', 'F';
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE).
# usedefault (TRUE/FALSE)- whether to use default parameter values (TRUE) or those specified by eta1 and nu1 as root correlations;
# onlydata (TRUE/FALSE) - whether only to generate the field data (TRUE) or also to calculate empirical correlations, maxima/minima locations and persistence diagrams (see below).


#Choose the parameters for the root correlation:
nu=nu1
eta=eta1
if(usedefault==TRUE){
if(dist=='Gauss') pp=c(0.5,20)
if(dist=='chisq1') pp=c(0.74,41)
if(dist=='chisq3') pp=c(0.54,42)
if(dist=='T') pp=c(0.58,22)
if(dist=='F') pp=c(0.54,50)
if(dist=='corrcheck') pp=c(0.65,12)
nu=pp[1]
eta=pp[2]
}

#There is an option for alternative (from default) degrees of freedom for chisq, T and F:
dfchisq1=1
dfchisq3=3
dfT=3
dfF=c(3,3)
#The total number of pixels:
n=d^2
#Setting the model up, with a routine from the RandomFields package:
model=RMmatern(nu=nu,var=1,scale=eta)
#Initialise
f=rep(0,length=n)
#Generate the data according to the distribution selected.  The routine RFsimulate is from RandomFields and generates a Gaussian random field:
if(dist=='Gauss'){
#Gauss is immediate:
if(eta==0) f=rnorm(n) else f=RFsimulate(model,x=1:d,y=1:d)$variable1
f=matrix(f,d,d)
}
if(dist=='chisq1'){
#chisq1 is the square of a Gauss:
for(jt in 1:dfchisq1){
#Using a separate method for exponential than Matern is a tiny bit quicker:
if(eta==0) z=rnorm(n) else z=RFsimulate(model,x=1:d,y=1:d)$variable1
f=f+z^2
}
#This performs the marginal transformation (similar below):
f=qnorm(pchisq(f,dfchisq1))
}
if(dist=='chisq3'){
#chisq3 is the sum of 3 independent squared Gauss:
for(jt in 1:dfchisq3){
if(eta==0) z=rnorm(n) else z=RFsimulate(model,x=1:d,y=1:d)$variable1
f=f+z^2
}
f=qnorm(pchisq(f,dfchisq3))
}
if(dist=='T'){
#T is from a Gauss divided by the square root of a scaled chisq: 
if(eta==0) znum=rnorm(n) else znum=RFsimulate(model,x=1:d,y=1:d)$variable1
zden=znum*0
for(jt in 1:dfT){
if(eta==0) z=rnorm(n) else z=RFsimulate(model,x=1:d,y=1:d)$variable1
zden=zden+z^2
}
f=znum*sqrt(dfT)/sqrt(zden)
f=qnorm(pt(f,dfT))
}
if(dist=='F'){
#F is the ratio of two chisq:
znum=zden=rep(0,length=n)
for(jt in 1:dfF[1]){
if(eta==0) z=rnorm(n) else z=RFsimulate(model,x=1:d,y=1:d)$variable1
znum=znum+z^2
}
for(jt in 1:dfF[2]){
if(eta==0) z=rnorm(n) else z=RFsimulate(model,x=1:d,y=1:d)$variable1
zden=zden+z^2
}
f=dfF[2]*znum/(dfF[1]*zden)
f=qnorm(pf(f,dfF[1],dfF[2]))
}
#Form a square matrix and standardise
fmat=matrix(f,d,d)
fmat=(fmat-mean(as.vector(fmat)))/sd(as.vector(fmat))
#The onlydat option is here because in some simulations we only want the field, in others we also want summaries (which slow down the generation):
if(onlydat==FALSE){
#Locations and numbers of local minima and maxima, using own routines:
fmin=getlocmin(fmat,usesquare=usesquare)
nmin=sum(as.vector(fmin))
ftem=-fmat
fmax=getlocmin(ftem,usesquare=usesquare)
nmax=sum(as.vector(fmax))
fmat=(fmat-mean(as.vector(fmat)))/sd(as.vector(fmat))
#The empirical correlations using own routine:
emp=getempcor1(fmat)
centres=emp$centres
empcor=emp$empcor
nbin=emp$nbin
#The persistence diagrams, using a routine from the TDA package:
pers=gridDiag(FUNvalues=fmat)$diagram}else{
#Set all these to be NULL if onlydata selected:
centres=empcor=pers=nbin=fmin=nmin=fmin=nmin=fmax=nmax=NULL}
list(fmat=fmat,centres=centres,empcor=empcor,pers=pers,nbin=nbin,fmin=fmin,nmin=nmin,fmax=fmax,nmax=nmax)
}

######################################################################################
gettab1=function(nugrid=c(0.5,0.5,1.5,1.5),etagrid=c(10,20,10,20),dgrid=c(32,64,256)){
#This gets Table 1 of the main paper, the expected number of local maxima.  
#
#The arguments are:
# nugrid, etagrid - vectors of Matern shape and scale parameters;
# dgrid - the sizes of fields to consider. 
out=NULL
#Loop through field sizes:
for(jt in 1:3){
d=dgrid[jt]
#Loop through parameter choices:
for(it in 1:4){
nu=nugrid[it]
eta=etagrid[it]
#Find the expected number of local maxima using each definition of connection:
en1=round(getenmax(d=d,eta=eta,,nu=nu, nbr='cross')[1],1)
en2=round(getenmax(d=d,eta=eta,,nu=nu, nbr='square')[1],1)
out=rbind(out,c(d,nu,eta,en1,en2))
}
}
lbs=c("d","nu","eta","ecross","esquare")
dimnames(out)=list(NULL,lbs)
out
}
###############################################################
gettab2=function(useseed=TRUE,seed=10,nsim=1000){
#This gets Table 2 of the main paper, the variance of the number of local maxima by simulation and by analytic approximation.
#
#The arguments are:
# useseed (TRUE/FALSE) - whether or not to use a seed for repeatability;
# seed - the seed to use, default 10 for the paper;
# nsim - the number of simulations.
if(useseed==TRUE) set.seed(seed)
#Set up.  Use three different dimensions.
tab2=matrix(0,3,4)
dd=c(32,64,256)
#Loop through two types of connection and the three dimensions:
for(it in 1:2){
for(jt in 1:3){
if(it==1) {nbruse='cross'
usq=FALSE} else {nbruse='square'
usq=TRUE}
duse=dd[jt]
#The analytic approximation using own routine:
theory=round(getenmax(d=duse,nbr=nbruse),1)
tab2[jt,2*it]=theory[2]
#The simulation results using own routine:
sims=round(ensims(d=duse,usesquare=usq,nsim=nsim)$nmax,1)
tab2[jt,2*it-1]=sims[2]
}}
list(tab2=tab2)
}

##########################################################
getcordat=function(nsim=50,d=256,lags=c(1:3,5,10,25,50)){
#This is to compare the correlations with target for the different distributions.
#
#The arguments are:
# nsim - the number of simulations;
# d   - for a d*d field;
# lags - a vector of distances to evaluate correlations at.
#

#Use a seed for repeatability.  Seed 100 for results in the paper.
set.seed(100)
#Loop through the five distributions:
dst=c('Gauss','chisq1','chisq3','T','F')
out=NULL
for(jt in 1:5){
dd=dst[jt]
out1=NULL
#Loop through the simulations, estimating and keeping the correlations:
for(it in 1:nsim){
cat("Model ",jt,"  Sim",it,"\n")
#Using own routines:
dat=simdat(d=d,dist=dd)$fmat
cc=getempcor1(dat,lags=lags)$empcor
out1=rbind(out1,cc)
}
#Summarise the correlations and make a nice table:
out=rbind(out,apply(out1,2,mean),apply(out1,2,sd))
}
dimnames(out)=list(rep(dst,each=2),as.character(lags))
round(out,3)
}

#######################################################################
geomsims=function(nsim=50,d=256,usesquare=FALSE,usedefault=TRUE){
#This routine simulates data under each of the five distributions and produces the geometric summary statistics for the associated persistence diagrams.
#
#The arguments are
# nsim - the number of simulations;
# d - the size of the simulated fields, ie d*d;
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE);
# usedefault (TRUE/FALSE) - whether to use default parameter values in simulations.
#
#Loop through the five distributions:
lbs=c('Gauss','chisq1','chisq3','T','F')
out=NULL
for(jt in 1:5){
#Loop through the simulations:
for(it in 1:nsim){
cat("Model ",jt,"  Sim",it,"\n")
#Generate the data:
dat=simdat(d=d,dist=lbs[jt],usesquare=usesquare,usedefault=usedefault)
#Loop through components (uty=0) and holes (uty=1): 
for(uty in 0:1){
#Loop through four different choices of convex peel sizes:
for(prop in c(0.99,0.95,0.9,0.85)){
smat=dat$pers[dat$pers[,1]==uty,2:3]
np=dim(smat)[1]
if(uty==0) nprh=dat$nmin else nprh=dat$nmax
#Calculate the geometric summaries, using own routine, but do not plot:
ch=convpeelplot(smat[,1],smat[,2],draw=FALSE,prop=prop)
gsums=c(uty,prop,jt,mean(smat[,1]),mean(smat[,2]),chsum(ch),np,nprh)
#Combine the results:
out=rbind(out,gsums)
}}
}
}
#The columns of the output are the feature type (betti=0 for components, 1 for holes), the size of the convex peel, the distribution number, the mean birth and death times, the centroids, area, perimeter and filamentarity of the final convex hull, and the counts of features.  
dimnames(out)=list(NULL,c("betti","prop","dist","mn1","mn2","cent1","cent2","area","perim","fil","npTDA","npRH"))
out
}

#################################################################
geomplot=function(d=256,nsim=50,prop=0.9,usesquare=FALSE,useplot=FALSE){
#This routine generates and plots the geometric summary statistics for simulated data.
#
#The arguments are:
# d - the size of the simulated fields, ie d*d;
# nsim - the number of simulations;
# prop - the size of the convex peels;
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE);
# useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.

#We use a seed (100) for repeatability. It is the same as in as in getcordat to ensure comparability.
set.seed(100)
if(useplot==TRUE) postscript("geomfig.ps")
#Get the simulations and summary statistics using own routine:
gsum=geomsims(nsim=nsim,d=d,usesquare=usesquare,usedefault=TRUE)
#The expected values and standard deviations for the default Gaussian random field:
if(usesquare==FALSE){
e=7786.5
s=66}else{
e=3929.5
s=45.2}
#Produce six plots as 2*3 array.  
par(mfrow=c(2,3))
#The first row is components and the second is holes:
for(uty in 0:1){
#Get the summaries and plot:
utys=gsum[,1]
props=gsum[,2]
ii=(utys==uty)&(props==prop)
dist=gsum[ii,3]
mn1=gsum[ii,4]
mn2=gsum[ii,5]
cent1=gsum[ii,6]  
cent2=gsum[ii,7]
area=gsum[ii,8]
perim=gsum[ii,9]
fil=gsum[ii,10]
np=gsum[ii,11]
nprh=gsum[ii,12]
all=gsum[ii,4:11]
i1=dist==1
i2=dist==2
i3=dist==3 
i4=dist==4
i5=dist==5
lbs=c('Gauss','chisq1','chisq3','T','F')
lev=range(c(cent1,cent2))
plot(cent1[i1],cent2[i1],ylim=lev,xlim=lev,ylab="Death level, Cd",xlab="Birth level, Cb",pch=4,cex=0.7)
title("Centroids")
a=min(lev)
b=max(lev)-min(lev)
if(uty==0) legend(a+0.2*b,a+0.3*b,lbs,pch=c(4,6,2,1,3),col=c(1,2,3,4,6),cex=1,bty='n')
points(cent1[i2],cent2[i2],col=2,pch=6,cex=0.7)
points(cent1[i3],cent2[i3],col=3,pch=20,cex=0.7)
points(cent1[i4],cent2[i4],col=4,pch=1,cex=0.7)
points(cent1[i5],cent2[i5],col=6,pch=3,cex=0.7)
plot(area[i1],perim[i1],xlim=range(area),ylim=range(perim),ylab="Perimeter, P",xlab="Area, A",pch=4,cex=0.7)
title("Area and perimeter")
points(area[i2],perim[i2],pch=6,col=2,cex=0.7)
points(area[i3],perim[i3],pch=20,col=3,cex=0.7)
points(area[i4],perim[i4],pch=1,col=4,cex=0.7)
points(area[i5],perim[i5],pch=3,col=6,cex=0.7)
if(uty==0) xl="No. components" else xl="No. holes"
plot(nprh[i1],fil[i1],xlim=range(nprh),ylim=range(fil),ylab="Filamentarity, F",xlab=xl,pch=4,cex=0.7)
title(paste(xl,"and filamentarity"))
points(nprh[i2],fil[i2],pch=6,col=2,cex=0.7)
points(nprh[i3],fil[i3],pch=20,col=3,cex=0.7)
points(nprh[i4],fil[i4],pch=1,col=4,cex=0.7)
points(nprh[i5],fil[i5],pch=3,col=6,cex=0.7)
abline(v=e)
abline(v=(e+2*s),lty=2)
abline(v=(e-2*s),lty=2)
}
par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()
}
######################################################
convpeelplot=function(x,y,prop=0.95,draw=TRUE){
#This routine produces the convex peel of a scatterplot.
#
#The arguments are:
# x,y - the scatterplot coordinates;
# prop - the size of the convexpeel;
# draw (TRUE/FALSE) - whether or not to plot the data and convex peel. 
#
#Plot the persistence diagram if requested:
yl=range(c(x,y))
if(draw==TRUE) {plot(x,y,xlim=yl,ylim=yl,xlab="Level",ylab="Level",pch=" ")
lines(yl,yl)}
#Get the successive convex peels until only proportion prop of points are left:
if(prop==1){
ch=getconvhull(cbind(x,y))
} else{
n=length(x)
nleft=n
dat1=cbind(x,y)
cl=c(2,4)
while(nleft>(prop*n)){
nleft=length(dat1[,1])
ch=getconvhull(dat1)
if(draw==TRUE) lines(ch$xch,ch$ych,col=grey(0.7))
#Discard the last peel:
if(length(ch$ch)<nleft) dat1=dat1[-ch$ch[-1],] else nleft=0
}
}
if(draw==TRUE) points(ch$xch,ch$ych,col=2,pch=20,type="b")
ch
}
##################################################################
chsum=function(ch){
#This routine calculates the summaries of convex hull (a polygon).
#The only argument is the output from convpeelplot.
#
#Get the coordinates of the vertices:
x=ch$xch
y=ch$ych
#The following calculates centroids, area, perimeter and filamentarity of the polygon:
n=length(x)
ii=1:(n-1)
sgnA=0.5*sum(x[ii]*y[ii+1]-x[ii+1]*y[ii])
#Area:
A=abs(sgnA)
#Centroids:
cx=sum( (x[ii]+x[ii+1])*(x[ii]*y[ii+1]-x[ii+1]*y[ii]))/(6*sgnA)
cy=sum( (y[ii]+y[ii+1])*(x[ii]*y[ii+1]-x[ii+1]*y[ii]))/(6*sgnA)
#Perimeter:
P=sum(  sqrt(  (x[ii+1]-x[ii])^2+(y[ii+1]-y[ii])^2))  
#Filamentarity:
F=(P^2-4*pi*A)/(P^2+4*pi*A) 
c(cx,cy,A,P,F)
}
#############################################################################
splitsims=function(simstart=1,simend=100,usesquare=FALSE, prop=0.9,useseed=TRUE){
#This routine allows simulations to be obtained in batches and then split into subregions.
#
#The arguments are:
# simstart, simend - where is this batch in a larger set of results;
# usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE);
# prop - the size of convexpeel to use;
# useseed (TRUE/FALSE) - whether or not to use a seed for repeatability.
#
#The seed is batch-specific:
if(useseed==TRUE) set.seed(simstart+simend)
#The five distributions:
lbs=c('Gauss','chisq1','chisq3','T','F')
#We will generate a 256*256 field and then select the 9 subregions defined by these coordinates (in both x and y directions):
imat=matrix(c(1:64,97:160,193:256),nrow=3,byrow=TRUE)
out=NULL
#Loop through the simulations:
for(it in simstart:simend){
print(it)
#Loop through the distributions:
for(jt in 1:5){
fmatbig=simdat(d=256,dist=lbs[jt],usesquare=usesquare,onlydat=TRUE)$fmat
subst=0
#Loop through the nine subregions:
for(kt in 1:3){
for(lt in 1:3){
subst=subst+1
#Select the subregion:
fmat=fmatbig[imat[kt,],imat[lt,]]
#get the local minima and maxima locations and count:
fmin=getlocmin(fmat,usesquare=usesquare)
nmin=sum(as.vector(fmin))
ftem=-fmat
fmin=getlocmin(ftem,usesquare=usesquare)
nmin=sum(as.vector(fmin))
#Get the persistence diagram for the subset using the TDA routine:
pers=gridDiag(FUNvalues=fmat)$diagram
#Loop through components, holes:
for(uty in 0:1){
#Get the geometric summaries for this persistence diagram:
smat=pers[pers[,1]==uty,2:3]
np=dim(smat)[1]
nprh=nmin
ch=convpeelplot(smat[,1],smat[,2],draw=FALSE,prop=prop)
gsums=c(uty,prop,jt,mean(smat[,1]),mean(smat[,2]),chsum(ch),np,nprh)
#Keep the results:
out=rbind(out,c(it,subst,gsums))
}
}}}}
dimnames(out)=list(NULL,c("sim","subset","betti","prop","dist","mn1","mn2","cent1","cent2","area","perim","fil","npTDA","npRH"))
data.frame(out)
}
#################################################################################
splittest=function(dist1=1,dist2=1,prop=0.9,adjsig=0.025, splits=read.table("simsums.dat",header=TRUE)){
#This routine examimes the performance of tests for differences between two fields, using simulated data.
#
#The arguments are:
#  dist1, dist2 - the distributions (1-5 for Gauss, chisq1, chisq3, T, F) being compared;
# prop - the choice of convex peel size to use;
# adjsig - the significance level to use for the test (0.025 for two tests);
# splits - simulation summaries, as produced by splitsims or (by default) read from file.
#
#Set up:
nsim=dim(splits)[1]/90
nrep=nsim/2
#Pick out the appropriate data from the larger results matrix, dividing into first and second halves of the simulations so no overlap.  First components, then holes:
splits10=splits[(splits$sim<=nrep)&(splits$dist==dist1)&(splits$prop==prop)&(splits$betti==0),]
splits20=splits[(splits$sim>nrep)&(splits$dist==dist2)&(splits$prop==prop)&(splits$betti==0),]
splits11=splits[(splits$sim<=nrep)&(splits$dist==dist1)&(splits$prop==prop)&(splits$betti==1),]
splits21=splits[(splits$sim>nrep)&(splits$dist==dist2)&(splits$prop==prop)&(splits$betti==1),]
pall=NULL
#Loop through half the simulations:
for(it in 1:nrep){
#Pick a pair of data sets to compare, again both components and holes: 
dat10=splits10[splits10$sim==it,]
dat20=splits20[splits20$sim==(nrep+it),]
dat11=splits11[splits11$sim==it,]
dat21=splits21[splits21$sim==(nrep+it),]
#Concatenat teady for testing
dat0=rbind(dat10,dat20)
dat1=rbind(dat11,dat21)
sim=dat0$sim
fil=dat1$fil
np=dat0$npRH
#This is a labelling:
gp=sim*0
gp[sim==(unique(sim)[2])]=1
#Wilcoxon tests for filamentarity for holes and number of features for components:
p1=wilcox.test(fil~gp)$p.value
p2=wilcox.test(np~gp)$p.value
pval=min(c(p1,p2))
pall=rbind(pall,c(p1,p2,pval))
}
#Power estimates:
lev=c(0.05,0.05,adjsig)
ns=rep(0,length=3)
for(it in 1:3){
ns[it]=sum(pall[,it]<lev[it])
}
pow=ns/nrep
pow
}
##############################################
splittests=function(simsums=read.table("simsums.dat",header=TRUE)){
#Routine to estimate power of the non-parametric method for comparing two fields.
#
#The argument is:
# simsums - simulation summaries, as produced by splitsims or (by default) read from file.
#The warnings about tied data in the Wilcoxon test can be ignored.
#
#Look at all pairs of distributions and use own routine to estimate size and power:  
pow=NULL
for(dist1 in 1:5){
for(dist2 in dist1:5){
pp=splittest(dist1=dist1,dist2=dist2,prop=0.9,splits=simsums)
print(c(dist1,dist2,pp))
pow=rbind(pow,c(dist1,dist2,pp))
}}
pow
}
##################################################################
getappfig=function(prop=0.9,useplot=FALSE){
#Get Figure 5 of the main paper, showing convex peels and cumulative counts for the three GASS regions.
#
#The arguments are;
# prop - the size of convex peel to use.  The paper used prop=0.9;
## useplot (TRUE/FALSE) - an indicator for whether the plot should be saved to file or not.
#
#Set up the plot and plot limits
xl=c(-4,4)
yl=c(-4,4)
lw=2
if(useplot==TRUE) postscript("appfig.ps",horizontal=FALSE)
par(mfrow=c(2,2))
#Find and plot the convex peels for components using own routines, and pre-stored persistence diagrams:
ch1=convpeelplot1(obs1$pers,usetype=1,prop=prop,draw=FALSE)
ch2=convpeelplot1(obs2$pers,usetype=1,prop=prop,draw=FALSE)
ch3=convpeelplot1(obs3$pers,usetype=1,prop=prop,draw=FALSE)
plot(ch1$xch,ch1$ych,type="l",xlab="Birth level",ylab="Death level",xlim=xl,ylim=yl,col=2,lwd=lw,main="(a) Components, 90% peel",cex.main=1)
lines(xl,yl)
lines(ch2$xch,ch2$ych,lty=2,col=3,lwd=lw)
lines(ch3$xch,ch3$ych,lty=3,col=4,lwd=lw)
legend(0,-1,c("Region 1","Region 2","Region 3"),lty=1:3,col=2:4,lwd=lw,bty='n')
#Repeat for holes:
ch1=convpeelplot1(obs1$pers,usetype=2,prop=prop,draw=FALSE)
ch2=convpeelplot1(obs2$pers,usetype=2,prop=prop,draw=FALSE)
ch3=convpeelplot1(obs3$pers,usetype=2,prop=prop,draw=FALSE)
plot(ch1$xch,ch1$ych,type="l",xlab="Birth level",ylab="Death level",xlim=xl,ylim=yl,col=2,lwd=lw, main="(b) Holes, 90% peel",cex.main=1)
lines(xl,yl)
lines(ch2$xch,ch2$ych,lty=2,col=3,lwd=lw)
lines(ch3$xch,ch3$ych,lty=3,col=4,lwd=lw)
#For each region, get a counting process of cumulative numbers of features, using own routine:
cp1=getcp(dat=obs1)
cp2=getcp(dat=obs2)
cp3=getcp(dat=obs3)
yl=c(0,6000)
#Plot these processes and, for Region 3, +/- 2 approx standard deviations:
plot(cp1$zmin,cp1$N1,type="s",xlab="Level",ylab="Cumulative components",xlim=xl,ylim=yl,lwd=lw,col=2,main="(c) Component count",cex.main=1 )
lines(cp2$zmin,cp2$N1,type="s",lwd=lw,col=3,lty=2)
lines(cp3$zmin,cp3$N1,type="s",lwd=lw,col=4,lty=3)
lines(cp3$zmin,cp3$U1,col=grey(0.5),lty=4,lwd=1)
lines(cp3$zmin,cp3$L1,col=grey(0.5),lty=4,lwd=1)
abline(v=0,col=grey(0.5))
plot(cp1$zmax,cp1$N2,type="s",xlab="Level",ylab="Cumulative holes",xlim=xl,ylim=yl,lwd=lw,col=2,main="(d) Hole count",cex.main=1)
lines(cp2$zmax,cp2$N2,type="s",lwd=lw,col=3,lty=2)
lines(cp3$zmax,cp3$N2,type="s",lwd=lw,col=4,lty=3)
lines(cp3$zmax,cp3$U2,col=grey(0.5),lty=4,lwd=1)
lines(cp3$zmax,cp3$L2,col=grey(0.5),lty=4,lwd=1)
abline(v=0,col=grey(0.5))
par(mfrow=c(1,1))
if(useplot==TRUE) dev.off()

}
##########################################################
getcp=function(dat=obs1){
#This routine gets the cumulative numbers of components and holes as the filtration proceeds, together with an approximate standard deviation.
#
#The argument is one of the lists produced by setrealdat.
#
#The field values and the locations of local minima and maxima:
z1=as.vector(dat$fmat)
imin=as.vector(dat$fmin)
imax=as.vector(dat$fmax)
zmin=sort(z1[imin])
zmax=sort(z1[imax])
zmin=c(min(zmin)-0.1,zmin,max(zmin)+0.1)
zmax=c(min(zmax)-0.1,zmax,max(zmax)+0.1)
nmin=length(zmin)
nmax=length(zmax)
#The counting processes and the estimated variances:
N1=0:(nmin-1)
N2=0:(nmax-1)
pv1=getestvar(zmin,N1)
pv2=getestvar(zmax,N2)
#Approximate 95% confidence intervals:
U1=N1+2*sqrt(pv1)
L1=N1-2*sqrt(pv1)
U2=N2+2*sqrt(pv2)
L2=N2-2*sqrt(pv2)
list(N1=N1,N2=N2,L1=L1,U1=U1,L2=L2,U2=U2,zmin=zmin,zmax=zmax,pv1=pv1,pv2=pv2)
}
######################################################
getestvar=function(z,N1,p=15,ng=50000){
#This routine estmates the variance of the counting process, based on a smoothed cumulative intensity with Poisson approximation.
#
#The arguments are:
# z - a vector of levels at which the counting process has increments;
# N1 - the counting process;
# p - the order of polynomial for smoothing;
# ng - the grid size for evaluating the smooth version.
#
#Form the polynonial terms:
x=rep(1,length=length(z))
for(it in 1:p){
x=cbind(x,(z^it))
}
#Smooth via a polynomial regression:
fit=lm(N1~-1+x)
y=fit$fitted.values
b=fit$coefficients
#Make a grid spanning the range of z:
zg=seq(min(z)-0.1,max(z)+0.1,length=ng)
#Create the smoothed values:
yd=b[2]*rep(1,length=length(ng))
for(it in 2:p){
yd=yd+it*b[it+1]*zg^(it-1)
}
#Increment and binary variance of increment in N1:
dg=zg[2]-zg[1]
ov=yd*dg*(1-yd*dg)
ov[ov<0]=0
#Cumulative variance:
ov=cumsum(ov)
#Another polynomial, over grid values
x=rep(1,length=length(zg))
for(it in 1:p){
x=cbind(x,zg^it)
}
#Fit to the cumulative variance:
fit=lm(ov~-1+x)
b=fit$coefficients
#Back to the original z instead of the grid:
x=rep(1,length=length(z))
for(it in 1:p){
x=cbind(x,z^it)
}
#Find the variance estimates:
pvsm=x%*%b
pvsm[pvsm<0]=0
pvsm
}
#####################################################
getcountsum=function(){
#This routine produces Table 5 of the main paper, a comparison of observed and expected counts.
#
#Get the expected counts based on empirical correlations:
out1=round(getenmax(empirical=TRUE,anis=FALSE,fmat=obs1$fmat),1)
out2=round(getenmax(empirical=TRUE,anis=FALSE,fmat=obs2$fmat),1)
out3=round(getenmax(empirical=TRUE,anis=FALSE,fmat=obs3$fmat),1)
#The observed counts were already obained during setrealdat()
out1=c(obs1$nmin,obs1$nmax,out1)
out2=c(obs2$nmin,obs2$nmax,out2)
out3=c(obs3$nmin,obs3$nmax,out3)
#Tabulate:
out=rbind(out1,out2,out3)
lbs=c("Components", "Holes", "Gauss-expected","Gauss-SD")
dimnames(out)=list(c("Region 1","Region 2","Region 3"),lbs)
out
}
################################################################
compregions=function(){
#This routine produces Table 6 of the main paper, a comparison between the three regions.
#
#The warning messages about ties in the data can be ignored.
#
#Bind together the geometric summaries of the nine sub-regions for each main region:
gsum=rbind(splitone(obs1,region=1),splitone(obs2,region=2),splitone(obs3,region=3))
#Test for differences:
qmat=splitex(gsum)
qmat
}
#################################################################
splitex=function(gsum){
#This routine calculates 
#
#The argument is:
# gsum - geometric summaries, produced by compregions.  Region is 
#Initial set up
gsum=gsum[,-1]
q12=q13=q23=matrix(0,2,2)
#Loop through components (uty=0) and holes (uty=1):
for(uty in 0:1){
#Extract the summary statistics:
utys=gsum[,1]
props=gsum[,2]
ii=(utys==uty)
region=gsum[ii,3]
mn1=gsum[ii,4]
mn2=gsum[ii,5]
cent1=gsum[ii,6]  
cent2=gsum[ii,7]
area=gsum[ii,8]
perim=gsum[ii,9]
fil=gsum[ii,10]
np=gsum[ii,11]
nprh=gsum[ii,12]
i12=(region!=3)
i13=(region!=2)
i23=(region!=1)
#Test region 1 v region2, region1 v region 3, and region 2 v region 3.
#Use filamentarity and number of features. Ignore any warning messages about ties:
p11=wilcox.test(fil[i12]~region[i12])$p.value
p21=wilcox.test(np[i12]~region[i12])$p.value
if(uty==0) q12[1,]=c(p11,p21) else q12[2,]=c(p11,p21)
p12=wilcox.test(fil[i13]~region[i13])$p.value
p22=wilcox.test(np[i13]~region[i13])$p.value
if(uty==0) q13[1,]=c(p12,p22) else q13[2,]=c(p12,p22)
p13=wilcox.test(fil[i23]~region[i23])$p.value
p23=wilcox.test(np[i23]~region[i23])$p.value
if(uty==0) q23[1,]=c(p13,p23) else q23[2,]=c(p13,p23)
all=gsum[ii,4:11]
}
#Tabulate p-values:
qmat=round(rbind(q12,q13,q23),3)
lbs=c("Components","Holes")
lbs=c(lbs,lbs,lbs)
dimnames(qmat)=list(lbs,c("Filamentarity","Number"))
qmat
}


#############################################################################
splitone=function(dat,region=1,usesquare=FALSE){
#This routine calculates the geometric summary statistics for the nine sub-regions of a field.
#
#The arguments are:
# dat - the field data under consideration:
# region - region
#  usesquare (TRUE/FALSE) -  to decide whether neighbouring values are connected if they share only an edge (FALSE) or either an edge or vertex (TRUE).
#
#This defines the sub-regions of a 256*256 field.
imat=matrix(c(1:64,97:160,193:256),nrow=3,byrow=TRUE)
out=NULL
fmatbig=dat$fmat
subst=0
#Loop through the sub-regions:
for(kt in 1:3){
for(lt in 1:3){
subst=subst+1
#Find the minima, maxima data and the persistence diagrams for each sub-region:
fmat=fmatbig[imat[kt,],imat[lt,]]
fmin=getlocmin(fmat,usesquare=usesquare)
nmin=sum(as.vector(fmin))
ftem=-fmat
fmin=getlocmin(ftem,usesquare=usesquare)
nmax=sum(as.vector(fmin))
pers=gridDiag(FUNvalues=fmat)$diagram
#Loop through components and holes:
for(uty in 0:1){
#Prop is the size of the convex peel.
#Only using one value now (previous versions had more):
for(prop in 0.9){
smat=pers[pers[,1]==uty,2:3]
np=dim(smat)[1]
if(uty==0) nprh=nmin else nprh=nmax
#Get the convex peel and the geometric summaries: 
ch=convpeelplot(smat[,1],smat[,2],draw=FALSE,prop=prop)
gsums=c(uty,prop,region,mean(smat[,1]),mean(smat[,2]),chsum(ch),np,nprh)
out=rbind(out,c(subst,gsums))
}}
}}
#Write the output as a list:
dimnames(out)=list(NULL,c("subset","betti","prop","region","mn1","mn2","cent1","cent2","area","perim","fil","npTDA","npRH"))
out
}
##################################################
getdistsum=function(){
#This routine finds Table 7 of the main paper, the bottleneck and Wasserstein distances between sub-regions.
#
#Get the persistence diagrams for the sub-regions:
obssplitpers=getrealpers()
#Get the bottleneck distances:
sumsb=getbots(obssplitpers)
#get the first and second Wasserstein distances:
sumsw1=getwass(obssplitpers,p=1)
sumsw2=getwass(obssplitpers,p=2)
#Collate the results:
all1=cbind(sumsb[,1:3],sumsw1[,2:3],sumsw2[,2:3])
all2=cbind(sumsb[,c(1,4,5)],sumsw1[,4:5],sumsw2[,4:5])
all=rbind(all1,all2)
round(all,2)
}
###########################################################################
getbots=function(dat){
#This routine gets the bottleneck distances between persistence diagrams of sub-regions.
#
#The argument is
# dat - the output of getrealpers
#
num=dat[,3]
uty=dat[,4]
smatbig=dat[,4:6]
out=NULL
#Loop through pairs from the nine subroutines in each of three regions:
for(it in 1:26){
ii=(num==it)
for(jt in (it+1):27){
print(c(it,jt))
jj=(num==jt)
#The persistence diagrams:
smat1=smatbig[ii,]
smat2=smatbig[jj,]
#Get the bottleneck distance using a routine from the TDA package, for each of components and holes:
b1=bottleneck(smat1,smat2,0)
b2=bottleneck(smat1,smat2,1)
out=rbind(out,c(it,jt,b1,b2))
}
}
#Format the output:
o1=o2=out[,1]*0+3
o1[out[,1]<=18]=2
o1[out[,1]<=9]=1
o2[out[,2]<=18]=2
o2[out[,2]<=9]=1
out=cbind(o1,o2,out)
sums=NULL
for(it in 1:3){
for(jt in it:3){
i1=out[,1]==it
i2=out[,2]==jt
b0=out[i1&i2,5]
b1=out[i1&i2,6]
sums=rbind(sums,c(it,jt,length(b0),mean(b0),sd(b0),length(b1),mean(b1),sd(b1)))
}}
print(sums)
lbs1=c("Region 1 Region 1","Region 1 Region 2","Region 1 Region 3","Region 2 Region 2","Region 2 Region 3","Region 3 Region 3")
lbs2=c("M","Mean","SD","Mean","SD")
sums=sums[,-c(1,2,6)]
sums=round(sums,2)
dimnames(sums)=list(lbs1,lbs2)
sums
}

###########################################################################
getwass=function(dat,p=1){
#This routine gets the pth bottleneck distances between persistence diagrams of sub-regions.
#
#The arguments are:
# dat - the output of getrealpers
# p - the power in the Wasserstein expression.
num=dat[,3]
uty=dat[,4]
smatbig=dat[,4:6]
out=NULL
#Loop through pairs from the nine subroutines in each of three regions:
for(it in 1:26){
ii=(num==it)
for(jt in (it+1):27){
print(c(it,jt))
jj=(num==jt)
#The persistence diagrams:
smat1=smatbig[ii,]
smat2=smatbig[jj,]
#Get the Wasserstein distance using a routine from the TDA package, for each of components and holes:
b1=wasserstein(smat1,smat2,p=p, dimension=0)
b2=wasserstein(smat1,smat2,p=p,dimension=1)
out=rbind(out,c(it,jt,b1,b2))
}
}
#Format the output:
o1=o2=out[,1]*0+3
o1[out[,1]<=18]=2
o1[out[,1]<=9]=1
o2[out[,2]<=18]=2
o2[out[,2]<=9]=1
out=cbind(o1,o2,out)
sums=NULL
for(it in 1:3){
for(jt in it:3){
i1=out[,1]==it
i2=out[,2]==jt
b0=out[i1&i2,5]
b1=out[i1&i2,6]
sums=rbind(sums,c(it,jt,length(b0),mean(b0),sd(b0),length(b1),mean(b1),sd(b1)))
}}
print(sums)
lbs1=c("Region 1 Region 1","Region 1 Region 2","Region 1 Region 3","Region 2 Region 2","Region 2 Region 3","Region 3 Region 3")
lbs2=c("M","Mean","SD","Mean","SD")
sums=sums[,-c(1,2,6)]
sums=round(sums,2)
dimnames(sums)=list(lbs1,lbs2)
sums
}
############################################################################
getrealpers=function(){
#This routine finds the persistence diagrams for the sub-regions of the three GASS regions.
#
#This defines the sub-regions of a 256*256 matrix:
imat=matrix(c(1:64,97:160,193:256),nrow=3,byrow=TRUE)
out=NULL
num=0
#Loop through the regions:
for(use in 1:3){
if(use==1) dat=obs1
if(use==2) dat=obs2
if(use==3) dat=obs3
fmatbig=dat$fmat
subst=0
#Loop through the sub-regions:
for(kt in 1:3){
for(lt in 1:3){
subst=subst+1
num=num+1
#Get the persistence diagram:
fmat=fmatbig[imat[kt,],imat[lt,]]
pers=gridDiag(FUNvalues=fmat)$diagram
tem=cbind(use,subst,num,pers)
out=rbind(out,tem)
}
}
}
out
}





