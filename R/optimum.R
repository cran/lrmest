optimum <-
function(formula,r,R,dpn,delt,aa1,aa2,aa3,k,d,press=FALSE,data=NULL,na.action,...)
{
k<-as.matrix(k)
d<-as.matrix(d)
j<-0
i<-0
for (j in 1:NROW(k))
{
for (i in 1:NROW(d))
{
fltn1<-file("ltn1.data","a+")
cat(k[j],d[i],"\n",file=fltn1,append=TRUE)
close(fltn1)
}
}
mat1<-read.table("ltn1.data")
unlink("ltn1.data")
mat1<-as.matrix(mat1)
varoglt <- function(formula, k, d, data = NULL,na.action, ...) 
{
  voglt1 <- oglt1(formula, k, d, data, na.action, ...)
  voglt2 <- oglt2(formula, k, d, data, na.action, ...)
  voglt3 <- oglt3(formula, k, d, data, na.action, ...)
  voglt1m <- matrix(voglt1, NROW(k) * NROW(d))
  mseoglt1 <- cbind(voglt1m, mat1)
  smoglt1 <- mseoglt1[order(mseoglt1[, 1L]), ]
  minoglt1 <- smoglt1[1L, ]
  voglt2m <- matrix(voglt2, NROW(k) * NROW(d))
  mseoglt2 <- cbind(voglt2m, mat1)
  smoglt2 <- mseoglt2[order(mseoglt2[, 1L]), ]
  minoglt2 <- smoglt2[1L, ]
  voglt3m <- matrix(voglt3, NROW(k) * NROW(d))
  mseoglt3 <- cbind(voglt3m, mat1)
  smoglt3 <- mseoglt3[order(mseoglt1[, 1L]), ]
  minoglt3 <- smoglt3[1L, ]
  voglt1cm <- as.matrix(minoglt1)
  voglt2cm <- as.matrix(minoglt2)
  voglt3cm <- as.matrix(minoglt3)
  aval <- rbind(minoglt1, minoglt2, minoglt3)
  colnames(aval) <- c("Optimum_mse", "Optimum_k", "Optimum_d")
  rownames(aval) <- c("OGLTE1", "OGLTE2", "OGLTE3")
  aval
}
OGLTE <- varoglt(formula, k, d, data, na.action)
varogalt<-function(formula,aa1,aa2,aa3,k,d,data=NULL,na.action,...)
{
  vogalt1<-ogalt1(formula,k,d,aa1,data,na.action,...)
  vogalt2<-ogalt2(formula,k,d,aa2,data,na.action,...)
  vogalt3<-ogalt3(formula,k,d,aa3,data,na.action,...)
  vogalt1m<-matrix(vogalt1,NROW(k)*NROW(d))
  mogalt1<-cbind(vogalt1m,mat1)
  smogalt1<-mogalt1[order(mogalt1[,1L]),]
  minogalte1<-smogalt1[1L,]
  vogalt2m<-matrix(vogalt2,NROW(k)*NROW(d))
  mogalt2<-cbind(vogalt2m,mat1)
  smogalt2<-mogalt2[order(mogalt2[,1L]),]
  minogalte2<-smogalt2[1L,]
  vogalt3m<-matrix(vogalt3,NROW(k)*NROW(d))
  mogalt3<-cbind(vogalt3m,mat1)
  smogalt3<-mogalt3[order(mogalt3[,1L]),]
  minogalte3<-smogalt3[1L,]         
  vogalt1mat<-as.matrix(minogalte1)
  vogalt2mat<-as.matrix(minogalte2)
  vogalt3mat<-as.matrix(minogalte3)
  aval<-rbind(minogalte1,minogalte2,minogalte3)
  colnames(aval)<-c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(aval)<-c("OGALTE1","OGALTE2","OGALTE3")
  aval
}
OGALTE<-varogalt(formula,aa1,aa2,aa3,k,d,data,na.action)
varlt<-function(formula,k,d,press=FALSE,data=NULL,na.action,...)
{
vlt1<-lte1(formula,k,d,press,data,na.action,...)
vlt2<-lte2(formula,k,d,press,data,na.action,...)
vlt3<-lte3(formula,k,d,press,data,na.action,...)
vlt1m<-matrix(vlt1,NROW(k)*NROW(d))
mselt1<-cbind(vlt1m,mat1)
smlt1<-mselt1[order(mselt1[,1L]),]
minlt1<-smlt1[1L,]
vlt2m<-matrix(vlt2,NROW(k)*NROW(d))
mselt2<-cbind(vlt2m,mat1)
smlt2<-mselt2[order(mselt2[,1L]),]
minlt2<-smlt2[1L,]
vlt3m<-matrix(vlt3,NROW(k)*NROW(d))
mselt3<-cbind(vlt3m,mat1)   
smlt3<-mselt3[order(mselt1[,1L]),] 
minlt3<-smlt3[1L,]          
vlt1cm<-as.matrix(minlt1)
vlt2cm<-as.matrix(minlt2)
vlt3cm<-as.matrix(minlt3)
aval<-rbind(minlt1,minlt2,minlt3)
colnames(aval)<-c("Optimum_mse","Optimum_k","Optimum_d")
rownames(aval)<-c("LTE1","LTE2","LTE3")
aval
}
LTE<-varlt(formula,k,d,press,data,na.action)
varalt<-function(formula,aa1,aa2,aa3,k,d,press=FALSE,data=NULL,na.action,...)
{
valt1<-alte1(formula,k,d,aa1,press,data,na.action,...)
valt2<-alte2(formula,k,d,aa2,press,data,na.action,...)
valt3<-alte3(formula,k,d,aa3,press,data,na.action,...)
valt1m<-matrix(valt1,NROW(k)*NROW(d))
malt1<-cbind(valt1m,mat1)
smalt1<-malt1[order(malt1[,1L]),]
minalte1<-smalt1[1L,]
valt2m<-matrix(valt2,NROW(k)*NROW(d))
malt2<-cbind(valt2m,mat1)
smalt2<-malt2[order(malt2[,1L]),]
minalte2<-smalt2[1L,]
valt3m<-matrix(valt3,NROW(k)*NROW(d))
malt3<-cbind(valt3m,mat1)
smalt3<-malt3[order(malt3[,1L]),]
minalte3<-smalt3[1L,]         
valt1mat<-as.matrix(minalte1)
valt2mat<-as.matrix(minalte2)
valt3mat<-as.matrix(minalte3)
aval<-rbind(minalte1,minalte2,minalte3)
colnames(aval)<-c("Optimum_mse","Optimum_k","Optimum_d")
rownames(aval)<-c("ALTE1","ALTE2","ALTE3")
aval
}
ALTE<-varalt(formula,aa1,aa2,aa3,k,d,press,data,na.action)
varrid<-function(formula,k,data=NULL,na.action,...)
{
vrid<-rid(formula,k,data,na.action,...)
msval<-vrid[1L,]
msval<-as.matrix(msval)
Optimum_k<-msval[1L,]
MSE<-msval[2L,]
Optimum_d<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("RE")
minmat
}
RE<-varrid(formula,k,data,na.action)
varogrid<-function(formula,k,data=NULL,na.action,...)
{
  vogrid<-ogre(formula,k,data,na.action,...)
  msval<-vogrid[1L,]
  msval<-as.matrix(msval)
  Optimum_k<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_d<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGRE")
  minmat
}
OGRE<-varogrid(formula,k,data,na.action)
varliu<-function(formula,d,data=NULL,na.action,...)
{
vliu<-liu(formula,d,data,na.action,...)
msval<-vliu[1L,]
msval<-as.matrix(msval)
Optimum_d<-msval[1L,]
MSE<-msval[2L,]
Optimum_k<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("LE")
minmat
}
LE<-varliu(formula,d,data,na.action)
varogliu<-function(formula,d,data=NULL,na.action,...)
{
  vliu<-ogliu(formula,d,data,na.action,...)
  msval<-vliu[1L,]
  msval<-as.matrix(msval)
  Optimum_d<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_k<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGLE")
  minmat
}
OGLE<-varogliu(formula,d,data,na.action)
varaure<-function(formula,k,data=NULL,na.action,...)
{
vaure<-aur(formula,k,data,na.action,...)
msval<-vaure[1L,]
msval<-as.matrix(msval)
Optimum_k<-msval[1L,]
MSE<-msval[2L,]
Optimum_d<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("AURE")
minmat
}
AURE<-varaure(formula,k,data,na.action)
varogaure<-function(formula,k,data=NULL,na.action,...)
{
  vaure<-ogaur(formula,k,data,na.action,...)
  msval<-vaure[1L,]
  msval<-as.matrix(msval)
  Optimum_k<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_d<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGAURE")
  minmat
}
OGAURE<-varogaure(formula,k,data,na.action)
varaule<-function(formula,d,data=NULL,na.action,...)
{
vaule<-aul(formula,d,data,na.action,...)
msval<-vaule[1L,]
msval<-as.matrix(msval)
Optimum_d<-msval[1L,]
MSE<-msval[2L,]
Optimum_k<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("AULE")
minmat
}
AULE<-varaule(formula,d,data,na.action)
varogaule<-function(formula,d,data=NULL,na.action,...)
{
  vaule<-ogaul(formula,d,data,na.action,...)
  msval<-vaule[1L,]
  msval<-as.matrix(msval)
  Optimum_d<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_k<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGAULE")
  minmat
}
OGAULE<-varogaule(formula,d,data,na.action)
varrl<-function(formula,r,R,delt,d,data=NULL,na.action,...)
{
vrliu<-rliu(formula,r,R,delt,d,data,na.action,...)
msval<-vrliu[1L,]
msval<-as.matrix(msval)
Optimum_d<-msval[1L,]
MSE<-msval[2L,]
Optimum_k<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("RLE")
minmat
}
RLE<-varrl(formula,r,R,delt,d,data,na.action)
varogrl<-function(formula,r,R,delt,d,data=NULL,na.action,...)
{
  vrliu<-ogrliu(formula,r,R,delt,d,data,na.action,...)
  msval<-vrliu[1L,]
  msval<-as.matrix(msval)
  Optimum_d<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_k<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGRLE")
  minmat
}
OGRLE<-varogrl(formula,r,R,delt,d,data,na.action)
varsrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
vsrre<-srre(formula,r,R,dpn,delt,k,data,na.action,...)
msval<-vsrre[1L,]
msval<-as.matrix(msval)
Optimum_k<-msval[1L,]
MSE<-msval[2L,]
Optimum_d<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("SRRE")
minmat
}
SRRE<-varsrr(formula,r,R,dpn,delt,k,data,na.action)
varogsrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
  vsrre<-ogsrre(formula,r,R,dpn,delt,k,data,na.action,...)
  msval<-vsrre[1L,]
  msval<-as.matrix(msval)
  Optimum_k<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_d<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGSRRE")
  minmat
}
OGSRRE<-varogsrr(formula,r,R,dpn,delt,k,data,na.action)
varsrliu<-function(formula,r,R,dpn,delt,d,data=NULL,na.action,...)
{
vsrliu<-srliu(formula,r,R,dpn,delt,d,data,na.action,...)
msval<-vsrliu[1L,]
msval<-as.matrix(msval)
Optimum_d<-msval[1L,]
MSE<-msval[2L,]
Optimum_k<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("SRLE")
minmat
}
SRLE<-varsrliu(formula,r,R,dpn,delt,d,data,na.action)
varogsrliu<-function(formula,r,R,dpn,delt,d,data=NULL,na.action,...)
{
  vsrliu<-ogsrliu(formula,r,R,dpn,delt,d,data,na.action,...)
  msval<-vsrliu[1L,]
  msval<-as.matrix(msval)
  Optimum_d<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_k<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGSRLE")
  minmat
}
OGSRLE<-varogsrliu(formula,r,R,dpn,delt,d,data,na.action)
varols<-function(formula,data=NULL,na.action,...)
{
cal<-match.call(expand.dots=FALSE)
mat<-match(c("formula","data","na.action"), names(cal))
cal<-cal[c(1L,mat)]
cal[[1L]]<-as.name("model.frame") 
cal<-eval(cal)
y<-model.response(cal)         
md<- attr(cal, "terms")
x<-model.matrix(md,cal,contrasts)
s<-t(x)%*%x
xin<-solve(s)
bb<-xin%*%t(x)%*%y
colnames(bb) <- c("Estimate")
ev<-(t(y)%*%y-t(bb)%*%t(x)%*%y)/(NROW(x)-NCOL(x))
ev<-diag(ev)
dbb<-ev*xin
mse<-sum(diag(dbb))
pr<-0
i<-1
m<-1
for (i in 1:nrow(x))
 {
               subsum<-0
bb<-c(x[i,]%*%t(x[i,]))
z<-solve(t(x)%*%x-bb)%*%(t(x)%*%y-x[i,]*y[i])
for (m in 1:ncol(x))      
subsum<-subsum+(x[i,m]%*%z[m])         
pr<-pr+(y[i]-subsum)^2
}
opre<-t(pr)
minpress<-opre
optimal_k<-NA
optimal_d<-NA
premat<-cbind(minpress,optimal_k,optimal_d)
colnames(premat)<-c("Optimun_pres","Optimum_k","Optimum_d")
rownames(premat)<-c("OLS")
optimum_k<-NA
optimum_d<-NA
msval<-cbind(mse,optimum_k,optimum_d)
rownames(msval)<-c("OLSE")
colnames(msval)<-c("Optimum_mse","Optimum_k","Optimum_d")
val<-list(msval,premat)
return(val)
}
mpols<-varols(formula,data,na.action)
msols<-mpols[[1L]]
prols<-mpols[[2L]]
varogols<-function(formula,data=NULL,na.action,...)
{
vogols<-ogols(formula,data,na.action,...)
mogols<-vogols[[2L]]
mogols<-as.matrix(mogols)
optk<-NA
optd<-NA
matval<-cbind(mogols,optk,optd)
colnames(matval)<-c("Optimun_mse","Optimum_k","Optimum_d")
rownames(matval)<-c("OGOLSE")
matval
}
OGOLSE<-varogols(formula,data,na.action)
varmix<-function(formula,r,R,dpn,delt,data=NULL,na.action,...)
{
vmixe<-mixe(formula,r,R,dpn,delt,data,na.action,...)
msval<-vmixe[[2L]]
msval<-as.matrix(msval)
Optimum_k<-NA
Optimum_d<-NA
minmat<-cbind(msval,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("MIXE")
minmat
}
MIXE<-varmix(formula,r,R,dpn,delt,data,na.action)
varogmix<-function(formula,r,R,dpn,delt,data=NULL,na.action,...)
{
  vmixe<-ogmix(formula,r,R,dpn,delt,data,na.action,...)
  msval<-vmixe[[2L]]
  msval<-as.matrix(msval)
  Optimum_k<-NA
  Optimum_d<-NA
  minmat<-cbind(msval,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGMIXE")
  minmat
}
OGMIXE<-varogmix(formula,r,R,dpn,delt,data,na.action)
varrls<-function(formula,r,R,delt,data=NULL,na.action,...)
{
vrls<-rls(formula,r,R,delt,data,na.action,...)
msval<-vrls[[2L]]
msval<-as.matrix(msval)
Optimum_k<-NA
Optimum_d<-NA
minmat<-cbind(msval,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("RLSE")
minmat
}
RLSE<-varrls(formula,r,R,delt,data,na.action)
varogrls<-function(formula,r,R,delt,data=NULL,na.action,...)
{
  vrls<-ogrls(formula,r,R,delt,data,na.action,...)
  msval<-vrls[[2L]]
  msval<-as.matrix(msval)
  Optimum_k<-NA
  Optimum_d<-NA
  minmat<-cbind(msval,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGRLSE")
  minmat
}
OGRLSE<-varogrls(formula,r,R,delt,data,na.action)
varrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
vrr<-rrre(formula,r,R,dpn,delt,k,data,na.action,...)
msval<-vrr[1L,]
msval<-as.matrix(msval)
Optimum_k<-msval[1L,]
MSE<-msval[2L,]
Optimum_d<-NA
minmat<-cbind(MSE,Optimum_k,Optimum_d)
colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
rownames(minmat)=c("RRRE")
minmat
}
RRRE<-varrr(formula,r,R,dpn,delt,k,data,na.action)
varogrr<-function(formula,r,R,dpn,delt,k,data=NULL,na.action,...)
{
  vrr<-ogrrre(formula,r,R,dpn,delt,k,data,na.action,...)
  msval<-vrr[1L,]
  msval<-as.matrix(msval)
  Optimum_k<-msval[1L,]
  MSE<-msval[2L,]
  Optimum_d<-NA
  minmat<-cbind(MSE,Optimum_k,Optimum_d)
  colnames(minmat)=c("Optimum_mse","Optimum_k","Optimum_d")
  rownames(minmat)=c("OGRRRE")
  minmat
}
OGRRRE<-varogrr(formula,r,R,dpn,delt,k,data,na.action)
mat1<-rbind(SRRE,OGSRRE,SRLE,MIXE,OGMIXE,msols,OGOLSE)
mat2<-rbind(RLE,OGRLE,RRRE,OGRRRE,RLSE,OGRLSE,msols)
mat3<-rbind(RE,OGRE,AURE,OGAURE,msols)
mat4<-rbind(LE,OGLE,AULE,OGAULE,LTE,OGLTE,ALTE,OGALTE,msols)
mat1<-round(mat1,digits<-4L)
mat2<-round(mat2,digits<-4L)
mat3<-round(mat3,digits<-4L)
mat4<-round(mat4,digits<-4L)
val<-list("Comparison of Stochastic Restricted Type Estimators"=mat1,"Comparison of Restricted Type Estimators"=mat2,"Comparison of Ridge Type Estimators"=mat3,"Comparison of Liu Type Estimators"=mat4)  
preval<-rbind(LTE,ALTE,prols)
preval<-round(preval,digits<-4L)
colnames(preval)<-c("Optimum_press","Optimum_k","Optimum_d")
prval<-list("Optimum Values of Prediction Sum of Square"=preval)
if (!press) val<-val else val<-prval
val
}
