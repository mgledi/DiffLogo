
Letter = function(x,y) {
   if(length(x) != length(y)) {
        stop("The length of vector x is different from length of vector y")
   }
   pts = list(x=x,y=y);
   class(pts) = "Letter";
   return(pts);
}
letters=list();
############## A
letters$A = Letter(
  c(0,  4,  6, 10,  8,  5, 2, 0, NA,2.2,2.6,7.4,7.8,2.2) * 0.1,
  c(0, 10, 10,  0,  0,7.5, 0, 0, NA,  3,  4,  4,  3,  3) * 0.1
)
############## P
tmpx.1=sin(seq(0, 1*pi, length = 80))*5 +5.0
tmpy.1=cos(seq(0, 1*pi, length = 80))*2.875 +7.125
tmpx.2=rev(sin(seq(0, 1*pi, length = 80))*3 +5.0)
tmpy.2=rev(cos(seq(0, 1*pi, length = 80))*1.375 +7.125)
letters$P = Letter(
  c(0, 2, 2,0,NA,0, 5, tmpx.1,5,0,0,5,tmpx.2,5,0)*0.1,
  c(10,10,0,0,NA,10,10,tmpy.1,4.25,4.25,5.75,5.75,tmpy.2,8.5,8.5)*0.1
)
#  c(0, 2, 2,0,NA,0, 5, 5,  0,  NA,0,   5,   5,   0,   NA,tmpx.1,tmpx.2)*0.1,
#  c(10,10,0,0,NA,10,10,8.5,8.5,NA,5.75,5.75,4.25,4.25,NA,tmpy.1,tmpy.2)*0.1

############## B
letters$B = Letter(
  c(0, 2, 2,0,NA,0, 5, 5,  0,  NA,0,   5,   5,   0,   NA,0,5,5,  0,  NA,tmpx.1,tmpx.2,NA,tmpx.1,tmpx.2)*0.1,
  c(10,10,0,0,NA,10,10,8.5,8.5,NA,5.75,5.75,4.25,4.25,NA,0,0,1.5,1.5,NA,tmpy.1,tmpy.2,NA,tmpy.1-4.25,tmpy.2-4.25)*0.1
)
############## R
letters$R = Letter(
  c(0, 2, 2,0,NA,0, 5, 5,  0,  NA,0,   5,   5,   0,   NA,5, 10, 8, 3,NA,tmpx.1,tmpx.2)*0.1,
  c(10,10,0,0,NA,10,10,8.5,8.5,NA,5.75,5.75,4.25,4.25,NA,5, 0 , 0, 5,NA,tmpy.1,tmpy.2)*0.1
)
############## C is copied from SeqLogo
angle1 = seq(0.3 + pi/2, pi, length = 40)
angle2 = seq(pi, 1.5 * pi, length = 40)
x.l1 = 0.5 + 0.5 * sin(angle1)
y.l1 = 0.5 + 0.5 * cos(angle1)
x.l2 = 0.5 + 0.5 * sin(angle2)
y.l2 = 0.5 + 0.5 * cos(angle2)
x.l = c(x.l1, x.l2)
y.l = c(y.l1, y.l2)
x = c(x.l, rev(x.l))
y = c(y.l, 1 - rev(y.l))
x.i1 = 0.5 + 0.35 * sin(angle1)
y.i1 = 0.5 + 0.35 * cos(angle1)
x.i1 = x.i1[y.i1 <= max(y.l1)]
y.i1 = y.i1[y.i1 <= max(y.l1)]
y.i1[1] = max(y.l1)
x.i2 = 0.5 + 0.35 * sin(angle2)
y.i2 = 0.5 + 0.35 * cos(angle2)
x.i = c(x.i1, x.i2)
y.i = c(y.i1, y.i2)
letters$C = Letter(
  c(x, rev(c(x.i, rev(x.i))))/max(x),
  c(y, rev(c(y.i, 1 - rev(y.i))))
)
############## G is copied from SeqLogo
letters$G = Letter(
  c(C$x, NA, 1,  0.5,0.5,0.8,0.8,1,1),
  c(C$y, NA, 0.4,0.4,0.3,0.3,0,  0,0.4)
)
############## D
tmpx.1=sin(seq(0, 1*pi, length = 80))*5 +5.0
tmpy.1=cos(seq(0, 1*pi, length = 80))*5 +5
tmpx.2=rev(sin(seq(0, 1*pi, length = 40))*3.5 +5.0)
tmpy.2=rev(cos(seq(0, 1*pi, length = 40))*3.5 +5.0)
letters$D = Letter(
  c(0, 2, 2,0,NA,0, 5, 5,  0,  NA,0,5,5, 0,   NA,tmpx.1,tmpx.2)*0.1,
  c(10,10,0,0,NA,10,10,8.5,8.5,NA,0,0,1.5,1.5,NA,tmpy.1,tmpy.2)*0.1
)
############## E
letters$E = Letter(
  c(0,  10, 10,  2,   2, 0, 0,  NA, 2, 9, 9,   2,   NA, 2,10,10  ,2  ) * 0.1,
  c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5, NA, 0,0,1.5,1.5) * 0.1
)
############## F
letters$F = Letter(
  c(0,  10, 10,  2,   2, 0, 0,  NA, 2, 8, 8, 2 ) * 0.1,
  c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5 ) * 0.1
)
############## H
letters$H = Letter(
  c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA,0,10,10,0) * 0.1,
  c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA,4,4,6,6) * 0.1
)
############## H
letters$I = Letter(
  c(3,  7,  7, 3) * 0.1,
  c(10, 10, 0, 0) * 0.1
)
############## J
tmpx.1=sin(seq(0.5*pi, 1.5*pi, length = 40))*3 +5
tmpy.1=cos(seq(0.5*pi, 1.5*pi, length = 40))*2.5 +2.5
tmpx.2=rev(sin(seq(.5*pi, 1.5*pi, length = 40))*1 +5)
tmpy.2=rev(cos(seq(.5*pi, 1.5*pi, length = 40))*1 +2.5)
letters$J = Letter(
  c(6,  8 , 8, 6,NA,tmpx.1,tmpx.2) * 0.1,
  c(10, 10, 2.5, 2.5,NA,tmpy.1,tmpy.2) * 0.1
)
############## K
letters$K = Letter(
	c(0,  2,  2, 0, NA, 0,8, 10, 0, NA, 5, 10, 8, 3) * 0.1,
	c(10, 10, 0, 0, NA, 4,10,10,2, NA, 6, 0 , 0, 6) * 0.1
)
############## L
letters$L = Letter(
  c(0,  2,  2,   10,  10, 0) * 0.1,
  c(10, 10, 1.5, 1.5, 0,  0 ) * 0.1
)
############## M
letters$M = Letter(
  c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, 1.5,4,6.0,3.5 , NA, 8.5,6,4.,6.5) * 0.1,
  c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10 ,1,1,  10, NA, 10, 1,  1,10) * 0.1
)
############## N
letters$N = Letter(
  c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, .5,2.5, 9.5,7.5 ) * 0.1,
  c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10, 10,  0,  0) * 0.1
)
############## O
tmpx=sin(seq(0, 2*pi, length = 80))/2+0.5
tmpy=cos(seq(0, 2*pi, length = 80))/2+0.5
letters$O = Letter(
  c(tmpx,rev(tmpx)*0.6+0.2),
  c(tmpy,tmpy*0.7+0.15)
)
############## Q
letters$Q = Letter(
  c(letters$O$x,NA,.5,.7,1.0,.8),
  c(letters$O$y,NA,.3,.3,.0,.0)
)
############## S
tmpx.1=(sin(seq(0.0*pi, 1.5*pi, length = 100))*5 +5)
tmpy.1=cos(seq(0.0*pi, 1.5*pi, length = 100))*2.825 +2.825
tmpx.2=rev(sin(seq(.0*pi, 1.5*pi, length = 100)))*3 +5
tmpy.2=rev(cos(seq(.0*pi, 1.5*pi, length = 100))*1.5 +2.825)

letters$S = Letter(
  c(tmpx.1,tmpx.2,NA,-tmpx.1+10,-tmpx.2+10) * 0.1,
  c(tmpy.1,tmpy.2,NA,-tmpy.1+10,-tmpy.2+10) * 0.1
)

############## T
letters$T = Letter(
  c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1,
  c(10, 10, 9, 9, 0, 0, 9, 9) * 0.1
)
############## U
tmpx.1=rev(sin(seq(0.5*pi, 1.5*pi, length = 80))*5 +5)
tmpy.1=cos(seq(0.5*pi, 1.5*pi, length = 80))*3 +3
tmpx.2=(sin(seq(.5*pi, 1.5*pi, length = 80))*3 +5)
tmpy.2=rev(cos(seq(.5*pi, 1.5*pi, length = 80))*1.75 +3)
letters$U = Letter(
  c(0,  0,tmpx.1,10, 10,8 ,8,tmpx.2,2,2 ) * 0.1,
  c(10, 3,tmpy.1,3, 10,10,3,tmpy.2,3,10) * 0.1
)
############## V
letters$V = Letter(
  c(0,4,6,10,8,5,2) * 0.1,
  c(10,0,0,10,10,2,10) * 0.1
)
############## W
letters$W = Letter(
  c(0, 2,4,5.5, 4, 3.0,1.5, NA,4.5, 6,8,10, 8.5,7,6) * 0.1,
  c(10,0,0,10, 10,  2,  10, NA,10,0,0,10,10,2,10) * 0.1
)
############## X
letters$X = Letter(
  c(0,2,10,8,NA,0,2,10,8) * 0.1,
  c(10,10,0,0,NA,0,0,10,10) * 0.1
)
############## Y
letters$Y = Letter(
  c(0,2,6,4,NA,4,6,10,8,NA,4,6,6,4) * 0.1,
  c(10,10,4.5,4.5,NA,4.5,4.5,10,10,NA,5,5,0,0) * 0.1
)
############## Z
letters$Z = Letter(
  c(0,2.5,10,7.5,NA,0,10,10,0,NA,0,10,10,0) * 0.1,
  c(1.5,1.5,8.5,8.5,NA,10,10,8.5,8.5,NA,0,0,1.5,1.5) * 0.1
)

pdf("tst.pdf");
par(mfcol=c(5,5),mar=c(0,0,0,0))
for(a in c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y")) {
  plot(NA,0,xlim=c(0,10),ylim=c(0,10),xaxt="n",yaxt="n"); 
  #text(0,0,a,adj=c(0,0),cex=10,font=2,family="sans")
  polygon(letters[[a]]$x*10,letters[[a]]$y*10,col="black",border="black")
}
dev.off();


alphabet = function(chars,cols) {
   lets = list();
   for ( i in 1:length(chars)) {
     lets[[chars[i]]] = letters[[chars[i]]];
     lets[[chars[i]]]$col = cols[i];
   }
   obj = list(chars=chars,cols=cols,letters=lets)
   class(obj)="alphabet"
   return(obj);
}

DNA = alphabet(c("A","C","G","T"),c("green4","blue","orange","red"));
