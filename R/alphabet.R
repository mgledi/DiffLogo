
angle1 = seq(0.3 + pi/2, pi, length = 100)
angle2 = seq(pi, 1.5 * pi, length = 100)

alphabet$A = list();
alphabet$A$x = c(0,  4,  6, 10,  8,  5, 2, 0, NA,2.2,2.6,7.4,7.8,2.2) * 0.1
alphabet$A$y = c(0, 10, 10,  0,  0,7.5, 0, 0, NA,  3,  4,  4,  3,  3) * 0.1

alphabet$B = list();
alphabet$B$x = c(0, 2, 2,0,NA) * 0.1
alphabet$B$y = c(10,10,0,0,NA) * 0.1

alphabet$E = list();
alphabet$E$x = c(0,  10, 10,  2,   2, 0, 0,  NA, 2, 9, 9,   2,   NA, 2,10,10  ,2  ) * 0.1
alphabet$E$y = c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5, NA, 0,0,1.5,1.5) * 0.1

alphabet$F = list();
alphabet$F$x = c(0,  10, 10,  2,   2, 0, 0,  NA, 2, 8, 8, 2 ) * 0.1
alphabet$F$y = c(10, 10, 8.5, 8.5, 0, 0, 10, NA, 4, 4, 5.5, 5.5 ) * 0.1

alphabet$H = list();
alphabet$H$x = c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA,0,10,10,0) * 0.1
alphabet$H$y = c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA,4,4,6,6) * 0.1

alphabet$I = list();
alphabet$I$x = c(3,  7,  7, 3) * 0.1
alphabet$I$y = c(10, 10, 0, 0) * 0.1

alphabet$K = list();
alphabet$K$x = c(0,  2,  2, 0, NA, 0,8, 10, 0, NA, 5, 10, 8, 3) * 0.1
alphabet$K$y = c(10, 10, 0, 0, NA, 4,10,10,2, NA, 6, 0 , 0, 6) * 0.1

alphabet$L = list();
alphabet$L$x = c(0,  2,  2,   10,  10, 0) * 0.1
alphabet$L$y = c(10, 10, 1.5, 1.5, 0,  0 ) * 0.1

alphabet$M = list();
alphabet$M$x = c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, 1.5,4,6.0,3.5 , NA, 8.5,6,4.,6.5) * 0.1
alphabet$M$y = c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10 ,0,0,  10, NA, 10, 0,  0,10) * 0.1

alphabet$N = list();
alphabet$N$x = c(0,  2,  2, 0, NA, 8, 10, 10, 8, NA, .5,2.5, 9.5,7.5 ) * 0.1
alphabet$N$y = c(10, 10, 0, 0, NA, 10, 10, 0, 0, NA, 10, 10,  0,  0) * 0.1

alphabet$O = list();
alphabet$O$x = c( (0.5 + 0.5 * sin(angle1)),(0.5 + 0.5 * sin(angle2)) ) 
alphabet$O$y = c( (0.5 + 0.5 * cos(angle1)),(0.5 + 0.5 * cos(angle2)) )

alphabet$T = list();
alphabet$T$x = c(0, 10, 10, 6, 6, 4, 4, 0) * 0.1
alphabet$T$y = c(10, 10, 9, 9, 0, 0, 9, 9) * 0.1

alphabet$V = list();
alphabet$V$x = c(0,4,6,10,8,5,2) * 0.1
alphabet$V$y = c(10,0,0,10,10,2,10) * 0.1

alphabet$W = list();
alphabet$W$x = c(0, 2,4,5.5, 4, 3.0,1.5, NA,4.5, 6,8,10, 8.5,7,6) * 0.1
alphabet$W$y = c(10,0,0,10, 10,  2,  10, NA,10,0,0,10,10,2,10) * 0.1

alphabet$X = list();
alphabet$X$x = c(0,2,10,8,NA,0,2,10,8) * 0.1
alphabet$X$y = c(10,10,0,0,NA,0,0,10,10) * 0.1

alphabet$Y = list();
alphabet$Y$x = c(0,2,6,4,NA,4,6,10,8,NA,4,6,6,4) * 0.1
alphabet$Y$y = c(10,10,4.5,4.5,NA,4.5,4.5,10,10,NA,5,5,0,0) * 0.1

alphabet$Z = list();
alphabet$Z$x = c(0,2.5,10,7.5,NA,0,10,10,0,NA,0,10,10,0) * 0.1
alphabet$Z$y = c(1.5,1.5,8.5,8.5,NA,10,10,8.5,8.5,NA,0,0,1.5,1.5) * 0.1

plot(0,0,xlim=c(-1,11),ylim=c(-1,11)); 
text(0,0,"O",adj=c(0,0),cex=33,font=2,family="sans")
polygon(alphabet$O$x*10,alphabet$O$y*10,col="green",border="green")
