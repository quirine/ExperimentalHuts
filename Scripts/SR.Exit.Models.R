
hut.movement.model <- function (t, x, params) {
  
  A <- x[1];  B <- x[2];  C <- x[3];   D <- x[4];   E <- x[5]
  X.a <- x[6];  X.b <- x[7];   X.c <- x[8];  X.d <- x[9];    X.e <- x[10];
  K.a <- x[11]; K.b <- x[12];  K.c <- x[13]; K.d <- x[14];   K.e <- x[15];
  q <- params[1:3];  p <- params[4:5];  r <- params[6:8];    k <- params[9:11]
  
  # differential equations
  dA <- (-q[1]-r[1]-k[1])*A + (p[1])*q[2]*B                   # (-qa - ra - ka)A + (pb)qb B             Rs
  dB <- (-q[2]-r[2]-k[2])*B + q[1]*A + (1-p[2])*q[3]*C          # (-qb - rb - kb)B + qa A + (1-pc)qc C         
  dC <- (-q[3]-r[3]-k[3])*C + (1-p[1])*q[2]*B + (1-p[1])*q[2]*D     # (-qc - rc - kc)C + (1-pb) qb B + (1-pb)qb D          
  dD <- (-q[2]-r[2]-k[2])*D + q[1]*E + p[2]*q[3]*C              # (-qb - rb - kb)D + qa E + pc qc C          
  dE <- (-q[1]-r[1]-k[1])*E + p[1]*q[2]*D                       # (-qa - r1 - ka)E + pb  qb D          
  dX.a <- r[1]*A;      dX.b <- r[2]*B;         dX.c <- r[3]*C;       dX.d <- r[2]*D;      dX.e <- r[1]*E;
  dK.a <- k[1]*A;      dK.b <- k[2]*B;         dK.c <- k[3]*C;       dK.d <- k[2]*D;      dK.e <- k[1]*E;
  
  list(c(dA,dB,dC,dD,dE,dX.a,dX.b,dX.c,dX.d,dX.e,dK.a,dK.b,dK.c,dK.d,dK.e))
}

hut.movement.model.ltfu <- function (t, x, params) {
  
  A <- x[1];  B <- x[2];  C <- x[3];   D <- x[4];   E <- x[5]
  X.a <- x[6];  X.b <- x[7];   X.c <- x[8];  X.d <- x[9];    X.e <- x[10];
  K.a <- x[11]; K.b <- x[12];  K.c <- x[13]; K.d <- x[14];   K.e <- x[15];
  U.a <- x[16]; U.b <- x[17];  U.c <- x[18]; U.d <- x[19];   U.e <- x[20];
  q <- params[1:3];  p <- params[4:5];  r <- params[6:8];    k <- params[9:11];   u <- params[12]; 
  
  # differential equations
  dA <- (-q[1]-r[1]-k[1]-u)*A + (p[1])*q[2]*B                   # (-qa - ra - ka)A + (pb)qb B             Rs
  dB <- (-q[2]-r[2]-k[2]-u)*B + q[1]*A + (1-p[2])*q[3]*C          # (-qb - rb - kb)B + qa A + (1-pc)qc C         
  dC <- (-q[3]-r[3]-k[3]-u)*C + (1-p[1])*q[2]*B + (1-p[1])*q[2]*D     # (-qc - rc - kc)C + (1-pb) qb B + (1-pb)qb D          
  dD <- (-q[2]-r[2]-k[2]-u)*D + q[1]*E + p[2]*q[3]*C              # (-qb - rb - kb)D + qa E + pc qc C          
  dE <- (-q[1]-r[1]-k[1]-u)*E + p[1]*q[2]*D                       # (-qa - r1 - ka)E + pb  qb D          
  dX.a <- r[1]*A;      dX.b <- r[2]*B;         dX.c <- r[3]*C;       dX.d <- r[2]*D;      dX.e <- r[1]*E;
  dK.a <- k[1]*A;      dK.b <- k[2]*B;         dK.c <- k[3]*C;       dK.d <- k[2]*D;      dK.e <- k[1]*E;
  dU.a <- u[1]*A;      dU.b <- u[1]*B;         dU.c <- u[1]*C;       dU.d <- u[1]*D;      dU.e <- u[1]*E;
  
  list(c(dA,dB,dC,dD,dE,dX.a,dX.b,dX.c,dX.d,dX.e,dK.a,dK.b,dK.c,dK.d,dK.e,dU.a,dU.b,dU.c,dU.d,dU.e))
}


hut.movement.model.no.abs <- function (t, x, params) {
  A <- x[1]
  B <- x[2]
  C <- x[3]
  D <- x[4]
  E <- x[5]
  q <- params[1:3]
  p <- params[4:5]
  dA <- -q[1]*A + p[1]*q[2]*B                             # -qaA + pb*qb*B
  dB <- -q[2]*B + q[1]*A + q[3]*p[2]*C                    # -qbB + qaA +        qc*pc*C
  dC <- -q[3]*C + q[2]*(1-p[1])*B + q[2]*p[1]*D           # -qcC + qB*(1-pb)B + qD*pD*D
  dD <- -q[2]*D + q[1]*E + q[3]*( 1-p[2] )*C              # -qdD + qeE +        qC*(1-pC)C  
  dE <- -q[1]*E + q[2]*(1-p[1])*D                         # -qeE + qd*(1-pd)Dß
  list(c(dA,dB,dC,dD,dE))
}
