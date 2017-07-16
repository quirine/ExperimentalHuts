
load('ExitData.ltfu.RData')

head(Data.end)               

Ind = which(Data$Dose==0)
kd.0 = sum(Data$k.kd.A[Ind]) + sum(Data$k.kd.B[Ind]) + sum(Data$k.kd.C[Ind]) + sum(Data$k.kd.D[Ind]) + sum(Data$k.kd.E[Ind])
Ind = which(Data.end$Dose==0) 
kd.0 = kd.0 + sum(Data.end$k.kd.A.end[Ind]) + sum(Data.end$k.kd.B.end[Ind]) + sum(Data.end$k.kd.C.end[Ind]) + sum(Data.end$k.kd.D.end[Ind]) + sum(Data.end$k.kd.E.end[Ind])

Ind = which(Data$Dose==0.0625)
kd.l = sum(Data$k.kd.A[Ind]) + sum(Data$k.kd.B[Ind]) + sum(Data$k.kd.C[Ind]) + sum(Data$k.kd.D[Ind]) + sum(Data$k.kd.E[Ind])
Ind = which(Data.end$Dose==0.0625) 
kd.l = kd.l + sum(Data.end$k.kd.A.end[Ind]) + sum(Data.end$k.kd.B.end[Ind]) + sum(Data.end$k.kd.C.end[Ind]) + sum(Data.end$k.kd.D.end[Ind]) + sum(Data.end$k.kd.E.end[Ind])

Ind = which(Data$Dose==0.125)
kd.h = sum(Data$k.kd.A[Ind]) + sum(Data$k.kd.B[Ind]) + sum(Data$k.kd.C[Ind]) + sum(Data$k.kd.D[Ind]) + sum(Data$k.kd.E[Ind])
Ind = which(Data.end$Dose==0.125) 
kd.h = kd.h + sum(Data.end$k.kd.A.end[Ind]) + sum(Data.end$k.kd.B.end[Ind]) + sum(Data.end$k.kd.C.end[Ind]) + sum(Data.end$k.kd.D.end[Ind]) + sum(Data.end$k.kd.E.end[Ind])
