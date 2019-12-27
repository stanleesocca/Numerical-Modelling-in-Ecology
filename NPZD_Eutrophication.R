
NPZDfunc_eutro <- function(t, x, p){
  # function called in NPZD_eutro.m

  T=T_water
  t  #round t to the next integer(t~=n);
  while (n > 365) {
    n = n - 365
  }
  # while (n>365) # count more than one year (for temp, data set only for one year)
  #   n=n-365
  # end
  n=n+(n==0)# day 0 does not exist
  #*****************************************************************
  #               P A R A M E T E R  L I S T
  #-----------------------------------------------------------------   
  a=.063 
  Epp=exp(a*T[n]) # Eppley factor
  # phytoplankton 
  k_N=.05              # half saturation constant
  r_max=1             # maximum  uptake rate of phytoplankton
  lPN=.1                # respiration / extracelluar release
  lPD=.1               # loss through sedimentation
  # zooplankton
  g_max=.5            # maximum grazing rate;
  I=3                      # Ivlev Konstante
  lZN=.3                # in upper layer recycled material
  lZD=.08              # loss rate of zooplankton (involves implicitly "higher" predators,e.g. fish)
  x0=.01                # background concentration
  #----------------   process control -------------------------
  tstart=75              # start of the vernal bloom
  d0=daylength(tstart)   # threshold for the length of the day
  dl=daylength(n)
  theta= dl > d0      #"switch on of "biology"
  A_mix=.5           # winter mixing rate, mixes upper and lower box 
  r_max=r_max*theta*dl*Epp
  k_N=k_N/Epp
  uptake=r_max*x[1]*x[1]/(k_N+x[1]*x[1])
  g_max=g_max*theta*dl*Epp 
  I=I*Epp*1.43*dl 
  grazing=g_max*(1-exp(-I*x[2]*x[2]))
  theti=1-theta      #"switch off" biology
  A_mix=theti*A_mix
  #
  LPN=lPN*uptake*(x[2] >= x0)
  LPD=(theta*lPD*exp(-(n-tstart)/30)+theti*lPD)*(x[2]>x0); # high sinking rates during the spring bloom
  LZN=lZN*grazing*(x[3]>x0)     # egestion and respiration proportiona to ingestion
  LZD=lZD*(1-dl)*(x[3]>x0)      # varying mortality, smaller during longer days (easier escape from predators)
  LD=.001*theta
  # external supply of nutrients and loss through sedimentation (burial)
  Q_ex=Q_ext
  Sed_ex=Sed_loss  # options 1 and 2
  if(OP_ti_on == 2){
    Q_ex=(1+(t<3*365))*Q_ext 
    Sed_ex=(1+(t>3*365))*Sed_loss # Eutrophication and Reversible Eutrophication  
  }
  # if OP_ti_on==2
  # Q_ex=(1+(t<3*365))*Q_ext; Sed_ex=(1+(t>3*365))*Sed_loss;    # Eutrophication and Reversible Eutrophication            
  # end
  #**********************************************
  #              M O D E L  E Q U A T I O N S
  #---------------------------------------------------------------------
  dN =-x[2]*uptake+LZN*x[3]+LPN*x[2]+A_mix*(x[4]-x[1])+LD*x[4]+Q_ex
  dP = x[2]*uptake-(LPN+LPD)*x[2]-grazing*x[3]
  dZ = x[3]*grazing-(LZD+LZN)*x[3]
  dD = LZD*x[3] + LPD*x[2] - A_mix*(x[4]-x[1])-LD*x[4]+Sed_ex;
  # xdot=xdot';
  list(c(dN, dP, dZ, dD))
}



# -------------------------------------------
#NPZD_eutro.m  
#       simple box model of a marine ecosystem similar to NPZDdyn_1
#       with time dependent, rate, temperature effects and external fluxes of nutrients
#       generates Fig. 2.18 and 2.19, for the options 1 and 2
#       options must be toggled 
#===============================================================
  # t0 initial time of simulation
# tf final time

# Prescribe annual cycle of temperature in the surface layer of the Arkona basin (Baltic Sea)
latitude=54
# Prescribe annual cycle of temperature in the surface layer of the Arkona basin (Baltic Sea)
TArkona= c(4, 3.5,  2.5,  2.4,  2,  2,  2,  2,  1.5,  1.5,  2,  2,  2,
           2,  2,  2.5,  3,  4,  5,  6,  7,  8,  9,  10,  11,
           12, 13,  14,  15, 15.5,  15.5,  16,  16.5,  16.5, 16,  15.5,  15,  14,  13,
           12,  11.5,  11,  10.5,  10,  9.5,  9,  8.5,  8,  7,  6,  5, 4.5, 4)
T_water = approx(0:52, TArkona, n = 365)$y   # daily resolution
T=T_water
# --------------------------------------------------------------
  # !!!  TOGGLE A OPTION with '#' symbol !!!
#OP_ti_on=1; disp(' N supply at the rate Q_ext=0.002/d for 3 years')

OP_ti_on = 1
if (OP_ti_on == 1){
  print(' N supply at the rate Q_ext=0.002/d for 3 years')
  tf=3*365
  Q_ext=0.002
  Sed_loss=0   
} else if (OP_ti_on == 2){
  print(' N supply at the rate Q_ext=0.002/d for 3 years and removal at the same rate for the next 3 years')
  tf=6*365
  Q_ext=0.002
  Sed_loss=-Q_ext 
}
#  if OP_ti_on==1 
#          tf=3*365; Q_ext=0.002; Sed_loss=0; 
# end
# if OP_ti_on==2 
#            tf=6*365; Q_ext=0.002; Sed_loss=-Q_ext; 
# end
#---------------------------------------------------------------    
t0=1   
tspan = seq(t0, tf, 1)
# state vector x= [N P Z D]
# initial vector 
x_in=c(N = 0.99, P = .01, Z = .01, D = .99)
#
library(deSolve)
out <- as.data.frame(ode(x_in, tspan, NPZDfunc_eutro, parms = 0))
# [t,x]=ode23('NPZDfunc_eutro',tspan,x_in); # ode- solver



if(OP_ti_on == 1){
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow = c(2, 1), oma = c(0, 0, 3, 0))
  plot(out$time, out$N, type = "l", lwd = 2, xlab = "", ylab = "State variables")
  lines(out$time, out$P, lwd = 2, lty = 2, col = "green")
  lines(out$time, out$Z, lwd = 2, lty = 2, col = "red")
  legend("bottomright", legend = c("N", "P", "Z"), lty = c(1, 2, 2), col = c("black", "red", "green"))
  # Draw the detritus 
  plot(out$time, out$D, lwd = 2, lty = 2, type = "l", col = "blue", xlab = "Time", ylab = "State variables")
  legend("bottomright", legend = c("D"), lty = 2, col = c("blue"))
  mtext("nutrient-phytoplankton-zooplankon cycle in upper box'", side = 3, outer = TRUE)
} else if (OP_ti_on == 2){
  par(mar=c(5.1,4.1,4.1,2.1))
  par(mfrow = c(2, 1), oma = c(0, 0, 3, 0))
  plot(out$time, out$N, type = "l", lwd = 2, xlab = "", ylab = "State variables")
  lines(out$time, out$P, lwd = 2, lty = 2, col = "green")
  lines(out$time, out$Z, lwd = 2, lty = 2, col = "red")
  legend("bottomright", legend = c("N", "P", "Z"), lty = c(1, 2, 2), col = c("black", "red", "green"))
  # Draw the detritus 
  plot(out$time, out$D, lwd = 2, lty = 2, type = "l", col = "blue", xlab = "Time", ylab = "State variables")
  legend("bottomright", legend = c("D"), lty = 2, col = c("blue"))
  mtext("nutrient-phytoplankton-zooplankon cycle in upper box'", side = 3, outer = TRUE)
}

