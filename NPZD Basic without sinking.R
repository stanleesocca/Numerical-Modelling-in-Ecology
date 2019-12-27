## SIMPLE MODEL 
# Translation of Wolfgang matlab code to R
# -----------------------------------------
#       SIMPLE NPZD Model
# -----------------------------------------
NPZD = function(t, x, p){
  # with(as.list(c(t, x)), {
  #*******************************************
    #               P A R A M E T E R  L I S T
  #-----------------------------------------------------------------
    # phytoplankton 
  k_N=.1                # half saturation constant
  r_max=1             # maximum  uptake rate of phytoplankton
  LPN=.01             # respiration / extracelluar release
  LPD=.02             # loss through sedimentation
  
  # zooplankton
  g_max=.5            # maximum grazing rate;
  I=1.2                   # Ivlev Konstante
  LZN=.01             # in upper layer recycled material
  LZD=.02             # loss rate of zooplankton [involves implicitly "higher" predators,e.g. fish]
  x0=.01                 # background level
  
  #----------------   process control -------------------------
    d0=daylength54(75)  # threshold for the length of the day
  dl=daylength54(t)
  theta=dl>d0            #"switch on of "biology"
A_mix=.5                  # winter mixing 
    r_max=r_max*theta
    uptake=r_max*x[1]*x[1]/(k_N+x[1]*x[1])
    g_max=g_max*theta
    grazing=g_max*(1-exp(-I*x[2]*x[2]))
    theti=1-theta      #"switch off " biology
A_mix=theti*A_mix
#*********************************************************************
#              M O D E L  E Q U A T I O N S
#---------------------------------------------------------------------
dN =-x[2]*uptake+LZN*(x[3]-x0)+LPN*(x[2]-x0)+A_mix*(x[4]-x[1]);
dP= x[2]*uptake-(LPN+LPD)*(x[2]-x0)-x[3]*grazing;
dZ=x[3]*grazing-(LZD+LZN)*(x[3]-x0);
dD=LZD*(x[3]-x0)+LPD*(x[2]-x0)-A_mix*(x[4]-x[1]);
return(list(c(dN, dP, dZ, dD)))
# })
}


daylength54 = function(t){
    # Calculation of the length of the day 
    # for the latitude 54° N/S,
    # based on astronomical formulas
    # ----------------------------------------------
      latitude=54
      lat=latitude*2*pi/360   # Latitude
      # Inclination seasonally varying angle between 
      # the equatorial plane and the sun
      # Fourier series (Spencer)
      A=c(.006918, -.399912, -.006758, -.002697);
      B=c(.070257, .000907, .001480);
      # Zenith angle zeta
      G=2*pi/365*t;
      inc=A[1]+A[2]*cos(G)+A[3]*cos(2*G)+A[4]*cos(3*G);
      inc=inc+B[1]*sin(G)+B[2]*sin(2*G)+B[3]*sin(3*G);
      # length of the day
      Delta_d=24/pi*acos(-sin(lat)*sin(inc)/cos(lat)/cos(inc)) #in hours
      Delta_d/24 # normalized
}

# NPZDdyn.m 
#        simple box model of a marine ecosystem
#        generates Fig. 2.11
#==============================================
# t0 initial time of simulation
# tf final time
# state vector [N P Z D]
# initial vector x_in

t0=0;tf=365;
tspan = seq(t0, tf)
x_in=c(N = 0.99, P = .01, Z = .01, D = .99)
# call ordinary differential equation solver
library(deSolve)
out = as.data.frame(ode(y= x_in, times =  tspan, func =  NPZD, parms = 0))

par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow = c(2, 1), oma = c(0, 0, 3, 0))
plot(out$time, out$N, type = "l", col = "red", ylab = "NPZ", xlab = "Time", lwd = 2)
lines(out$time, out$P, col = "green", lwd = 2)
lines(out$time, out$Z, col = "blue", lwd = 2)
plot(out$time, out$D, ylab = "Detritus", xlab = "Time", type = "l", lwd = 2)
mtext("NPZD Model", side = 3, outer = TRUE)

