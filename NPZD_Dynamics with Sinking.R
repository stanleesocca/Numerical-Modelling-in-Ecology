#NPZDdyn_1.m  
#       simple box model of a marine ecosystem
#       with time dependent, rate, temperature effects and external fluxes of nutrients,
#       generates Fig. 2.14
#===============================================================
  # t0 initial time of simulation
# tf final time


NPZDfunc_1 <- function(t, x, p){
      # function corresponding to NPZDdyn_1
      # 
      #*******************************************
      #               P A R A M E T E R  L I S T
      #----------------------------------------------------------------
      
      #global T_water Q_ext Sed_loss 
      T=T_water
      n=round(t)  #round t to the next integer(t~=n)
      a=.063
      Epp=exp(a*T[n]) # Eppley factor
      
      # phytoplankton 
      k_N=.05              # half saturation constant
      r_max=1             # maximum  uptake rate of phytoplankton
      lPN=.1                # respiration / extracelluar release
      lPD=.1                # loss through sedimentation
      
      # zooplankton
      g_max=.5            # maximum grazing rate
      I=3                      # Ivlev Konstante
      lZN=.3                # in upper layer recycled material
      lZD=.08              # loss rate of zooplankton (involves implicitly "higher" predators,e.g. fish)
      x0=.01                # background concentration
      
      #----------------   process control -------------------------
      tstart=75              # start of the vernal bloom
      d0=daylength(tstart)   # threshold for the length of the day
      dl=daylength(t)
      theta=dl>d0      #"switch on of "biology"
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
      LPN=lPN*uptake*(x[2]>=x0)
      LPD=(theta*lPD*exp(-(t-tstart)/30)+theti*lPD)*(x[2]>x0) # high sinking rates during the spring bloom
      LZN=lZN*grazing*(x[3]>x0)     # egestion and respiration proportional to ingestion
      LZD=lZD*(1-dl)*(x[3]>x0)       # varying mortality, smaller during longer days (easier escape from predators)
      LD=0.001*theta                          # vertical diffusive transport
      
      #**********************************************
      #              M O D E L  E Q U A T I O N S
      #---------------------------------------------------------------------
      dN=-x[2]*uptake+LZN*x[3]+LPN*x[2]+A_mix*(x[4]-x[1])+LD*x[4]+Q_ext
      dP= x[2]*uptake-(LPN+LPD)*x[2]-grazing*x[3]
      dZ=x[3]*grazing-(LZD+LZN)*x[3]
      dD=LZD*x[3]+LPD*x[2]-A_mix*(x[4]-x[1])-LD*x[4]+Sed_loss
  list(c(dN, dP, dZ, dD)) 
}



daylength <- function(t){
  # Calculation of the length of the day 
  # based on astronomical formulas
  # ----------------------------------------------
  #
  #global latitude;
  lat=latitude*2*pi/360    
  # Inclination seasonally varying angle between 
  # the equatorial plane and the sun
  # Fourier series (Spencer)
  A=c(.006918, -.399912, -.006758, -.002697)
  B=c(.070257, .000907, .001480)
  #
  # Zenith angle zeta
  #
  G=2*pi/365*t
  inc=A[1]+A[2]*cos(G)+A[3]*cos(2*G)+A[4]*cos(3*G)
  inc=inc+B[1]*sin(G)+B[2]*sin(2*G)+B[3]*sin(3*G)
  #
  # length of the day
  Delta_d=24/pi*acos(-sin(lat)*sin(inc)/cos(lat)/cos(inc)) #in hours
  Delta_d=Delta_d/24 # normalized
  
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
  A=c(.006918, -.399912, -.006758, -.002697)
  B=c(.070257, .000907, .001480)
  # Zenith angle zeta
  G=2*pi/365*t
  inc=A[1]+A[2]*cos(G)+A[3]*cos(2*G)+A[4]*cos(3*G)
  inc=inc+B[1]*sin(G)+B[2]*sin(2*G)+B[3]*sin(3*G)
  # length of the day
  Delta_d=24/pi*acos(-sin(lat)*sin(inc)/cos(lat)/cos(inc)) #in hours
  Delta_d/24 # normalized
}


## ---------------------------

# global latitude T_water Q_ext Sed_loss
latitude=54
# Prescribe annual cycle of temperature in the surface layer of the Arkona basin (Baltic Sea)
TArkona= c(4, 3.5,  2.5,  2.4,  2,  2,  2,  2,  1.5,  1.5,  2,  2,  2,
         2,  2,  2.5,  3,  4,  5,  6,  7,  8,  9,  10,  11,
         12, 13,  14,  15, 15.5,  15.5,  16,  16.5,  16.5, 16,  15.5,  15,  14,  13,
         12,  11.5,  11,  10.5,  10,  9.5,  9,  8.5,  8,  7,  6,  5, 4.5, 4)
# T_water=spline(0:52,TArkona,0:1/7:52) # daily resolution
T_water = approx(0:52, TArkona, n = 365)$y
T=T_water
# external supply of nutrients and loss through sedimentation (burial)
Q_ext=0.005 
Sed_loss=-0.005 
t0=1   
tf=365 
tspan = seq(t0, tf, 1)
#tspan=c(t0, tf)
# state vector x= [N P Z D]
# initial vector x_in
x_in=c(N = 0.99, P = .01, Z = .01, D = .99)

library(deSolve)
# ---  ordinary differential equation solver ---
out = as.data.frame(ode(x_in, tspan, NPZDfunc_1, parms = 0))

#[t,x]=ode23('NPZDfunc_1',tspan,x_in) 
# visualize the results 
par(mar=c(5.1,4.1,4.1,2.1))
par(mfrow = c(2, 1), oma = c(0, 0, 3, 0))
plot(out$time, out$N, type = "l", lwd = 2, xlab = "", ylab = "State variables")
lines(out$time, out$P, lwd = 2, lty = 2, col = "green")
lines(out$time, out$Z, lwd = 2, lty = 2, col = "red")
legend("bottomright", legend = c("N", "P", "Z"), lty = c(1, 2, 2), col = c("black", "red", "green"))
# Draw the detritus 
plot(out$time, out$D, lwd = 2, lty = 2, type = "l", col = "blue", xlab = "Time", ylab = "State variables")
legend("bottomright", legend = c("D"), lty = 2, col = c("blue"))
mtext("NPZD Model", side = 3, outer = TRUE)

       