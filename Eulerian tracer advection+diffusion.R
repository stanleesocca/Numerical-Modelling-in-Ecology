## Eulerian tracer advection+diffusion

## 1. Plot the current field ----
  
 # load current field
library(R.matlab)
vortex <- readMat("MEMS_to_R/vortex.mat")
# x, y in m
# U, W in m/d
U <- vortex$U
W <- vortex$W
x <- vortex$x
z <- vortex$z

source("filled-contour3.R")
op <- par(mfcol = c(1, 2))
filled.contour3(x, z, U, main = 'U = From coast', ylab = "Depth", xlab = "Distance")
filled.contour3(x, z, W, main = 'W = From coast', ylab = "Depth", xlab = "Distance")
par(op)

# plot it as a vector field
plot(x, z, ylim = c(200, 0), xlim = c(0, 200), type = "n")
scaling <- 1
arrows(x[row(U)], z[col(U)], x[row(U)] + scaling*U, z[col(U)] + scaling*W, length = .12, angle = 45)


## ---------------------------------------------------
##             DEVELOPING THE 2D MODEL 
# ----------------------------------------------------
## 2. Define domain and initial conditions ----

# NB: To keep Navier-Stokes fluid equilibrium, the `x` and `z` dimensions
#     need to be the same.

# `z` dimension = lines = depth
dz = z[2]-z[1]
nz = dim(z)[2]

# `x` dimension = columns = distance from the coast
dx = x[2]-x[1]
nx = dim(x)[2]

# define the time vector
t0 = 0    # start
tf = 100  # end
# integration time step
dt = 0.05
t = seq(t0, tf, by = dt)
nt = length(t)

# prepare a matrix to store concentrations
C = array(0, dim = c(nz, nx, nt))
# and initialise it with dye in one cell
# (at 70m depth and 70m from the coast)
C[which(z == 70), which(z == 70), 1] = 10

# plot it
filled.contour(x, z, C[,,1])

# diffusivity coefficient
K = 1 # in m2/d here
  
adv_dif_2D = function(C, W, U, dz, dx, z_flag, x_flag, K){
  # vertical 
  if(z_flag == 0){   # In the model domain
    # Advection
    dcz = max(0, -W[1,2]) * C[1,2]/dz + # in from the top
          max(0, W[3,2])* C[3,2]/dz  - # in from the bottom
          abs(W[2,2]) * C[2,2]/dz +   # out from box
      # Diffusion
      K * C[1,2]/dz/dz/2 + # in from top
      K * C[3,2]/dz/dz/2 - # in from bottom
      K * C[2,2]/dz/dz     # out from box
    
  } else if(z_flag == 1){ # top of the model domain
      # Advection
    dcz = max(0, W[3,2]) * C[3,2]/dz -
          max(0, -W[2,2]) * C[2,2]/dz + 
      # Diffusion
      K * C[3,2]/dz/dz/2 -
      K * C[2,2]/dz/dz 
  } else if(z_flag == 2){ # bottom of the model domain
      # Advection
      dcz = max(0, -W[1,2]) * C[1,2]/dz - 
            max(0, W[2,2]) * C[2,2]/dz +
        # Diffusion
        K * C[1,2]/dz/dz/2 -
        K * C[2,2]/dz/dz/2
  } else {
      stop("Unknown position flag!!")
  }
  
  # Horizontal component
  # index can be interpreted as i,j so moving from celll j to cell i
  if(x_flag == 0){
    dcx = max(0, U[2,1]) * C[2,1]/dx + # in from left (from cell 1 into cell 2)
          max(0, -U[2,3]) * C[2,3]/dx -  # in from right (from cell 3 into 2, -U because velocity is moving toward the coast)
          abs(U[2,2]) * C[2,2] +  # out (total movment out from cell 2)
      K * C[2,1]/dx/dx/2 +  # in from left
      K * C[2,3]/dx/dx/2 -  # in from right
      K * C[2,2]/dx/dx   # out
    
  } else if(x_flag == 1){                 # left border
    dcx = max(0, -U[2,3]) * C[2,3]/dx -   #  in from right
          max(0, U[2,2]) * C[2,2]/dx +    # out
      K * C[2,3]/dx/dx/2 +
      K * C[2,2]/dx/dx/2
  } else if(x_flag == 2){                 # right border
    dcx = max(0, U[2,1]) * C[2,1]/dx -    # in from left
          max(0, -U[2,2]) * C[2,2]/dx +   # out
      K * C[2,1]/dx/dx/2 -   # in from left
      K * C[2,2]/dx/dx/2    # out
  } else {
    stop("Unkown flag position!!!")
  }
  
  return(dcz + dcx)
}


## Compute advection+diffusion of tracer ----

# Using Euler forwards for the numerical integration

# for all time steps
for (it in 2:nt) {
  for (iz in 1:nz) {
    for (ix in 1:nx) {
      # Extract a 3x3 neighbourhood of the current point, in which to
      # compute advection+diffusion
      # NB: on the sides, just repeat the first/last points, to keep
      #     the shape 3x3; they won't be used anyways (periodic boundary)
      
      if(iz == 1){ # surface
        z_flag = 1
        z_range = matrix(c(iz, iz, iz+1),1,3)
      } else if(iz == nz){
        z_flag = 2
        z_range = matrix(c(iz-1, iz, iz),1,3)
      } else {
        z_flag = 0
        z_range = matrix(c(iz-1, iz, iz+1),1,3)
      }
      
      if(ix == 1){
        x_flag = 1
        x_range = matrix(c(ix, ix, ix+1),1,3)
      } else if(ix == nx){
        x_flag = 2
        x_range = matrix(c(ix-1, ix, ix),1,3)
      } else {
        x_flag = 0
        x_range = matrix(c(ix-1, ix, ix+1),1,3)
      }
      
      # extract the neighbourhood
      Cc = C[z_range, x_range, it-1]
      Wc = W[z_range, x_range]
      Uc = U[z_range, x_range]
      
      # compute advection+diffusion from center cell
      dc = adv_dif_2D(Cc, Wc,Uc, dz,dx, z_flag,x_flag, K)
      
      # add to current concentration field
      C[iz,ix,it] = C[iz,ix,it-1] + dt*dc
      
    }
  }
  
  # plot once per day
  if((t[it]%% 1) == 0) {
    filled.contour(x, z, C[, , it], main = paste("Time", it), col =topo.colors(20))
  }
}
# Yay!