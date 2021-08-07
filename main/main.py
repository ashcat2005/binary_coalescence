# Binary Coalescence

import numpy as np
import matplotlib.pyplot as plt


# Orbital parameters
def semilatus(L):
  '''
  Returns the semilatus rectum of the orbit 
  from the angular momentum L
  '''
  return L**2/GM


def eccentricity(E, L):
  '''
  Returns the eccentricity of the orbit from the 
  energy E and the angular momentum L
  '''
  return np.sqrt(1 + 2*E*L**2/(GM**2))


def phi(L, r, ang0):
  '''
  Returns the true anomaly obtained from E, L and r
  '''
  ang1 = ang0 + dt*L/r**2
  return ang1

def ellipse(p, ecc, f):
  '''
  Returns the coordinates of the point and 
  the radius at angle f
  '''
  r = p/(1-ecc*np.cos(f - omega))
  x = r*np.cos(f)
  y = r*np.sin(f)
  z = 0.
  return x, y, z, r

# Quadrupole tensor components
def qij(x):
  '''
  Returns the components of the 3x3 quadrupole tensor
  '''
  return mu*(np.outer(x,x) - np.identity(3)*np.dot(x,x)/3)

def second_derivative(f,h):
  '''
  Second derivative of a function
  '''
  return (11*f[0] - 56*f[1] + 114*f[2] - 104*f[3] + 35*f[4])/(12*h**2)


def third_derivative(f,h):
  '''
  Third derivative of a function
  '''
  return (3*f[0] - 14*f[1] + 24*f[2] - 18*f[3] + 5*f[4])/(2*h**3)

# Time evolution
def advance_time(i):
  global dt
  global r_crit
  ecc = eccentricity(E,L)
  p = semilatus(L)
  # Update the state variables
  time[i+1] = time[i] + dt
  angle[i+1] = phi(L, position[i,3], angle[i])
  position[i+1] = ellipse(p, ecc, angle[i+1])
  q[i+1] = qij(position[i,0:3])

  # Check the separation radius to change the time-step
  if position[i+1,3]>1E-1:
    dt = 1E-3 # time-step back to the original
    r_crit = 1E-1 # critical radius back to the original
  elif position[i+1,3]<r_crit:
    dt = 0.01*dt # diminish the time-step
    r_crit = 0.01*r_crit # diminish the critical radius
  elif position[i+1,3]>r_crit:
    dt = 100*dt # increase the time-step
    r_crit = 100*r_crit # increase the critical radius


G = 4*np.pi**2 # Gravitational constant
c = 63197.8 # Speed of light in units of au/yr

# Masses  
m1 = 10. # Solar masses
m2 = 7.  # Solar masses

M = m1 + m2 # Total mass (solar masses)
mu = m1*m2/M # Reduced mass
GM = G*M # gravitational parameter

# Radii
r1 = 2*G*m1/c**2 # Schwarzschild radius for m1
r2 = 2*G*m2/c**2 # Schwarzschild radius for m2
r_merge = r1 + r2 # Separation distance for fusion/collision

# Initial Values
E = -300. # energy
L = 15. # angular momentum
omega = np.pi/3 # argument of the pericenter


# Time grid definition
n = 500000 # number of steps
time = np.zeros(n)
dt = 1E-3 # timestep

position = np.zeros([n,4]) # Position
angle = np.zeros(n) # True Anomaly
q = np.zeros([n,3,3]) # Quadrupole tensor
d3qdt = np.zeros([3,3]) # Third derivative of the Quadrupole tensor
d2qdt = np.zeros([3,3]) # Second derivative of the Quadrupole tensor
h_plus = np.zeros(n) # h_+ polarization component of the GW
h_cross = np.zeros(n) # h_x polarization component of the GW

# Initial condition in the grid
angle[0] = 0.
position[0] = ellipse(semilatus(L), eccentricity(E,L), angle[0])
q[0] = qij(position[0,0:3])

# Critical radius to modify dt
r_crit = 1E-1


# First four steps without lose of E and L
for i in range(4):
  advance_time(i)


# Main Loop
for i in range(4,n-1):
  # Check for the merge of the binary system
  if position[i,3]<r_merge:
    print('Binary system merged after {:.5f} years using {:.0f} steps'.format(time[i],i))
    print('\n The final separation distance is {:.5e} au.\n'.format(position[-1,3]))
    break

  advance_time(i)
  # Third derivative of the quadrupole tensor
  d3qdt = third_derivative(q[i-4:i+2],dt)
  dE = (G/(5*c**5))*(d3qdt[0,0]**2 + d3qdt[1,1]**2 + d3qdt[2,2]**2 +\
                     2*d3qdt[0,1]**2 + 2*d3qdt[0,2]**2 + 2*d3qdt[1,2]**2)
  # Second derivative of the quadrupole tensor
  d2qdt = second_derivative(q[i-4:i+2],dt)
  dL = (2*G/(5*c**5))*(d2qdt[0,0]*d3qdt[1,0] + d2qdt[0,1]*d3qdt[1,1] +\
                       d2qdt[0,2]*d3qdt[1,2] - d2qdt[1,0]*d3qdt[0,0] -\
                       d2qdt[1,1]*d3qdt[0,1] - d2qdt[1,2]*d3qdt[0,2])
  h_plus[i] = (d2qdt[0,0] - d2qdt[1,1])/position[i,3]
  h_cross[i] = 2*d2qdt[0,1]/position[i,3]
  
  #Energy and Angular Momentum change (Amplified version)
  L = L - 1e8*dL*dt
  E = E - 1e8*dE*dt
  if L < 0:
    L = 0.



plt.figure(figsize=(10,7))
plt.plot(position[:,0],position[:,1])
plt.axhline(color='black',alpha=0.3)
plt.axvline(color='black',alpha=0.3)
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.title('Binary System Coalescence')
plt.savefig('BinarySystemCoalescence.jpg')
plt.show()