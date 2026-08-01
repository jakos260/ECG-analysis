"""
HodgkinHuxley.py  :   Hodgkin-Huxley model of the ventricular action potential.

Calculation of a Hodgkin-Huxley action potential.
The Euler method is used to integrate the first order non-linear
differential equations that describe the transmembrane potential and
the gating variables.
     Im = Cm(dVm/dt) + gNa*(Vm-ENa) + gK*(Vm-EK) + gL*(Vm-EL)
     gNa = gNa_max * m^3 * h
     gK  = gK_max * n^4
     gL  = gL
     dm/dt = am*(1-m) - bm*m
     dh/dt = ah*(1-h) - bh*h
     dn/dt = ah*(1-n) - bn*n
where the rate coefficients am, ah, an, bm, bh, bn are given by nonlinear
functions of the membrane potential.

AUTHOR
    Johannes J. Struijk
    Cardiotechnology Research Group
    Dept. Health Science and Technology
    Aalborg University
    Denmark
    jjs@hst.aau.dk

DATE
   14 March 2022
"""

import numpy as np
from matplotlib import pyplot as plt 

# Initialization

# Parameters for numerical integration
TotalTime = 10.0e-3            # Total time (s)
delta_t = 1.0e-6               # Integration time step (s)
N = round(TotalTime/delta_t)+1 # Number of time steps

# Stimulus current
IAmpl = 0.2              # Stimulus current (uA/cm^2)
PD = 0.2e-3              # Stimulus pulse duration (s)
N_on = round(PD/delta_t) # Number of time steps that the current is on
Istim = np.zeros(N)      # Initialize stimulus current
Istim[0:N_on-1] = IAmpl  # Stim. pulse is on during first N_on time steps

# Membrane parameters
Cm = 1.0e-6      # Membrane capacitance (F/cm^2)
gNa_max = 120e-3 # Maxiumum sodium conductance (S/cm^2)
gK_max = 36e-3   # Maximum potassium conductance (S/cm^2)
gL = 0.3e-3      # Leakage conductance (S/cm^2)
ENa = 50         # Sodium Nernst potential(mV)
EK = -77         # Potassium Nernst potential(mV)
EL = -54.4       # Leakage Nernst potential(mV)
Vrest = -65      # Resting membrane potential (mV)

# Calculate membrane potential (Note: all voltages are in mV, time is in s)
Vm = np.zeros(N) # Reserve memory for membrane potential
Vm[0] = Vrest    # Start with the resting membrane potential
m = np.zeros(N)  # Reserve memory for the sodium activation variable m
m[0] = (2.5/(np.exp(2.5)-1)) / (4+(2.5/(np.exp(2.5)-1)))  # Initialize m
h = np.zeros(N)  # Reserve memory for the sodium inactivation variable h
h[0] = 0.07/(0.07+(1/(np.exp(3)+1)))                      # Initialize h
n = np.zeros(N)   # Reserve memory for the potassium activation variable n
n[0] = (0.1/(np.exp(1)-1)) / (0.125+(0.1/(np.exp(1)-1)))  # Initialize n
gNa = gNa_max*(m[0]**3)*h[0]
gK = gK_max*n[0]**4

for i in range(0,N-1):
    # Update the membrane potential:
    dVmdt = (Istim[i]-gNa*(Vm[i]-ENa)-gK*(Vm[i]-EK)-gL*(Vm[i]-EL))/Cm
    Vm[i+1] = Vm[i] + delta_t*dVmdt
    # Update the ionic conductances:
    # - calculate the transfer rate coefficients 
    Vp = Vm[i]-Vrest;  # Vp is used eight times in this loop
    am = 0.1*(25-Vp)/(np.exp(0.1*(25-Vp))-1)
    bm = 4/np.exp(Vp/18)
    ah = 0.07/np.exp(Vp/20)
    bh = 1/(np.exp(0.1*(30-Vp))+1)
    an = 0.01*(10-Vp)/(np.exp(0.1*(10-Vp))-1)
    bn = 0.125/np.exp(Vp/80)
    # - update the gating variables and conductances
    dmdt = am*(1-m[i]) - bm*m[i] # in 1/ms, so a factor 1000 is missing
    dhdt = ah*(1-h[i]) - bh*h[i]
    dndt = an*(1-n[i]) - bn*n[i]
    m[i+1] = m[i] + delta_t*dmdt*1e3 # correction for the factor 1000
    h[i+1] = h[i] + delta_t*dhdt*1e3
    n[i+1] = n[i] + delta_t*dndt*1e3
    gNa = gNa_max*(m[i+1]**3)*h[i+1]
    gK = gK_max*n[i+1]**4


# Plot
t = np.linspace(0,(N-1)*delta_t,N) * 1000; # ms
fig1, axes = plt.subplots(nrows=4, ncols=1, sharex=True)
axes[0].plot(t,Vm)
axes[0].set_ylabel("Vm (mV)")
axes[1].plot(t,m)
axes[1].set_ylabel("m (Na activation)")
axes[2].plot(t,h)
axes[2].set_ylabel("h (Na inactivation)")
axes[3].plot(t,n)
axes[3].set_ylabel("n (K activation)")
axes[3].set_xlabel("time (ms)")
plt.show()

# EOF
