# INFO

"""
Potassium Heat Pipe Cooled Microreactor
Emir Yarkın Yaman, 2200386009, Hacettepe University Nuclear Engineering Dep.

Fuel: UCO (HALEU TRISO)
Cold Start: Am-Be
Coolant: K
Control Rod: B4C
Reflector: BeO
Moderator: ZrH
"""

# LIBRARIES
import numpy as np
import matplotlib.pyplot as plt

# CONSTANT VARIABLES
beta = 0.0064
lambda_decay = 0.080
lambda_generation = 1e-3

# EN2TECH PROJECT SPECIFIC VARIABLES
am_be_source = 1e10

initial_rho = -914.4e-5

dt = 1e-3
t_max = 265.0
reactivity_rise_per_second = 5.0e-5

n_0 = (-am_be_source*lambda_generation)/initial_rho
c_0 = (beta*n_0)/(lambda_decay*lambda_generation)

# TIME INTERVALS
t_interval_1 = 0
t_interval_2 = t_max
t_interval_3 = 375
t_interval_4 = 500

final_rho = initial_rho+reactivity_rise_per_second*(2*t_max-t_interval_3)
# CONSTANT FUNCTIONS
def get_reactivity(t0):
    if t_interval_1 <= t0 < t_interval_2:
        return initial_rho+(reactivity_rise_per_second*t0)
    elif t_interval_2 <= t0 < t_interval_3:
        return initial_rho+(reactivity_rise_per_second*t0)*(-1.0)
    else:
        return final_rho

def get_point_kinetics(t0,y,beta,lambda_decay,lambda_generation):
    n_solution,c_solution = y
    reactivity = get_reactivity(t0)

    dn_dt = ((reactivity-beta)/lambda_generation)*n_solution+(lambda_decay*c_solution)+am_be_source
    dc_dt = (beta/lambda_generation)*n_solution-(lambda_decay*c_solution)

    return np.array([dn_dt,dc_dt])

# RUNGE-KUTTA SOLUTION
def RK(t1,dt):
    float_to_int_conversion_2 = int(t1/dt)+1
    time_dots = np.linspace(0,t1,float_to_int_conversion_2)

    n_list = [n_0]
    c_list = [c_0]

    Y = np.array([n_0,c_0])

    for t in time_dots[:-1]:
        K1 = get_point_kinetics(t,Y,beta,lambda_decay,lambda_generation)
        K2 = get_point_kinetics(t+0.5*dt,Y+0.5*dt*K1,beta,lambda_decay,lambda_generation)
        K3 = get_point_kinetics(t+0.5*dt,Y+0.5*dt*K2,beta,lambda_decay,lambda_generation)
        K4 = get_point_kinetics(t+dt,Y+dt*K3,beta,lambda_decay,lambda_generation)

        Y_ = Y+(dt/6.0)*(K1+2*K2+2*K3+K4)

        n_list.append(Y_[0])
        c_list.append(Y_[1])

        Y = Y_

    return time_dots,np.array(n_list),np.array(c_list)

# SOLUTION OUTPUT
t_final = t_interval_4
t_runge,n_runge,c_runge = RK(t_final,dt)

plt.style.use("dark_background")
plt.figure(figsize=(12,8),dpi=100)

am_be_source_1 = f"{round(am_be_source/1e10,1)}e10 n/cm^3.s"
title = f"Initial Reactivity: {initial_rho*1e5} pcm.\nReactivity Rise/Drop Per Second: {int(reactivity_rise_per_second*1e5)} pcm."

# PLOT I: NEUTRON DENSITY
plt.subplot(2,1,1)
plt.plot(t_runge,n_runge,linewidth=1.67,linestyle="-",color="#FF00FF",alpha=0.67)
plt.ylabel("Neutron Density n(t) n/cm^3")
plt.title(f"Runge-Kutta Method On Microreactor's Cold Start — n(t)\nTimestep: {dt} secs. | Max. Time: {t_final} secs.\nAm-Be Source: {am_be_source_1}\n{title}")
plt.grid(True,which="both",linestyle="--",alpha=0.7)

# PLOT II: PRECURSOR CONCENTRATION
plt.subplot(2,1,2)
plt.plot(t_runge,c_runge,linewidth=1.67,linestyle="-",color="#00FFFF",alpha=0.67)
plt.ylabel("Precursor Concentration c(t) a/cm^3")
plt.xlabel("Time (s)")

plt.title(f"Runge-Kutta Method On Microreactor's Cold Start — c(t)")
plt.grid(True,which="both",linestyle="--",alpha=0.7)

plt.tight_layout()
plt.show()
plt.close()