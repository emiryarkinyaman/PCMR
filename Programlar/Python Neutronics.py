# # # # 0.0 INFO

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

# # # # 1.0 LIBRARY IMPORT
import math

# # # # 2.0 VARIABLE DEFINITIONS
# Universal Variables
N_a = 6.0221e23 # Number
N_a_barn = round(N_a*1e-24,4) # Barn Optimized Ratio

# Volume Specific Variables
fuel_volume_ratio_in_one_lattice = 0.4 # Ratio
moderator_volume_ratio_in_one_lattice = 0.4 # Ratio
coolant_volume_ratio_in_one_lattice = 0.2 # Ratio # Method: Coolant Displacement. Either control rod or the coolant tubes will be functional.
rod_volume_ratio_in_one_lattice = 0.1 # Ratio # Method: Coolant Displacement. Either control rod or the coolant tubes will be functional.

# Control Rod Specific Variables
rod_insertion_ratio = 0.74 # Ratio
rod_count = 7 # Number

# Fuel Specific Variables
fuel_enrichment = 0.1975 # Ratio
fuel_packing_factor = 0.4000 # Ratio

# Fuel Specific Industry Standards
p_resonance = 0.8300 # Ratio
fast_fission_factor = 1.04 # Ratio
nu = 2.4300 # Ratio

# Cube Core Specific Variables
one_lattice_width = 10.0 # cm.
outer_reflector_expansion = 6.5 # cm.

active_core_height = 110.0 # cm.
active_core_width = active_core_height # cm.
active_core_volume = round(active_core_height**3,1) # cm^3

total_core_height = round(active_core_height+(2*outer_reflector_expansion),1) # cm.
total_core_width = total_core_height # cm.
total_core_volume = round(total_core_height**3,1) # cm^3
one_lattice_volume = round((one_lattice_width*one_lattice_width*active_core_height),6) # cm^3

lattice_count = int((active_core_height*active_core_width)//(one_lattice_width*one_lattice_width)) # Number

# Volume Calculations
total_fuel_volume = round(fuel_volume_ratio_in_one_lattice*lattice_count*one_lattice_volume,1) # cm^3
total_moderator_volume = round(moderator_volume_ratio_in_one_lattice*lattice_count*one_lattice_volume,1) # cm^3
total_coolant_volume = round(coolant_volume_ratio_in_one_lattice*(lattice_count-rod_count)*one_lattice_volume,1) # cm^3
total_rod_volume = round(rod_volume_ratio_in_one_lattice*rod_count*one_lattice_volume,1) # cm^3

fuel_volume_ratio = round(total_fuel_volume/active_core_volume,6) # Ratio
moderator_volume_ratio = round(total_moderator_volume/active_core_volume,6) # Ratio
coolant_volume_ratio = round(total_coolant_volume/active_core_volume,6) # Ratio
rod_volume_ratio = round((total_rod_volume*rod_insertion_ratio)/active_core_volume,6) # Ratio
vacuum_volume_ratio = round((active_core_volume-total_fuel_volume-total_moderator_volume-total_coolant_volume-total_rod_volume)/active_core_volume,6) # Ratio
reflector_volume_ratio = round((total_core_volume-active_core_volume)/total_core_volume,6) # Ratio

# # # # 3.0 DENSITY CALCULATIONS
# B-10
B_abundance = 0.199 # Ratio
B4C_molar_mass = round((4*(10.81))+(1*(12.01)),6) # g/mol

B4C_density_1 = 2.520 # g/cm^3
B4C_density_2 = ((B4C_density_1*N_a_barn)/B4C_molar_mass) # molecules/barn-cm

B10 = round((B4C_density_2*4)*B_abundance*rod_volume_ratio,6) # molecules/barn-cm

# K-39
K_molar_mass = 39.10 # g/mol
K_density_1 = 0.740 # g/cm^3
K39 = round(((K_density_1*N_a_barn)/K_molar_mass)*coolant_volume_ratio*0.9326,6)

# U-235 & U-238 & C-12 & O-16
UCO_molar_mass = round(238.0+(0.5*12.01)+(1.5*16.0),6) # g/mol

UCO_density_1 = 10.50 # g/cm^3
UCO_density_2 = ((UCO_density_1*N_a_barn)/UCO_molar_mass) # molecules/barn-cm

fuel_density = UCO_density_2*fuel_volume_ratio*fuel_packing_factor # molecules/barn-cm

U235 = round(fuel_density*fuel_enrichment,6) # molecules/barn-cm
U238 = round(fuel_density*(1-fuel_enrichment),6) # molecules/barn-cm
O16  = round(fuel_density*1.5,6) # molecules/barn-cm
C12 = round((B4C_density_2*1*rod_volume_ratio)+(fuel_density*0.5),6) # molecules/barn-cm

# H-1 & Zr-90
ZrH2_molar_mass = round(91.22+(1.7*1.008),3) # g/mol

ZrH2_density_1 = 5.610 # g/cm^3
ZrH2_density_2 = ((ZrH2_density_1*N_a_barn)/ZrH2_molar_mass) # molecules/barn-cm

Zr90 = round(ZrH2_density_2*1.0*moderator_volume_ratio,6) # molecules/barn-cm
H1   = round(ZrH2_density_2*1.7*moderator_volume_ratio,6) # molecules/barn-cm

atom_densities = {
    "B-10":B10,
    "K-39":K39,
    "U-238":U238,
    "U-235":U235,
    "C-12":C12,
    "O-16": O16,
    "H-1":H1,
    "Zr-90":Zr90
}

cross_sections = {
    "B-10":{"T_c":3840.00000,"T_f":0.0,"T_s":2.3,"F_c":0.50,"F_f":0.0,"F_s":1.4},
    "K-39":{"T_c":2.10000,"T_f":0.0,"T_s":1.4,"F_c":0.00,"F_f":0.0,"F_s":3.0},
    "U-238":{"T_c":2.70000,"T_f":0.0,"T_s":10.0,"F_c":0.07,"F_f":0.3,"F_s":5.0},
    "U-235":{"T_c":100.00000,"T_f":583.0,"T_s":10.0,"F_c":0.09,"F_f":1.0,"F_s":4.0},
    "C-12":{"T_c":0.00340,"T_f":0.0,"T_s":4.8,"F_c":0.00,"F_f":0.0,"F_s":4.8},
    "O-16":{"T_c":0.00019,"T_f":0.0,"T_s":3.7,"F_c":0.00,"F_f":0.0,"F_s":3.7},
    "H-1":{"T_c":0.34000,"T_f":0.0,"T_s":80.0,"F_c":0.00,"F_f":0.0,"F_s":4.0},
    "Zr-90":{"T_c":0.01400,"T_f":0.0,"T_s":6.4,"F_c":0.00,"F_f":0.0,"F_s":6.0}
}

energy_decrements = {"B-10":0.187,"K-39":0.050,"U-238":0.008,"U-235":0.008,"C-12":0.157,"O-16":0.120,"H-1":1.000,"Zr-90":0.023}

# # # # 4.0 THEORETICAL FISSION
macro_sigma = {
    "Thermal Absorption":0.0,
    "Thermal Fission":0.0,
    "Thermal Scatter":0.0,
    "Fast Absorption":0.0,
    "Fast Fission":0.0,
    "Fast Scatter":0.0,
    "Fuel Thermal Absorption":0.0,
    "Energy Decrement":0.0
}

for iso, density in atom_densities.items():
    if iso in cross_sections:
        T_absorption = cross_sections[iso]["T_c"]+cross_sections[iso]["T_f"]
        F_absorption = cross_sections[iso]["F_c"]+cross_sections[iso]["F_f"]

        macro_sigma["Thermal Absorption"] += density*T_absorption
        macro_sigma["Thermal Fission"] += density*cross_sections[iso]["T_f"]
        macro_sigma["Thermal Scatter"] += density*cross_sections[iso]["T_s"]

        macro_sigma["Fast Absorption"] += density*F_absorption
        macro_sigma["Fast Fission"] += density*cross_sections[iso]["F_f"]
        macro_sigma["Fast Scatter"] += density*cross_sections[iso]["F_s"]

        if iso in ["U-235","U-238"]:
            macro_sigma["Fuel Thermal Absorption"] += density*T_absorption

        macro_sigma["Energy Decrement"] += density*cross_sections[iso]["F_s"]*energy_decrements[iso]

# # # #

B = round(3*((math.pi/active_core_width)**2),6)

thermal_utilization = round(macro_sigma["Fuel Thermal Absorption"]/macro_sigma["Thermal Absorption"],6)
reproduction = round((nu*atom_densities["U-235"]*cross_sections["U-235"]["T_f"])/macro_sigma["Fuel Thermal Absorption"],6)

thermal_diffusion = round(1.0/(3.0*(macro_sigma["Thermal Absorption"]+macro_sigma["Thermal Scatter"])),6)
thermal_diffusion_length_squared = round(thermal_diffusion/macro_sigma["Thermal Absorption"],6)

fast_diffusion = round(1.0/(3.0*(macro_sigma["Fast Absorption"]+macro_sigma["Fast Scatter"])),4)
slowing_down_spectra = round(fast_diffusion/macro_sigma["Energy Decrement"],5)

migration_area = round(thermal_diffusion_length_squared+slowing_down_spectra,5)

fast_phi = round(math.exp(-B*slowing_down_spectra),6)
thermal_phi = round(1.0/(1.0+(B*thermal_diffusion_length_squared)),6)

K_INF = round(thermal_utilization*reproduction*p_resonance*fast_fission_factor,6)
K = round(K_INF*fast_phi*thermal_phi,6)

RHO = round((K-1.0)/K,6)
PCM = round(RHO*100000,2)

# # # # 5.0 OUTPUT I
print("Microreactor: Specs")
print(f"Fuel: {fuel_packing_factor*100}% Compressed {fuel_enrichment*100}% Enriched UCO (TRISO)\nModerator: ZrH")
print(f"Cold Start: Am-Be\nCoolant: K\nControl Rod: B4C\nReflector: BeO")
print(f"Rod Situation: {rod_insertion_ratio*100}% inserted total of {rod_count} control rods.")
print(f"Total Core Dimensions: {total_core_height} x {total_core_height} x {total_core_width} cm.")
print(f"Total Core Volume: {total_core_volume*0.001} L.")
print(f"Active Core Dimensions: {active_core_height} x {active_core_height} x {active_core_width} cm.")
print(f"Active Core Volume: {active_core_volume*0.001} L.")
print(f"Vacuum/Active Volume: {round(vacuum_volume_ratio*100,6)}%")
print(f"Fuel/Active Volume: {round(fuel_volume_ratio*100,5)}%")
print(f"Moderator/Active Volume: {round(moderator_volume_ratio*100,5)}%")
print(f"Coolant/Active Volume: {round(coolant_volume_ratio*100,5)}%")
print(f"Rod/Active Volume: {round(rod_volume_ratio*100,6)}%")
print(f"Reflector/Total Volume: {round(reflector_volume_ratio*100,6)}%")
print()
print("Microreactor: Analysis")
print(f"Selected Rod Configuration: {rod_insertion_ratio*100}% inserted total of {rod_count} control rods.")
# print(f"Geometric Buckling (B^2): {B} cm^-2")
print(f"Thermal Diffusion Length Squared: {thermal_diffusion_length_squared} cm^2")
print(f"Slowing Down Spectra: {slowing_down_spectra} cm^2")
print(f"Migration Area: {migration_area} cm^2")
print(f"Thermal Phi: {thermal_phi}")
print(f"Fast Phi: {fast_phi}")
print(f"Thermal Utilization: {thermal_utilization}")
print(f"Reproduction: {reproduction}")
print(f"K Inf: {K_INF}")
print(f"K: {K}")
print(f"Reactivity: {PCM} pcm.")

# # # # 6.0 NEUTRON FLUX (RK4)
beta_eff = 0.0064
am_be_neutron_flux = 1e10

prompt_generation = 1e-4
prompt_decay = 0.08

dt = 1e-3

def get_diff_eq(phi,c,rho_val):
    d_phi_dt = ((rho_val-beta_eff)/prompt_generation)*phi+prompt_decay*c+am_be_neutron_flux
    d_c_dt = (beta_eff/prompt_generation)*phi-prompt_decay*c
    return d_phi_dt,d_c_dt

rho_rise_per_second = 0.000050

initial_rho = RHO
initial_pcm = PCM

current_rho = initial_rho
current_pcm = initial_pcm
current_phi = 0
current_c = 0
current_t = 0

print(f"\nMicroreactor: Cold Start | t={current_t} secs.")
print(f"Initial Reactivity: {current_pcm} pcm.")
print(f"Reactivity Rise Per Second: {round(rho_rise_per_second*100000,6)} pcm.")
print(f"dt: {dt} s.")
print(f"Am-Be Source: {am_be_neutron_flux} n/cm^2-s")

criticality_counter = 0
while current_rho < 410.6/1e5:
    k1_phi,k1_c = get_diff_eq(current_phi,current_c,current_rho)
    k2_phi,k2_c = get_diff_eq(current_phi+0.5*dt*k1_phi,current_c+0.5*dt*k1_c,current_rho)
    k3_phi,k3_c = get_diff_eq(current_phi+0.5*dt*k2_phi,current_c+0.5*dt*k2_c,current_rho)
    k4_phi,k4_c = get_diff_eq(current_phi+dt*k3_phi,current_c+dt*k3_c,current_rho)

    current_phi += (dt/6.0)*(k1_phi+2*k2_phi+2*k3_phi+k4_phi)
    current_c += (dt/6.0)*(k1_c+2*k2_c+2*k3_c+k4_c)

    current_t += dt
    current_rho = initial_rho+(rho_rise_per_second*current_t)
    criticality_counter += 1

# # # #

brayton_eff = 0.43
total_macro_fission = macro_sigma["Thermal Fission"]+macro_sigma["Fast Fission"]

MeV_per_fission = 200
J_per_fission = MeV_per_fission*1.6e-13
capacity_factor = 0.97

MWt = active_core_volume*total_macro_fission*current_phi*J_per_fission/1e6
MWe = MWt*brayton_eff*capacity_factor

print(f"\nMicroreactor: Delayed Supercritical | t={round(current_t,4)} secs.")
print(f"Final Reactivity: {round(current_rho*100000,4)} pcm.")
print(f"Final Flux: {round(current_phi/1e11,1)}e11 n/cm^2-s.")
print(f"Expected Thermal Power: {round(MWt,1)} MWt")
print(f"Expected Electrical Output: {round(MWe,1)} MWe")