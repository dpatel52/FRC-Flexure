# main_script.py
from Mk_Functions import zone111, zone211,zone212 ,zone221 ,zone222 ,zone311,zone312 ,zone321 ,zone322 ,zone411,zone412 ,zone421 ,zone422,zone4222    # Import specific functions
try:
    from plotBetaVsMandK import plotBetaVsMandK
except ImportError:
    import sys
    sys.path.append('path_to_directory_containing_plotBetaVsMandK')
    from plotBetaVsMandK import plotBetaVsMandK
import numpy as np

from Envelope_Calc import calculate_envelope_new_2 # type: ignore

from deflection_functions import calculate_deflection, calculate_deflection_3PB,plot_deflection # type: ignore

from draw_doubly_reinforced_beam import draw_doubly_reinforced_beam

import matplotlib.pyplot as plt

hasFlex = 0  # Set to 1 if experimental data is available, otherwise set to 0
exp_load_data = []
exp_deflection_data = []

if hasFlex == 1:
    try:
        data = pd.read_excel('U3-15M-D.xlsx')
        exp_load_data = data.iloc[:, 1].to_numpy()      # Second column for load data
        exp_deflection_data = data.iloc[:, 0].to_numpy()  # First column for deflection data
    except Exception as e:
        print("Error: Could not read the Excel file. Please check the file and try again.")
        print(e)
else:
    print("No experimental data will be plotted.")


# ======================= Geometry and Loading ===========================

pointBend = 3  # Number of point loads applied
L = 1000  # Length of the beam (in inches)
b = 120  # Width of the beam (in inches)
h = 175  # Height of the beam (in inches)
S2 = 125  # Spacing between Loading Points (4PB) only
Lp = 125  # Length of the span (in inches)
PostLp = Lp  # Half of the height (in inches)
cover = 25  # Concrete cover (in inches)
alpha = (h - cover) / h  # Ratio of effective height

# ==================== Tension Model Parameters ==========================

E = 34000  # Modulus of elasticity (psi)
epsilon_cr = 0.00012  # Cracking strain
mu_1 = 0.44  # Multiplier for post-crack strength at beta 1
mu_2 = 2.5   # Multiplier for post-crack strength at beta 2
mu_3 = 0.33  # Multiplier for post-crack strength at beta 3
beta_1 = 3  # Strain parameter for zone 1
beta_2 = 30  # Strain parameter for zone 2
beta_3 = 290  # Strain parameter for zone 3
eta_1 = (mu_1 - 1) / (beta_1 - 1)  # Slope for zone 1
eta_2 = (mu_2 - mu_1) / (beta_2 - beta_1)  # Slope for zone 2
eta_3 = (mu_3 - mu_2) / (beta_3 - beta_2)  # Slope for zone 3

# =================== Compression Model Parameters ======================


xi = 1.01  # Ex/E ratio (>1.01)
omega = 7.3  # Transition point for compression model
mu_c = 1  # Multiplier for compression
ecu = 0.003  # Ultimate strain in compression
lambda_cu = ecu / epsilon_cr  # Ultimate strain in compression normalized
eta_c = (mu_c - 1) / (lambda_cu - omega)  # Slope for compression model

# ======================== Steel Properties =============================

Es = 200000  # Modulus of elasticity of steel
n = Es / E  # Es/E ratio
kappa = 0.002 / epsilon_cr  # Nominal yield strain esy/ecr
mu_s = 1  # Multiplier for steel
chi_su = 6 * kappa  # Ultimate strain in steel
eta_s = (mu_s - 1) / (chi_su - kappa)  # Slope for steel model

# ======================= Reinforcement Details =========================

topDiameter = 10  # Diameter of top bars (inches) - #3 bar
topCount = 2  # Number of top bars
botDiameter = 12  # Diameter of bottom bars (inches) - #4 bar
botCount = 2  # Number of bottom bars

# ====================== Reinforcement Areas ============================

s_area_bot = botCount * (botDiameter ** 2 * np.pi / 4)  # Bottom steel area
rho_t = s_area_bot / (b * h)  # Reinforcement ratio in tension
s_area_top = topCount * (topDiameter ** 2 * np.pi / 4)  # Top steel area
rho_c = s_area_top / (b * h)  # Reinforcement ratio in compression

# ======================== Plots =============================

# Tension Model (Matrix)
strainT = [0, epsilon_cr, epsilon_cr * beta_1, epsilon_cr * beta_2, epsilon_cr * beta_3]
stressT = [0, epsilon_cr * E, mu_1 * epsilon_cr * E, mu_2 * epsilon_cr * E, mu_3 * epsilon_cr * E]


# Compression Model (Matrix)
strainC = [0, omega * epsilon_cr, epsilon_cr * lambda_cu]
stressC = [0, epsilon_cr * E * xi * omega, epsilon_cr * E * xi * omega * mu_c]

# Reinforcement Model
strainR = [0, kappa * epsilon_cr, epsilon_cr * chi_su]
stressR = [0, epsilon_cr * E * n * kappa, epsilon_cr * E * n * kappa * mu_s]

'''
# Plot Tension Model
plt.figure()
plt.plot(strainT, stressT, '-o', linewidth=2, color='r', label='Tension Model')
plt.title('Tension Model (Matrix)', fontsize=14, fontweight='bold')
plt.xlabel('Strain (mm/mm)', fontsize=14, fontweight='bold')
plt.ylabel('Stress (MPa)', fontsize=14, fontweight='bold')
plt.grid(True)
plt.legend(loc='best')
plt.show()

# Plot Compression Model
plt.figure()
plt.plot(strainC, stressC, '-o', linewidth=2, color='b')
plt.title('Compression Model (Matrix)', fontsize=14, fontweight='bold')
plt.xlabel('Strain (mm/mm)', fontsize=14, fontweight='bold')
plt.ylabel('Stress (MPa)', fontsize=14, fontweight='bold')
plt.grid(True)
plt.show()

# Plot Reinforcement Model
plt.figure()
plt.plot(strainR, stressR, '-o', linewidth=1, color='g')
plt.title('Reinforcement Model')
plt.xlabel('Strain')
plt.ylabel('Stress')
plt.grid(True)
plt.show()
'''

# ======================== Cracking Properties =============================

# Obtain k from zone111 function
kcr, Mcr = zone111(1, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)


# Calculate cracking moment, curvature, and flexural rigidity
M_cr = (epsilon_cr * E * b * h**2) / (12 * (1 - kcr))
Phi_cr = epsilon_cr / ((1 - kcr) * h)
EI_cr = M_cr / Phi_cr

'''

# Display the results
print("Cracking Moment (M_cr):", M_cr)
print("Cracking Curvature (Phi_cr):", Phi_cr)
print("Flexural Rigidity (EI_cr):", EI_cr)

'''

# Beta Arrays Initialization
beta_z1 = np.linspace(0, 1, 1000)
beta_z2 = np.linspace(1, beta_1, 2000)
beta_z3 = np.linspace(beta_1, beta_2, 2000)
beta_z4 = np.linspace(beta_2, beta_3, 2000)

# ======================== Envelope Calculation =============================

k111, M111 = zone111(beta_z1, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k211, M211 = zone211(beta_z2, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k212, M212 = zone212(beta_z2, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k221, M221 = zone221(beta_z2, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k222, M222 = zone222(beta_z2, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k311, M311 = zone311(beta_z3, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k312, M312 = zone312(beta_z3, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k321, M321 = zone321(beta_z3, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k322, M322 = zone322(beta_z3, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k411, M411 = zone411(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k412, M412 = zone412(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k421, M421 = zone421(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k422, M422 = zone422(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)
k4222, M4222 = zone4222(beta_z4, L, b, h, alpha, E, epsilon_cr, beta_1, beta_2, beta_3, eta_1, eta_2, eta_3, xi, omega, eta_c, n, kappa, eta_s, rho_c, rho_t)


# Initialize zone indicators
T = 1
C = 1
R = 1
RC = 1

# Process each beta zone sequentially. Notice that the updated T, C, R, RC
# are passed from one call to the next.
Envelope1, T, C, R, RC = calculate_envelope_new_2(
    rho_c,rho_t,kappa, omega, epsilon_cr, beta_z1,
    k111, M111,
    k211, M211, k212, M212,
    k221, M221, k222, M222,
    k311, M311, k312, M312,
    k321, M321, k322, M322,
    k411, M411, k412, M412,
    k421, M421, k422, M422,
    k4222, M4222,
    beta_1, beta_2, beta_3, alpha,
    T, C, R, RC
)

Envelope2, T, C, R, RC = calculate_envelope_new_2(
    rho_c,rho_t,kappa, omega, epsilon_cr, beta_z2,
    k111, M111,
    k211, M211, k212, M212,
    k221, M221, k222, M222,
    k311, M311, k312, M312,
    k321, M321, k322, M322,
    k411, M411, k412, M412,
    k421, M421, k422, M422,
    k4222, M4222,
    beta_1, beta_2, beta_3, alpha,
    T, C, R, RC
)

Envelope3, T, C, R, RC = calculate_envelope_new_2(
    rho_c,rho_t,kappa, omega, epsilon_cr, beta_z3,
    k111, M111,
    k211, M211, k212, M212,
    k221, M221, k222, M222,
    k311, M311, k312, M312,
    k321, M321, k322, M322,
    k411, M411, k412, M412,
    k421, M421, k422, M422,
    k4222, M4222,
    beta_1, beta_2, beta_3, alpha,
    T, C, R, RC
)

Envelope4, T, C, R, RC = calculate_envelope_new_2(
    rho_c,rho_t,kappa, omega, epsilon_cr, beta_z4,
    k111, M111,
    k211, M211, k212, M212,
    k221, M221, k222, M222,
    k311, M311, k312, M312,
    k321, M321, k322, M322,
    k411, M411, k412, M412,
    k421, M421, k422, M422,
    k4222, M4222,
    beta_1, beta_2, beta_3, alpha,
    T, C, R, RC
)

# Concatenate the four envelopes vertically
Envelope = np.vstack((Envelope1, Envelope2, Envelope3, Envelope4))
print(Envelope)
# Combine all beta arrays into one long vector (if needed)
beta_all = np.concatenate((beta_z1, beta_z2, beta_z3, beta_z4))

# Combine beta_all and Envelope into one array
Envelope_beta = np.column_stack((beta_all, Envelope))

#print(Envelope)

strain_bot = beta_all * epsilon_cr
NA_from_bot = (1 - Envelope[:, 0]) * h
Curvature = strain_bot / NA_from_bot
Moment = Envelope[:, 1]

plt.figure(figsize=(10, 6))
plt.plot(Curvature, Moment, color='darkblue', linewidth=2, label='Moment-Curvature')
plt.xlabel('Curvature', fontsize=14, fontweight='bold')
plt.ylabel('Moment', fontsize=14, fontweight='bold')
plt.title('Moment vs. Curvature', fontsize=16, fontweight='bold')
plt.grid(True)
plt.legend(fontsize=12)
plt.tight_layout()
plt.show()

mu = (mu_1 + mu_2 + mu_3) / 3  # A factor to adjust the number of fibers

# Call the function without an axis
draw_doubly_reinforced_beam(h, b, cover, topDiameter, topCount, botDiameter, botCount, mu)



plt.figure()
# In Python, array indexing is zero-based so the second column is index 1
plt.plot(beta_all, Envelope[:, 1], '-', linewidth=2, color='k')
plt.xlabel('beta')   # Optional: Label for the x-axis
plt.ylabel('Moment')  # Optional: Label for the y-axis
plt.title('Beta vs Moment')  # Optional: Title for the plot
plt.show()

# Example plotting code with markers and legend
plt.figure()

# Zone 1
plt.plot(beta_z1, M111, 'ko-', linewidth=2, markersize=6, label='Zone111')

# Zone 2 (four curves with different markers)
plt.plot(beta_z2, M211, 'r^-', linewidth=2, markersize=6, label='Zone211')
plt.plot(beta_z2, M212, 'rs-', linewidth=2, markersize=6, label='Zone212')
plt.plot(beta_z2, M221, 'rv-', linewidth=2, markersize=6, label='Zone221')
plt.plot(beta_z2, M222, 'rd-', linewidth=2, markersize=6, label='Zone222')

# Zone 3 (four curves with different markers)
plt.plot(beta_z3, M311, 'go-', linewidth=2, markersize=6, label='Zone311')
plt.plot(beta_z3, M312, 'g^-', linewidth=2, markersize=6, label='Zone312')
plt.plot(beta_z3, M321, 'gs-', linewidth=2, markersize=6, label='Zone321')
plt.plot(beta_z3, M322, 'gv-', linewidth=2, markersize=6, label='Zone322')

# Zone 4 (four curves with blue markers and one with cyan)
plt.plot(beta_z4, M411, 'bo-', linewidth=2, markersize=6, label='Zone411')
plt.plot(beta_z4, M412, 'b^-', linewidth=2, markersize=6, label='Zone412')
plt.plot(beta_z4, M421, 'bs-', linewidth=2, markersize=6, label='Zone421')
plt.plot(beta_z4, M422, 'bd-', linewidth=2, markersize=6, label='Zone422')
plt.plot(beta_z4, M4222, 'c*-', linewidth=2, markersize=8, label='Zone4222')
#plt.plot(beta_all, Envelope[:, 1], '-', linewidth=2, color='k')

plt.title('Moment vs Beta', fontsize=14, fontweight='bold')
plt.xlabel('Beta', fontsize=14, fontweight='bold')
plt.ylabel('Moment M', fontsize=14, fontweight='bold')
plt.grid(True)
plt.legend()  # Display the legend with labels
plt.show()


# Example plotting code with markers and legend
plt.figure()

# Zone 1
plt.plot(beta_z1, k111, 'ko-', linewidth=2, markersize=6, label='Zone111')

# Zone 2 (four curves with different markers)
plt.plot(beta_z2, k211, 'r^-', linewidth=2, markersize=6, label='Zone211')
plt.plot(beta_z2, k212, 'rs-', linewidth=2, markersize=6, label='Zone212')
plt.plot(beta_z2, k221, 'rv-', linewidth=2, markersize=6, label='Zone221')
plt.plot(beta_z2, k222, 'rd-', linewidth=2, markersize=6, label='Zone222')

# Zone 3 (four curves with different markers)
plt.plot(beta_z3, k311, 'go-', linewidth=2, markersize=6, label='Zone311')
plt.plot(beta_z3, k312, 'g^-', linewidth=2, markersize=6, label='Zone312')
plt.plot(beta_z3, k321, 'gs-', linewidth=2, markersize=6, label='Zone321')
plt.plot(beta_z3, k322, 'gv-', linewidth=2, markersize=6, label='Zone322')

# Zone 4 (four curves with blue markers and one with cyan)
plt.plot(beta_z4, k411, 'bo-', linewidth=2, markersize=6, label='Zone411')
plt.plot(beta_z4, k412, 'b^-', linewidth=2, markersize=6, label='Zone412')
plt.plot(beta_z4, k421, 'bs-', linewidth=2, markersize=6, label='Zone421')
plt.plot(beta_z4, k422, 'bd-', linewidth=2, markersize=6, label='Zone422')
plt.plot(beta_z4, k4222, 'c*-', linewidth=2, markersize=8, label='Zone4222')

# Plot settings
plt.title('K vs Beta', fontsize=14, fontweight='bold')
plt.xlabel('Beta', fontsize=14, fontweight='bold')
plt.ylabel('K', fontsize=14, fontweight='bold')
plt.grid(True)
plt.legend()  # Display the legend with labels
plt.show()

# --------------------- Calculate Load -----------------------
if pointBend == 4:
    S1 = (L - S2) / 2
    load_4pb = (Envelope[:, 1] * 2) / S1
elif pointBend == 3:
    load_4pb = (Envelope[:, 1] * 4) / L
else:
    load_4pb = (Envelope[:, 1] * 4) / L

# -------------- Moment–Area Method Calculations --------------
# MC_data is a stacked array: [beta_all, Envelope(:,2), Curvature]
MC_data = np.column_stack((beta_all, Envelope[:, 1], Curvature))
Mmax = np.max(MC_data[:, 1])
idx_array = np.where(MC_data[:, 1] == Mmax)[0]
if idx_array.size > 0:
    idx = idx_array[0]
    Cmax = MC_data[idx, 2]
else:
    raise ValueError("Mmax not found in MC_data")

mom = MC_data[:, 1]
cv = MC_data[:, 2]

# --------- Calculate deflection based on bending type ---------
if pointBend == 4:
    delta_total = calculate_deflection(mom, cv, M_cr, Cmax, Phi_cr,
                                        L, Lp, PostLp, Mmax, S2)
elif pointBend == 3:
    delta_total = calculate_deflection_3PB(mom, cv, M_cr, Cmax, Phi_cr,
                                            L, Lp, PostLp, Mmax, S2)
else:
    delta_total = calculate_deflection(mom, cv, M_cr, Cmax, Phi_cr,
                                        L, Lp, PostLp, Mmax, S2)

# ---------------- Plot the Load-Deflection Curve ----------------
plot_deflection(kappa, epsilon_cr, chi_su, omega, lambda_cu,
                rho_t, rho_c, beta_all, delta_total, load_4pb,
                alpha, Envelope, hasFlex, exp_deflection_data, exp_load_data)


