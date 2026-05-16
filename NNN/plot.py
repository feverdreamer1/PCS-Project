import numpy as np 
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==============================================================================
# GLOBAL STYLE (paper-ready)
# ==============================================================================
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "legend.fontsize": 8,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.dpi": 120
})

# ==============================================================================
# 1. ONSAGER EXACT SOLUTION (2D ISING)
# ==============================================================================
Tc_onsager = 2.0 / np.log(1.0 + np.sqrt(2.0))

def mag_onsager(T):
    if T < Tc_onsager:
        return (1.0 - np.sinh(2.0 / T)**(-4))**(1.0 / 8.0)
    return 0.0

T_teorica = np.linspace(1.5, 6.0, 500)
M_teorica = [mag_onsager(t) for t in T_teorica]


# ==============================================================================
# 2. DATA LOADING
# ==============================================================================
def load_data(file):
    data = np.loadtxt(file)
    return {
        'N': data[:, 0],
        'T': data[:, 1],
        'M': data[:, 2],
        'err_M': data[:, 3],
        'E': data[:, 4],
        'err_E': data[:, 5],
        'C': data[:, 6],
        'err_C': data[:, 7]
    }

d_nn = load_data('termodinamica.dat')
d_nnn = load_data('termodinamica-next-next.dat')


# ==============================================================================
# 3. THERMODYNAMIC PLOTS (SEPARATED PANELS)
# ==============================================================================
def plot_thermodynamics(data, title, filename,
                       include_onsager=False,
                       Tc_sim=None,
                       xlim=(1.5, 6.0)):
    
    # Extraemos el nombre base sin el .png para añadir sufijos
    base_name = filename.replace('.png', '')

    # --------------------------------------------------------------------------
    # FIGURA 1: MAGNETIZACIÓN
    # --------------------------------------------------------------------------
    plt.figure(figsize=(6, 4))
    for n in np.unique(data['N']):
        mask = (data['N'] == n)
        plt.errorbar(data['T'][mask], data['M'][mask],
                     yerr=data['err_M'][mask],
                     fmt='o-', ms=4, capsize=2,
                     label=rf'$N={int(n)}$')

    if include_onsager:
        plt.plot(T_teorica, M_teorica, 'k--', lw=1.8,
                 label=r'Onsager ($N \to \infty$)')
        plt.axvline(Tc_onsager, color='black', linestyle=':', linewidth=1)

    if Tc_sim is not None:
        plt.axvline(Tc_sim, color='gray', linestyle='--', linewidth=1.2,
                    label=rf'$T_c^{{\mathrm{{sim}}}} = {Tc_sim:.3f}$')

    plt.ylabel(r'$|m|$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nMagnetization")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{base_name}_mag.png", dpi=300, bbox_inches='tight')
    plt.show()

    # --------------------------------------------------------------------------
    # FIGURA 2: ENERGÍA
    # --------------------------------------------------------------------------
    plt.figure(figsize=(6, 4))
    for n in np.unique(data['N']):
        mask = (data['N'] == n)
        plt.errorbar(data['T'][mask], data['E'][mask],
                     yerr=data['err_E'][mask],
                     fmt='s-', ms=4, capsize=2,
                     label=rf'$N={int(n)}$')

    if include_onsager:
        plt.axvline(Tc_onsager, color='black', linestyle=':', linewidth=1)

    if Tc_sim is not None:
        plt.axvline(Tc_sim, color='gray', linestyle='--', linewidth=1.2)

    plt.ylabel(r'$e$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nEnergy")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{base_name}_ene.png", dpi=300, bbox_inches='tight')
    plt.show()

    # --------------------------------------------------------------------------
    # FIGURA 3: CALOR ESPECÍFICO
    # --------------------------------------------------------------------------
    plt.figure(figsize=(6, 4))
    for n in np.unique(data['N']):
        mask = (data['N'] == n)
        plt.errorbar(data['T'][mask], data['C'][mask],
                     yerr=data['err_C'][mask],
                     fmt='^-', ms=4, capsize=2,
                     label=rf'$N={int(n)}$')

    if include_onsager:
        plt.axvline(Tc_onsager, color='black', linestyle=':', linewidth=1)

    if Tc_sim is not None:
        plt.axvline(Tc_sim, color='gray', linestyle='--', linewidth=1.2)

    plt.ylabel(r'$c_V$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nSpecific Heat")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(f"{base_name}_cv.png", dpi=300, bbox_inches='tight')
    plt.show()


# ==============================================================================
# 4. FINITE-SIZE SCALING
# ==============================================================================
def extract_tc(data):
    Tc_list = []
    invN = []

    for n in np.unique(data['N']):
        mask = (data['N'] == n)
        T = data['T'][mask]
        C = data['C'][mask]

        Tc_list.append(T[np.argmax(C)])
        invN.append(1.0 / n)

    # Use linregress to get standard error of the intercept
    res = linregress(invN, Tc_list)
    return invN, Tc_list, res.slope, res.intercept, res.intercept_stderr

invN_nn, Tc_nn, slope_nn, Tc_inf_nn, err_Tc_inf_nn = extract_tc(d_nn)
invN_nnn, Tc_nnn, slope_nnn, Tc_inf_nnn, err_Tc_inf_nnn = extract_tc(d_nnn)


# ==============================================================================
# 5. FIGURES
# ==============================================================================

# NN → ONLY ONSAGER
plot_thermodynamics(
    d_nn,
    r'2D Ising Model (Nearest Neighbors)',
    'fig_nn.png',
    include_onsager=True,
    Tc_sim=None,
    xlim=(1.4, 3.6)
)

# NNN → ONLY SIMULATED Tc
plot_thermodynamics(
    d_nnn,
    r'2D Ising Model (Nearest + Next-Nearest Neighbors)',
    'fig_nnn.png',
    include_onsager=False,
    Tc_sim=Tc_inf_nnn,
    xlim=(5.8, 11.1)
)


# ==============================================================================
# 6. Tc EXTRAPOLATION
# ==============================================================================
plt.figure(figsize=(6, 4.5))

plt.plot(invN_nn, Tc_nn, 'o', label='NN data')
plt.plot(invN_nn, slope_nn*np.array(invN_nn) + Tc_inf_nn,
         '-', label=rf'NN fit ($T_c={Tc_inf_nn:.3f} \pm {err_Tc_inf_nn:.3f}$)')

plt.plot(invN_nnn, Tc_nnn, 's', label='NN+NNN data')
plt.plot(invN_nnn, slope_nnn*np.array(invN_nnn) + Tc_inf_nnn,
         '-', label=rf'NN+NNN fit ($T_c={Tc_inf_nnn:.3f} \pm {err_Tc_inf_nnn:.3f}$)')

plt.plot([0], [Tc_onsager], 'k*', ms=10,
         label=r'Onsager ($T_c=2.269$)')

plt.xlabel(r'$1/N$')
plt.ylabel(r'$T_c(N)$')
plt.title(r'Finite-size scaling extrapolation')
plt.legend(frameon=False)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('fig_tc_extrapolation.png', dpi=300)
plt.show()


# ==============================================================================
# 7. RESULTS
# ==============================================================================
print("\n" + "=" * 60)
print("CRITICAL TEMPERATURE ESTIMATES")
print("=" * 60)
print(f"Onsager exact (NN):           Tc = {Tc_onsager:.4f}")
print(f"Simulation NN (extrapolated): Tc = {Tc_inf_nn:.4f} +/- {err_Tc_inf_nn:.4f}")
print(f"Simulation NN+NNN (extra.):   Tc = {Tc_inf_nnn:.4f} +/- {err_Tc_inf_nnn:.4f}")
print("=" * 60)