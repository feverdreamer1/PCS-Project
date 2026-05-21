import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.optimize import curve_fit
plt.rcParams.update({
    "font.size": 12,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "figure.dpi": 120
})

Tc_wannier = 4.0 / np.log(3.0) 
beta_teorico = 0.125
nu_teorico = 1.000

def load_thermo_data(file):
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

def load_corr_data(file):
    data = np.loadtxt(file)
    return data[:, 0], data[:, 1], data[:, 2]
def plot_thermodynamics(data, title, base_name,
                        Tc_sim=None,
                        Tc_exact=None,
                        xlim=(2.0, 6.0)):
    plt.figure(figsize=(6, 4))

    for n in np.unique(data['N']):

        mask = (data['N'] == n)

        plt.errorbar(
            data['T'][mask],
            data['M'][mask],
            yerr=data['err_M'][mask],
            fmt='o-',
            ms=4,
            capsize=2,
            label=rf'$N={int(n)}$'
        )

    if Tc_exact is not None:
        plt.axvline(
            Tc_exact,
            color='black',
            linestyle=':',
            linewidth=1.5,
            label=r'Wannier exact'
        )

    if Tc_sim is not None:
        plt.axvline(
            Tc_sim,
            color='gray',
            linestyle='--',
            linewidth=1.2,
            label=rf'$T_c^{{sim}}={Tc_sim:.3f}$'
        )

    plt.ylabel(r'Magnetization $|m|$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nMagnetization")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(f"{base_name}_mag.png", dpi=300)
    plt.show()

    plt.figure(figsize=(6, 4))

    for n in np.unique(data['N']):

        mask = (data['N'] == n)

        plt.errorbar(
            data['T'][mask],
            data['E'][mask],
            yerr=data['err_E'][mask],
            fmt='s-',
            ms=4,
            capsize=2,
            label=rf'$N={int(n)}$'
        )

    if Tc_exact is not None:
        plt.axvline(Tc_exact, color='black', linestyle=':', linewidth=1.5)

    if Tc_sim is not None:
        plt.axvline(Tc_sim, color='gray', linestyle='--', linewidth=1.2)

    plt.ylabel(r'Energy $e$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nEnergy")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(f"{base_name}_ene.png", dpi=300)
    plt.show()

    # --------------------------------------------------------------------------
    # SPECIFIC HEAT
    # --------------------------------------------------------------------------
    plt.figure(figsize=(6, 4))

    for n in np.unique(data['N']):

        mask = (data['N'] == n)

        plt.errorbar(
            data['T'][mask],
            data['C'][mask],
            yerr=data['err_C'][mask],
            fmt='^-',
            ms=4,
            capsize=2,
            label=rf'$N={int(n)}$'
        )

    if Tc_exact is not None:
        plt.axvline(Tc_exact, color='black', linestyle=':', linewidth=1.5)

    if Tc_sim is not None:
        plt.axvline(Tc_sim, color='gray', linestyle='--', linewidth=1.2)

    plt.ylabel(r'Specific Heat $c_V$')
    plt.xlabel(r'$k_B T / J$')
    plt.title(f"{title}\nSpecific Heat")
    plt.grid(alpha=0.3)
    plt.xlim(xlim)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(f"{base_name}_cv.png", dpi=300)
    plt.show()

def extract_tc(data):

    Tc_list = []
    invN = []

    for n in np.unique(data['N']):

        mask = (data['N'] == n)

        T = data['T'][mask]
        C = data['C'][mask]

        Tc_list.append(T[np.argmax(C)])

        # Since ν=1 assumed
        invN.append(n**(-1.0))

    res = linregress(invN, Tc_list)

    return (
        invN,
        Tc_list,
        res.slope,
        res.intercept,
        res.intercept_stderr
    )

def calculate_beta(data, Tc, T_min, title, filename):

    mask_N = (data['N'] == 512)

    T = data['T'][mask_N]
    M = data['M'][mask_N]

    # Critical region below Tc
    mask_fit = (T >= T_min) & (T < Tc - 0.005)

    T_fit = T[mask_fit]
    M_fit = M[mask_fit]

    x = np.log(Tc - T_fit)
    y = np.log(M_fit)

    slope, intercept, _, _, _ = linregress(x, y)

    beta = slope

    plt.figure(figsize=(6, 4.5))

    plt.plot(
        x,
        y,
        'o',
        color='red',
        label='Simulation (N=512)'
    )

    plt.plot(
        x,
        slope*x + intercept,
        'r-',
        label=rf'Fit: $\beta={beta:.3f}$'
    )

    plt.title(f"{title}\nExponent $\\beta$")
    plt.xlabel(r'$\ln(T_c-T)$')
    plt.ylabel(r'$\ln(m)$')

    plt.grid(alpha=0.3)
    plt.legend(frameon=False)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()

    return beta

def corr_model(r, A, xi, eta):
    return A * np.exp(-r/xi) / (r**eta)

def calculate_nu(T_all,
                 r_all,
                 f_all,
                 Tc,
                 T_fit_list,
                 r_min,
                 r_max,
                 title,
                 filename):

    fig, axs = plt.subplots(1, 2, figsize=(11, 4.5))

    xi_list = []
    T_valid = []

    colors = plt.cm.plasma(
        np.linspace(0, 0.85, len(T_fit_list))
    )

    for i, T_target in enumerate(T_fit_list):

        mask = (
            np.isclose(T_all, T_target, atol=1e-6)
            & (r_all >= r_min)
            & (r_all <= r_max)
        )

        r_fit = r_all[mask]
        f_fit = f_all[mask]

        # Remove noisy values
        valid = (f_fit > 1e-6)

        r_fit = r_fit[valid]
        f_fit = f_fit[valid]

        if len(r_fit) < 6:
            continue

        try:

            # Initial parameters
            p0 = [1.0, 10.0, 0.25]

            popt, _ = curve_fit(
                corr_model,
                r_fit,
                f_fit,
                p0=p0,
                maxfev=10000
            )

            A_fit, xi_fit, eta_fit = popt

            if xi_fit > 0:

                xi_list.append(xi_fit)
                T_valid.append(T_target)

                r_smooth = np.linspace(
                    r_fit.min(),
                    r_fit.max(),
                    300
                )

                axs[0].plot(
                    r_fit,
                    np.log(f_fit),
                    'o',
                    color=colors[i],
                    markersize=4
                )

                axs[0].plot(
                    r_smooth,
                    np.log(corr_model(r_smooth, *popt)),
                    '-',
                    color=colors[i],
                    label=rf'$T={T_target}$ ($\xi={xi_fit:.2f}$)'
                )

        except RuntimeError:
            print(f"Fit failed for T={T_target}")

    axs[0].set_title(
        f"{title}\nCorrelation Function Fits"
    )

    axs[0].set_xlabel(r"Distance $r$")
    axs[0].set_ylabel(r"$\ln(f(r))$")

    axs[0].grid(alpha=0.3)
    axs[0].legend(frameon=False, fontsize=8)

    T_valid = np.array(T_valid)
    xi_list = np.array(xi_list)

    x_nu = np.log(T_valid - Tc)
    y_nu = np.log(xi_list)

    slope_nu, intercept_nu, _, _, _ = linregress(
        x_nu,
        y_nu
    )

    nu = -slope_nu

    axs[1].plot(
        x_nu,
        y_nu,
        's',
        color='black',
        label='Simulation'
    )

    axs[1].plot(
        x_nu,
        slope_nu*x_nu + intercept_nu,
        'r-',
        label=rf'Fit: $\nu={nu:.3f}$'
    )

    axs[1].set_title(r"Critical Exponent $\nu$")
    axs[1].set_xlabel(r"$\ln(T-T_c)$")
    axs[1].set_ylabel(r"$\ln(\xi)$")

    axs[1].grid(alpha=0.3)
    axs[1].legend(frameon=False)

    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.show()

    return nu

d_tri = load_thermo_data(
    r'C:\Users\shara\Desktop\Fisica Sistemas Complejos\TrabajoIising\Triangulo\termodinamica_triangular.dat'
)

invN_tri, Tc_tri, slope_tri, Tc_inf_tri, err_Tc_inf_tri = extract_tc(d_tri)

plot_thermodynamics(
    d_tri,
    r'2D Ising Model (Triangular Lattice)',
    'fig_tri',
    Tc_sim=Tc_inf_tri,
    Tc_exact=Tc_wannier,
    xlim=(2.0, 6.0)
)

plt.figure(figsize=(6, 4.5))

plt.plot(
    invN_tri,
    Tc_tri,
    'go',
    label='Simulation'
)

plt.plot(
    invN_tri,
    slope_tri*np.array(invN_tri) + Tc_inf_tri,
    'g-',
    label=rf'Fit ($T_c={Tc_inf_tri:.4f}\pm{err_Tc_inf_tri:.4f}$)'
)

plt.plot(
    [0],
    [Tc_wannier],
    'k*',
    ms=10,
    label=rf'Wannier exact ($T_c={Tc_wannier:.4f}$)'
)

plt.xlabel(r'$N^{-1}$')
plt.ylabel(r'$T_c(N)$')

plt.title(r'Finite-Size Scaling')

plt.legend(frameon=False)
plt.grid(alpha=0.3)
plt.ylim(3.5,3.8)
plt.tight_layout()
plt.savefig('fig_tri_tc_extrapolation.png', dpi=300)
plt.show()

calc_beta_tri = calculate_beta(
    d_tri,
    Tc_inf_tri,
    T_min=3.40,
    title="Triangular Lattice",
    filename='fig_tri_exponent_beta.png'
)

try:

    T_corr_all, r_corr_all, f_corr_all = load_corr_data(
        r'C:\Users\shara\Desktop\Fisica Sistemas Complejos\TrabajoIising\Triangulo\correlacion_triangular.dat'
    )

    # Temperatures above Tc
    T_fit_tri = [
        3.66,
        3.68,
        3.70,
        3.75,
        3.80,
        3.90
    ]

    calc_nu_tri = calculate_nu(
        T_corr_all,
        r_corr_all,
        f_corr_all,
        Tc_inf_tri,
        T_fit_tri,
        r_min=4,
        r_max=40,
        title="Triangular Lattice",
        filename='fig_tri_exponent_nu.png'
    )

except OSError:

    print(" correlation file not found.")

    calc_nu_tri = None

print("\n" + "="*60)
print("CRITICAL TEMPERATURE ESTIMATE (TRIANGULAR)")
print("="*60)

print(f"Wannier exact:                  Tc = {Tc_wannier:.4f}")
print(f"Simulation extrapolated:        Tc = {Tc_inf_tri:.4f} +/- {err_Tc_inf_tri:.4f}")

relative_error = abs(
    Tc_inf_tri - Tc_wannier
) / Tc_wannier * 100

print(f"Relative error:                 {relative_error:.3f}%")

print("="*60)

print("\n" + "="*60)
print("CRITICAL EXPONENTS")
print("="*60)

print(f"BETA (Theory = {beta_teorico:.3f})")
print(f"  Simulated beta = {calc_beta_tri:.3f}")

print("-"*60)

if calc_nu_tri is not None:

    print(f"NU (Theory = {nu_teorico:.3f})")
    print(f"  Simulated nu = {calc_nu_tri:.3f}")

print("="*60)
