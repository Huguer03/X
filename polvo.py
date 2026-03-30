from os import error
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks
from scipy.constants import N_A
from pybaselines import Baseline
import utiles as ut
import scienceplots
plt.style.use(['science'])  

def plot_difractograma(beta:np.ndarray, R_0:np.ndarray, R_suave:np.ndarray, R_neto:np.ndarray, beta_picos:np.ndarray, R_picos:np.ndarray, fondo:np.ndarray) -> None:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    ax1.plot(beta, R_0, label='Original', alpha=0.6)
    ax1.plot(beta, fondo, 'r--', label='Fondo (AsLS)')
    ax1.set_ylabel(r'$R(\frac{1}{s})$')
    ax1.legend()

    ax2.plot(beta, R_neto, alpha=0.3, label='Señal con ruido')
    ax2.plot(beta, R_suave, color='black', linewidth=1.5, label='Señal Suavizada')
    for i,r in zip(beta_picos, R_picos):
        ax2.axvline(i, color="red", linestyle="--")
        ax2.text(i + 0.2, r, fr'${i:.1f}$º', color='black', fontsize=9, fontweight='bold', verticalalignment='bottom')
    
    ax2.set_xlabel(r'$\beta$ (º)')
    ax2.set_ylabel("Intensidad")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("FIGURAS/polvo.png", dpi=300)

def plot_indexado(beta:np.ndarray, beta_picos:np.ndarray, R_picos:np.ndarray, R_suave:np.ndarray, miller_index:list[str], path:str="imagen.png") -> None:
    plt.figure(figsize=(10,8))
    plt.plot(beta, R_suave, color='black', linewidth=1.5)
    for i,r,mi in zip(beta_picos, R_picos, miller_index):
        plt.axvline(i, color="red", linestyle="--")
        plt.text(i + 0.2, r, mi, color='black', fontsize=12, fontweight='bold', verticalalignment='bottom')
    
    plt.xlabel(r'$\beta$ (º)')
    plt.ylabel("Intensidad")

    plt.tight_layout()
    plt.savefig(path, dpi=300)

def indexaccion(beta_picos:np.ndarray, longituz_onda:float=1.5406, estructura:str="fcc") -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    theta = np.radians(beta_picos)
    sin2  = np.sin(theta)**2

    if estructura.lower() == "fcc":
        f_norm = sin2[0]/3
    elif estructura.lower() == "bcc":
        f_norm = sin2[0]/2
    else: 
        f_norm = sin2[0]

    s_raw = sin2 / f_norm
    s = np.round(s_raw).astype(int)
    a = (longituz_onda / 2) * (np.sqrt(s) / np.sin(theta))

    miller_index = ["" for _ in beta_picos]

    for i in range(len(s)):
        limite = int(np.sqrt(s[i])) + 1
        for h in range(limite, -1, -1):
            for k in range(h, -1, -1):
                for l in range(k, -1, -1):
                    if h**2 + k**2 + l**2 == s[i]:
                        miller_index[i] = f"({h}{k}{l})"

    return sin2, s, a, miller_index

def analisis():
    rho = 2.165 
    P_m = 58.5

    data_df = pd.read_excel("DATOS/polvo.xlsx")

    beta = data_df["Beta"].to_numpy() + 0.4
    R_0  = data_df["R_0"].to_numpy()

    baseline_fitter = Baseline(x_data=beta)
    
    fondo, _ = baseline_fitter.asls(R_0, lam=1e5, p=0.2)
    
    R_neto = R_0 - fondo

    R_suave = savgol_filter(R_neto, 11, 4)

    indices_picos, _ = find_peaks(R_suave[beta>=13], height=0.5, prominence=0.2, distance=10, width=3)

    beta_picos, R_picos = beta[beta>=13][indices_picos], R_suave[beta>=13][indices_picos]

    plot_difractograma(beta, R_0, R_suave, R_neto, beta_picos, R_picos, fondo)

    sin2, s, a, mi = indexaccion(beta_picos)

    a_m       = np.sum(a)/len(a)
    error_str = np.sqrt(np.sum((a-a_m)**2) / ((len(a)-1) * len(a)))

    Z     = N_A * rho * (a_m*1e-8)**3 / P_m
    Z_err = 3.0 * Z * error_str / a_m

    resultados_df = pd.DataFrame()
    resultados_df[r"$\theta$"]            = np.around(np.radians(beta_picos),3)
    resultados_df[r"$\sin^2(\theta)$"]    = np.around(sin2,3)
    resultados_df[r"$h^2 + k^2 + l^2$"]   = s
    resultados_df[r"$(hkl)$"]             = mi
    resultados_df[r"$a(\mathring{\text{A}})$"] = np.around(a,3)
    resultados_df.at[0,r"$\bar{a}(\mathring{\text{A}})$"], resultados_df.at[0,r"$\Delta a(\dot{\text{A}})$"] = ut.redondear_escalares(a_m, error_str)
    resultados_df.at[0, r"$Z$"], resultados_df.at[0,r"$\Delta Z$"] = ut.redondear_escalares(Z,Z_err)
    resultados_df.to_excel("RESULTADOS/polvo.ods", engine="odf")

    plot_indexado(beta, beta_picos, R_picos, R_suave, mi, "FIGURAS/indexado.png")

    sin2, s, a, mi = indexaccion(beta_picos[1:], estructura="sc")

    a_m       = np.sum(a)/len(a)
    error_str = np.sqrt(np.sum((a-a_m)**2) / ((len(a)-1) * len(a)))

    Z     = N_A * rho * (a_m*1e-8)**3 / P_m
    Z_err = 3.0 * Z * error_str / a_m

    resultados_df = pd.DataFrame()
    resultados_df[r"$\theta$"]            = np.around(np.radians(beta_picos[1:]),3)
    resultados_df[r"$\sin^2(\theta)$"]    = np.around(sin2,3)
    resultados_df[r"$h^2 + k^2 + l^2$"]   = s
    resultados_df[r"$(hkl)$"]             = mi
    resultados_df[r"$a(\mathring{\text{A}})$"] = np.around(a,3)
    resultados_df.at[0,r"$\bar{a}(\mathring{\text{A}})$"], resultados_df.at[0,r"$\Delta a(\dot{\text{A}})$"] = ut.redondear_escalares(a_m, error_str)
    resultados_df.at[0, r"$Z$"], resultados_df.at[0,r"$\Delta Z$"] = ut.redondear_escalares(Z,Z_err)
    resultados_df.to_excel("RESULTADOS/polvo_defectuoso.ods", engine="odf")

    plot_indexado(beta, beta_picos[1:], R_picos[1:], R_suave, mi, "FIGURAS/indexado_defectuoso.png")

if __name__ == "__main__":
    analisis()