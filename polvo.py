import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter, find_peaks
from pybaselines import Baseline
import utiles as ut
import scienceplots
plt.style.use(['science'])  

def analisis():
    data_df = pd.read_excel("DATOS/polvo.xlsx")

    beta = data_df["Beta"].to_numpy() + 0.4
    R_0  = data_df["R_0"].to_numpy()

    baseline_fitter = Baseline(x_data=beta)
    
    fondo, _ = baseline_fitter.asls(R_0, lam=1e5, p=0.2)
    
    R_neto = R_0 - fondo

    R_suave = savgol_filter(R_neto, 11, 4)

    indices_picos, _ = find_peaks(R_suave[beta>=13], height=0.5, prominence=0.2, distance=10, width=3)

    beta_picos, R_picos = beta[beta>=13][indices_picos], R_suave[beta>=13][indices_picos]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10), sharex=True)

    ax1.plot(beta, R_0, label='Original', alpha=0.6)
    ax1.plot(beta, fondo, 'r--', label='Fondo (AsLS)')
    ax1.set_ylabel(r'$R(\frac{1}{s})$')
    ax1.legend()

    ax2.plot(beta, R_neto, alpha=0.3, label='Señal con ruido')
    ax2.plot(beta, R_suave, color='black', linewidth=1.5, label='Señal Suavizada')
    for i,r in zip(beta_picos, R_picos):
        ax2.axvline(i, color="red", linestyle="--")
        ax2.text(i + 0.2, r, fr'${i:.1f}º', color='black', fontsize=9, fontweight='bold', verticalalignment='bottom')
    
    ax2.set_xlabel(r'$\beta$ (º)')
    ax2.set_ylabel("Intensidad")
    ax2.legend()
    ax2.grid(True)

    plt.tight_layout()
    plt.savefig("FIGURAS/polvo.png")
    plt.show()

if __name__ == "__main__":
    analisis()