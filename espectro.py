import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
import utiles as ut
import scienceplots
plt.style.use(['science'])  

def analisis():
    d = 2.8201
    data_df = pd.read_excel("DATOS/espectro.xlsx")

    beta = data_df["Beta"].to_numpy() 
    R_0  = data_df["R_0"].to_numpy()
    R_1  = data_df["R_1"].to_numpy()
    R_2  = data_df["R_2"].to_numpy()
    R_3  = data_df["R_3"].to_numpy()

    Q        = np.asarray([np.max(R_0[beta<=5]), np.max(R_1[beta<=5]), np.max(R_2[beta<=5.7]), np.max(R_3[beta<=6.6])])
    beta_min = np.asarray([beta[R_0==Q[0]][0],beta[R_1==Q[1]][0],beta[R_2==Q[2]][0],beta[R_3==Q[3]][0]])

    resultados_df = pd.DataFrame()
    resultados_df["V_kV"] = [35, 30, 25, 20]
    resultados_df["Q"] = Q
    resultados_df.to_excel("DATOS/kramer_datos.ods", engine="odf")

    lambda_min = 2 * d * np.sin(beta_min * np.pi / 180)
 
    resultados_df = pd.DataFrame()
    resultados_df["V_kV"] = [35, 30, 25, 20]
    resultados_df["lambda_A"] = lambda_min
    resultados_df.to_excel("DATOS/planck_datos.ods", engine="odf")

    indices_picos, _ = find_peaks(R_0, height=200)

    K   = beta[indices_picos]
    R_k = R_0[indices_picos]

    plt.figure(figsize=(10, 8))
    plt.scatter(beta_min, Q, color="purple", s=25, marker="x", label=r"$\beta_{min}$")
    plt.plot(beta, R_0, color="lightseagreen", label = r"$V=35kV$")
    plt.plot(beta, R_1, color="green", label = r"$V=30kV$") 
    plt.plot(beta, R_2, color="orange", label = r"$V=25kV$") 
    plt.plot(beta, R_3, color="red", label = r"$V=20kV$")
    plt.axvline(K[0], color="black", ls="--", alpha=0.3)
    plt.text(K[0] - 0.6, R_k[0], r'$K_{\beta}$', color='black', fontsize=9, fontweight='bold', verticalalignment='bottom')
    plt.axvline(K[1], color="black", ls="--", alpha=0.3)
    plt.text(K[1] + 0.2, R_k[1], r'$K_{\alpha}$', color='black', fontsize=9, fontweight='bold', verticalalignment='bottom')
    plt.xlabel(r'$\beta$ (º)')
    plt.ylabel(r'$R(\frac{1}{s})$')
    plt.legend()
    plt.grid(True)
    plt.savefig("FIGURAS/espectro.png", dpi=300)
    plt.show()

if __name__ == "__main__":
    analisis()