import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.constants import e, c
import utiles as ut
import scienceplots
plt.style.use(['science'])  

def analisis():
    delta_beta = 0.2

    data_df = pd.read_excel("DATOS/espectro.xlsx")

    beta = data_df["Beta"].to_numpy()
    R_0  = data_df["R_0"].to_numpy()
    R_1  = data_df["R_1"].to_numpy()
    R_2  = data_df["R_2"].to_numpy()
    R_3  = data_df["R_3"].to_numpy()

    resultados_df = pd.DataFrame()
    resultados_df["V_kV"] = [35, 30, 25, 20]

    Q = np.zeros(4)
    for i in zip(R_0,R_1,R_2,R_3):
        Q += np.asarray(i) * delta_beta
    
    resultados_df["Q"] = Q

    resultados_df.to_excel("DATOS/kramer_datos.ods", engine="odf")

    plt.figure(figsize=(8, 6))
    plt.plot(beta, R_0, label = r"$V=35kV$")
    plt.plot(beta, R_1, label = r"$V=30kV$") 
    plt.plot(beta, R_2, label = r"$V=25kV$") 
    plt.plot(beta, R_3, label = r"$V=20kV$")  
    plt.xlabel(r'$\beta$ (º)')
    plt.ylabel(r'$R(\frac{1}{s})$')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    analisis()