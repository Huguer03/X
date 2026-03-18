import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.constants import e, c
import utiles as ut
import scienceplots
plt.style.use(['science'])    

def f(x, m, n):
    """
    Función lineal para el ajuste: y = m*x + n
    """
    return m * x + n

def coef_correlacion_lineal(y,x,popt):
    scr = np.sum((y - f(x, *popt))**2)
    sct = np.sum((y - np.mean(y))**2)
    return 1 - scr/sct

def analisis():
    data_df = pd.read_excel("DATOS/planck_datos.ods", engine="odf")

    V_inbersa  = 1 / (data_df["V_kV"].to_numpy() * 1E3)
    lambda_min = data_df["lambda_A"].to_numpy() * 1E-10

    
    popt, pcov = curve_fit(f,V_inbersa, lambda_min)

    m, n         = popt
    m_err, n_err = np.sqrt(np.diag(pcov))

    h     = m * e / c
    h_err = m_err * e / c

    resultados_df = pd.DataFrame()
    resultados_df.at[0,"m"], resultados_df.at[0,"dm"] = ut.redondear_escalares(m, m_err)
    resultados_df.at[0,"n"], resultados_df.at[0,"dn"] = ut.redondear_escalares(n, n_err)
    resultados_df.at[0,"R2"] = coef_correlacion_lineal(lambda_min,V_inbersa,popt)
    resultados_df.at[0,"h"], resultados_df.at[0,"h_err"] = ut.redondear_escalares(h, h_err)

    resultados_df.to_excel("RESULTADOS/planck_resultados.ods", engine="odf")

    plt.figure(figsize=(8, 6))
    plt.scatter(V_inbersa, lambda_min, label='Datos experimentales', color='red')
    plt.plot(V_inbersa, f(V_inbersa, *popt), label='Ajuste lineal') 
    plt.xlabel(r'$\frac{1}{V}$ (V)')
    plt.ylabel(r'$\lambda$ (m)')
    plt.legend()
    plt.grid(True)
    plt.savefig("FIGURAS/planck.png")
    plt.show()

if __name__ == "__main__":
    analisis()