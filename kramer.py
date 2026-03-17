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

def coef_correlacion_lineal(x,y,popt):
    scr = np.sum((y - f(x, *popt))**2)
    sct = np.sum((y - np.mean(y))**2)
    return 1 - scr/sct

def analisis():
    data_df = pd.read_excel("DATOS/kramer_datos.ods", engine="odf")

    lnV = np.log(data_df["V_kV"].to_numpy() * 1E3)
    lnQ = np.log(data_df["Q"].to_numpy()) 

    
    popt, pcov = curve_fit(f,lnV, lnQ)

    m, n         = popt
    m_err, n_err = np.sqrt(np.diag(pcov))

    resultados_df = pd.DataFrame()
    resultados_df.at[0,"m"], resultados_df.at[0,"dm"] = ut.redondear_escalares(m, m_err)
    resultados_df.at[0,"n"], resultados_df.at[0,"dn"] = ut.redondear_escalares(n, n_err)
    resultados_df.at[0,"R2"] = coef_correlacion_lineal(lnV,lnQ,popt)

    resultados_df.to_excel("RESULTADOS/kramer_resultados.ods", engine="odf")

    plt.figure(figsize=(8, 6))
    plt.scatter(lnV, lnQ, label='Datos experimentales', color='red')
    plt.plot(lnV, f(lnV, *popt), label='Ajuste lineal') 
    plt.xlabel(r'$\ln(V)$')
    plt.ylabel(r'$\ln(A)$')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    analisis()