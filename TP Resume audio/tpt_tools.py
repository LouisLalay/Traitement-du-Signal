import matplotlib.pyplot as plt
import numpy as np

def F_plot1(x_v, y_v, labelX, labelY):
    plt.plot(x_v, y_v)
    plt.xlabel(labelX)
    plt.ylabel(labelY)
    plt.grid(True)
    return

def F_plot2(data_m, col_v=np.zeros(0), row_v=np.zeros(0), labelCol='', labelRow=''):
    plt.imshow(data_m, origin='lower', aspect='auto', extent=[row_v[0], row_v[-1], col_v[0], col_v[-1]], interpolation='nearest')
    plt.colorbar()
    plt.set_cmap('gray_r')
    plt.xlabel(labelRow)
    plt.ylabel(labelCol)
    plt.grid(True)
    return

def nextpow2(i):
    n = 1
    while n < i:
        n *= 2
    return n
