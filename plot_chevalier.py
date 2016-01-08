import numpy as np
import matplotlib.pyplot as plt
from array_io import *

fname = "chevalier.dimensionfull.txt"

x, M, u, P, rho, T = read_six_arrays(fname)

mp = 1.6737236e-24
n = rho / mp

xmin = 0.01
xmax = 2000
ymin = -15.0
ymax =  10.0

fsize = 5.0

bg_size = 0.15
tg_size = 0.02
lg_size = 0.15
rg_size = 0.05

wx = 1.0 - lg_size - rg_size
wy = 1.0 - bg_size - tg_size

plt.figure(figsize=(fsize,fsize*wx/wy))
a0  = plt.axes([lg_size,bg_size,wx,wy])
plt.xlim(xmin,xmax)
plt.ylim(ymin,ymax)
plt.xlabel('r [pc]')
plt.ylabel('$log_{10}\,y$')

log_x   = np.log10(x)
log_M   = np.log10(M)
log_u   = np.log10(u)
log_P   = np.log10(P)
log_n   = np.log10(n)
log_T   = np.log10(T)
plt.plot(x,log_M,'-',color="0.0")
plt.text(450,-0.5,'$log_{10}\,M$',color='0.0')
plt.plot(x,log_u,'-',color="red")
plt.text(450,2.5,'$log_{10}\,(u)$ [km/s]',color='red')
plt.plot(x,log_P,'-',color="blue")
plt.text(450,-9.5,'$log_{10}\,(P)$ [dyne/$cm^{2}$]',color='blue')
plt.plot(x,log_n,'-',color="green")
plt.text(450,-4,r'$log_{10}\,(n)$ [cm$^{-3}$]',color='green')
plt.plot(x,log_T,'-',color="purple")
plt.text(450,6,r'$log_{10}\,(T)$ [K]',color='purple')

s = 'chevalier.png'
plt.savefig(s,bbox_inches='tight')
plt.show()



