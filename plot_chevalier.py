import numpy as np
import matplotlib.pyplot as plt
from array_io import *



fname = "chevalier.txt"
r, M, us, Ps, rhos = read_five_arrays(fname)

fname = "chevalier.dimensionfull.txt"
x, M, u, P, rho, T = read_six_arrays(fname)

mp = 1.6737236e-24
n = rho / mp

log_r    = np.log10(r)
log_rhos = np.log10(rhos)
log_us   = np.log10(us)
log_Ps   = np.log10(Ps)
log_x    = np.log10(x)
log_M    = np.log10(M)
log_u    = np.log10(u)
log_P    = np.log10(P)
log_n    = np.log10(n)
log_T    = np.log10(T)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2)

ax1.plot(log_r, log_rhos)
ax1.plot(log_r, log_us)
ax1.plot(log_r, log_Ps)
ax1.set_xlim([-0.5, 0.5])
ax1.set_ylim([-4.0, 1.0])
ax1.set_xlabel(r'$log_{10}(r/R_{\ast})$')
ax1.text(0.15,-1.3,r'$log_{10}(\rho_{\ast})$',color='blue')
ax1.text(0.15,0.3,r'$log_{10}(u_{\ast})$',color='green')
ax1.text(0.15,-3.5,r'$log_{10}(P_{\ast})$',color='red')

ax2.plot(x,log_n,'-',color='black')
ax2.set_xlim([0, 2000])
ax2.set_ylim([-3.0, 0.0])
ax2.set_xlabel('r [pc]')
ax2.set_ylabel('$log_{10}(n) [cm^{-3}]$')
ax2.vlines(1000, -3.0, 0.0, linestyle='dashed')

ax3.plot(x,log_u,'-',color='black')
ax3.set_xlim([0, 2000])
ax3.set_ylim([0, 4.0])
ax3.set_xlabel('r [pc]')
ax3.set_ylabel('$log_{10}(u) [km/s]$')
ax3.vlines(1000, 0.0, 4.0, linestyle='dashed')

ax4.plot(x,log_T,'-',color='black')
ax4.set_xlim([0, 2000])
ax4.set_ylim([6.0, 8.0])
ax4.set_xlabel('r [pc]')
ax4.set_ylabel('$log_{10}(T) [K]$')
ax4.vlines(1000, 6.0, 8.0, linestyle='dashed')

fig.subplots_adjust(bottom=0.1)
fig.subplots_adjust(top=0.95)
fig.subplots_adjust(wspace=0.3)
fig.subplots_adjust(hspace=0.3)


s = 'chevalier.png'
plt.savefig(s,bbox_inches='tight')
plt.show()



