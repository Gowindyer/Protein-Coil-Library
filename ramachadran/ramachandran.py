import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
import sys

def getProbabilityDistribution(grid, phi, psi, sigma):
    psi_, phi_ = grid[:,0], grids[:, 1]
    number = np.zeros_like(psi_)
    for i in range(len(phi)):
        for j in range(len(phi_)): 
            delta_psi = min(abs(psi[i] - psi_[j]), 360.0 - abs(psi[i] - psi_[j]))
            delta_phi = min(abs(phi[i] - phi_[j]), 360.0 - abs(phi[i] - phi_[j]))
            number[j] += np.exp(-0.5*(delta_psi**2 + delta_phi**2)/sigma**2)
    Probability = (number + 0.2)/np.sum(number)
    return Probability

def plot(fig, X, Y, Z, num_levels=8):
    plt.xlabel(r'$\phi (^\circ)$',fontsize=24)
    plt.ylabel(r'$\psi (^\circ)$',fontsize=24)
    levels = np.linspace(2.7, 12, num_levels)
    fc = plt.contourf(X, Y, Z, cmap=cm.gist_gray, levels=levels)
    cbar = plt.colorbar(fc,fraction=0.046, pad=0.04)
    cbar.ax.get_yaxis().labelpad = 20
    #cbar.set_label(label=colorBarLabel, rotation=270,size='18')
    cbar.ax.tick_params(labelsize=12)
    return fig


residuename = str(sys.argv[1])
resPhiPsi = []
PhiPsi = np.load('../proteinCoilLibrary/AllResiduePhiPsi.npy',allow_pickle='TRUE').item()
resPhiPsi = PhiPsi[residuename]
phiedges = np.linspace(-180,180,20)
psiedges = np.linspace(-180,180,20)
resPhiPsi = np.array(resPhiPsi)*180/np.pi
phi = resPhiPsi[:,0]
psi = resPhiPsi[:,1]

num_grids = 36
num_levels = 9
sigma = 10
psi_grids = np.linspace(-180, 180, num_grids)
phi_grids = np.linspace(-180, 180, num_grids)
psi_grids = (psi_grids[:-1] + psi_grids[1:])/2
phi_grids = (phi_grids[:-1] + phi_grids[1:])/2
grids = np.array(np.meshgrid(psi_grids, phi_grids)).T.reshape((-1, 2))

Probability = getProbabilityDistribution(grids, phi, psi, sigma)
Z = Probability.reshape(psi_grids.shape[0], phi_grids.shape[0])
X, Y = np.meshgrid(phi_grids, psi_grids)
fig, ax = plt.subplots(figsize=(4.6,4))
#plt.title('Ramachandran plot for {}'.format(residuename),fontsize=20)
plot(fig, X, Y, -np.log(Z), num_levels=num_levels)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.set_xticks(np.arange(-160,200,40))
ax.set_yticks(np.arange(-160,200,40))
ax.tick_params(axis='x', which='minor')
ax.tick_params(axis='y', which='minor')
ax.xaxis.set_minor_locator(MultipleLocator(20))
ax.yaxis.set_minor_locator(MultipleLocator(20))
#ax.minorticks_on()
ax.tick_params(which='both' , width=1.0)
#plt.plot(phi, psi, 'o', color='black',alpha=0.1)
plt.savefig('ramachandran_%s.pdf'%residuename,dpi=600,bbox_inches='tight')

