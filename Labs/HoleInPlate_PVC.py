##### STRESS AND STRAIN FOR A HOLE IN A PLATE #####
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.ion()
plt.close('all')
mpl.rcParams.update({'font.size': 15})

## Define sample properties
E = 3.2e9 #Young's Modulus [N/m^2]
nu = 0.4  #Poisson's Ratio []
G = E/(2*(1+nu))

## Define Sample Dimensions
W = .03 	#ROI Width [m]
H = .12     #ROI Height [m]
T = 0.00313 #Thickness[m]
a = 0.003125 #Hole radius [m]

## Set Applied Load
Load = 5000 #Applied load [N]
so = Load/(W*T)	#Applied stress [N/m^2]

## Set Gauge Positions
G1 = [0, 0.043]  #Gauge 1 (x,y) [m]
G2 = [0.0054, 0]  #Gauge 2 (x,y) [m]
G3 = [0, 0.0051]  #Gauge 3 (x,y) [m]
G4 = [0.01, 0]   #Gauge 4 (x,y) [m]

## Create Meshgrid
n = 1000 	#Spatial Resolution
x = np.linspace(-W/2,W/2,n)	#x grid
y = np.linspace(-H/2,H/2,n)	#y grid
[X,Y] = np.meshgrid(x,y)
r = (X**2 + Y**2)**0.5		#r grid
th = np.arctan(Y/X) + np.pi/2		#theta grid
twoPi = np.linspace(0,2*np.pi,100)
xr = a*np.cos(twoPi)	#Hole x
yr = a*np.sin(twoPi)	#Hole y

#Convert Gauge Positions to Indices
P = [[int(n/2 + x[0]/W*n), int(n/2 + x[1]/H*n)] for x in [G1,G2,G3,G4]]

#%% Calculate Stress and Strain
## Calculate Hoop and Radial Stresses
srr = so/2*((1-a**2/r**2) + (1+3*a**4/r**4-4*a**2/r**2)*np.cos(2*th))
srr[r<=a] = 0
stt = so/2*((1+a**2/r**2) - (1+3*a**4/r**4)*np.cos(2*th))
stt[r<=a] = 0
srt =-so/2*(1 + 2*a**2/r**2 -3*a**4/r**4)*np.sin(2*th)
srt[r<=a] = 0

## Convert to Cartesian Stress
sxx = srr*np.cos(th)**2 + stt*np.sin(th)**2 - srt*np.sin(2*th)
syy = srr*np.sin(th)**2 + stt*np.cos(th)**2 + srt*np.sin(2*th)
sxy = np.sin(th)*np.cos(th)*(srr-stt) + srt*np.cos(2*th)

## Convert to Strain
exx = (sxx - nu*syy)/E
eyy = (syy - nu*sxx)/E
exy = sxy/G

#%% Plot Strain Response at Each Gauge
fig,axs = plt.subplots(1,2,figsize=(8,6))
for i,p in enumerate(P):
    axs[0].plot([0, exx[p[1],p[0]]*100],[0,Load],label = 'Gauge '+str(i+1))
    axs[1].plot([0, eyy[p[1],p[0]]*100],[0,Load],label = 'Gauge '+str(i+1))
axs[0].set_title('Axial')
axs[0].set_xlabel('Strain [%]')
axs[0].set_ylabel('Load (N)')
axs[1].set_title('Transverse')
axs[1].set_xlabel('Strain [%]')
axs[1].set_ylabel('Load [N]')
plt.legend();
plt.tight_layout()

#%% Plot Strain Field
Strain = [exx, eyy, exy]
eType = ['{xx}','{yy}','{xy}']
tOff = 0.001
fig, axs = plt.subplots(1,3,figsize=(10,8))
for i,ax in enumerate(axs):
    eC = ax.contourf(X,Y,Strain[i]*100, vmin=-so/E*100, vmax=3*so/E*100, cmap='jet', levels=100)
    ax.plot(xr,yr,'k--')
    ax.set_title('$\epsilon_'+eType[i]+'$')
    ax.set_aspect('equal', 'box')
    for j,g in enumerate([G1,G2,G3,G4]):
        ax.plot(g[0],g[1],'kx',ms=10,mew=3)
        ax.text(g[0]+tOff,g[1]+tOff,str(j+1),fontsize=15)
fig.colorbar(eC, label='Strain (%)')
plt.tight_layout()

#%% Plot Stress Field
Stress = [srr,stt,srt]
sType = ['{rr}','{θθ}','{rθ}']

fig, axs = plt.subplots(1,3,figsize=(10,8))
for i,ax in enumerate(axs):
    sC = ax.contourf(X,Y,Stress[i]*1e-6, vmin=-so*1e-6, vmax=3*so*1e-6, cmap='jet', levels=20)
    ax.plot(xr,yr,'k--')
    ax.set_title('$\sigma_'+sType[i]+'$')
    ax.set_aspect('equal', 'box')
fig.colorbar(sC, label='Stress (MPa)')
plt.tight_layout()