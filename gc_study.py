import glob
import sys
import re
import matplotlib.pyplot   as plt
import numpy               as np
from   numpy import linalg as la

#----------------------------------------------------------
#-                     ------------------------------------
#-Define tool functions------------------------------------
#-                     ------------------------------------
#----------------------------------------------------------
def getIntFromFilename(filename):
    """ Return the 1st integer from filename string
    filename : filename string
    """
    return int(re.search(r'\d+', filename).group(0))

#----------------------------------------------------------
#-                         --------------------------------
#-Read barycenter variables--------------------------------
#-                         --------------------------------
#----------------------------------------------------------
baryList = glob.glob("./barycenters_vrbls_*")
baryList.sort(key=getIntFromFilename)
nLvl = int(input("How many grid levels to take into account (>=3)? "))
if nLvl > len(baryList):
    nMax = len(baryList) - 1
elif nLvl >=3 and nLvl <= len(baryList):
    nMax = nLvl - 1
else:
    sys.exit("ERROR : number of grid levels must be >= 3 !")
print("Number of grid levels = "+str(nMax+1))

ratio = 2                    #-grid refinement ratio
Nx    = []                   #-list of normalized cell numbers
for i in range(nMax):        #-include all grid levels except the finest
    Nx.append(pow(ratio, i))    # = [(2^dim)^i]^(1/dim)
NxLim = [0.5, Nx[-1]*2]      #-xlim for log-log plots

rhoAll = []
pAll   = []
TAll   = []
VxAll  = []
VyAll  = []
for baryFile in baryList[:nMax+1]:
    with open(baryFile, "r") as f:
        rho = []
        p   = []
        T   = []
        Vx  = []
        Vy  = []
        for line in f:
            vrbls = [float(vrbl) for vrbl in line.split()]
            rho.append(vrbls[0])
            p  .append(vrbls[1])
            T  .append(vrbls[2])
            Vx .append(vrbls[3])
            Vy .append(vrbls[4])
        rhoAll.append(rho)
        pAll  .append(p)
        TAll  .append(T)
        VxAll .append(Vx)
        VyAll .append(Vy)

rhoAll = np.array(rhoAll)
pAll   = np.array(pAll)
TAll   = np.array(TAll)
VxAll  = np.array(VxAll)
VyAll  = np.array(VyAll)

#----------------------------------------------------------
#-                                                   ------
#-Read barycenter coordinates for error scatter plots------
#-                                                   ------
#----------------------------------------------------------
with open("barycenters", "r") as f:
    nps   = int(f.readline())
    xList = []
    yList = []
    for line in f:
        co = [float(val) for val in line.split()]
        xList.append(co[0])
        yList.append(co[1])

#----------------------------------------------------------
#-                            -----------------------------
#-Compute infinity & 2nd norms-----------------------------
#-                            -----------------------------
#----------------------------------------------------------
rhoInfnorm = []
rho2ndnorm = []
pInfnorm   = []
p2ndnorm   = []
TInfnorm   = []
T2ndnorm   = []
VxInfnorm  = []
Vx2ndnorm  = []
VyInfnorm  = []
Vy2ndnorm  = []
for i in range(nMax):
    rhoInfnorm.append(la.norm(np.absolute(rhoAll[i] - rhoAll[nMax]), np.inf))
    rho2ndnorm.append(la.norm(np.absolute(rhoAll[i] - rhoAll[nMax]), 2))
    pInfnorm.append  (la.norm(np.absolute(pAll  [i] - pAll  [nMax]), np.inf))
    p2ndnorm.append  (la.norm(np.absolute(pAll  [i] - pAll  [nMax]), 2))
    TInfnorm.append  (la.norm(np.absolute(TAll  [i] - TAll  [nMax]), np.inf))
    T2ndnorm.append  (la.norm(np.absolute(TAll  [i] - TAll  [nMax]), 2))
    VxInfnorm.append (la.norm(np.absolute(VxAll [i] - VxAll [nMax]), np.inf))
    Vx2ndnorm.append (la.norm(np.absolute(VxAll [i] - VxAll [nMax]), 2))
    VyInfnorm.append (la.norm(np.absolute(VyAll [i] - VyAll [nMax]), np.inf))
    Vy2ndnorm.append (la.norm(np.absolute(VyAll [i] - VyAll [nMax]), 2))

#----------------------------------------------------------
#-                   --------------------------------------
#-Error scatter plots--------------------------------------
#-                   --------------------------------------
#----------------------------------------------------------
cm   = plt.cm.get_cmap('RdYlBu')
pRho = input("Plot Density error scatter? Y or N: ")
if pRho == "Y" or pRho == "y":
    eRho = []
    for i in range(nMax):
        eRho.append(np.absolute(rhoAll[i] - rhoAll[nMax]))
        sRho = plt.scatter(xList, yList, c=eRho[i], vmin=0, vmax=max(eRho[0]), cmap=cm)
        plt.colorbar(sRho)
        plt.title("Density Error at Grid Level "+str(i))
        plt.figure()
    plt.show()
pP   = input("Plot Pressure error scatter? Y or N: ")
if pP == "Y" or pP == "y":
    eP   = []
    for i in range(nMax):
        eP.append(np.absolute(pAll[i] - pAll[nMax]))
        sP = plt.scatter(xList, yList, c=eP[i], vmin=0, vmax=max(eP[0]), cmap=cm)
        plt.colorbar(sP)
        plt.title("Pressure Error at Grid Level "+str(i))
        plt.figure()
    plt.show()
pT   = input("Plot Temperature error scatter? Y or N: ")
if pT == "Y" or pT == "y":
    eT   = []
    for i in range(nMax):
        eT.append(np.absolute(TAll[i] - TAll[nMax]))
        sT = plt.scatter(xList, yList, c=eT[i], vmin=0, vmax=max(eT[0]), cmap=cm)
        plt.colorbar(sT)
        plt.title("Temperature Error at Grid Level "+str(i))
        plt.figure()
    plt.show()
pVx  = input("Plot x-Velocity error scatter? Y or N: ")
if pVx == "Y" or pVx == "y":
    eVx  = []
    for i in range(nMax):
        eVx.append(np.absolute(VxAll[i] - VxAll[nMax]))
        sVx = plt.scatter(xList, yList, c=eVx[i], vmin=0, vmax=max(eVx[0]), cmap=cm)
        plt.colorbar(sVx)
        plt.title("x-Velocity Error at Grid Level "+str(i))
        plt.figure()
    plt.show()

#----------------------------------------------------------
#-              -------------------------------------------
#-Compute orders-------------------------------------------
#-              -------------------------------------------
#----------------------------------------------------------
rhoInfq = []
rho2ndq = []
pInfq   = []
p2ndq   = []
TInfq   = []
T2ndq   = []
VxInfq  = []
Vx2ndq  = []
VyInfq  = []
Vy2ndq  = []
for i in range(1,nMax):
    rhoInfq.append(-(np.log(rhoInfnorm[i]) - np.log(rhoInfnorm[i-1]))/np.log(ratio))
    rho2ndq.append(-(np.log(rho2ndnorm[i]) - np.log(rho2ndnorm[i-1]))/np.log(ratio))
    pInfq.append  (-(np.log(pInfnorm  [i]) - np.log(pInfnorm  [i-1]))/np.log(ratio))
    p2ndq.append  (-(np.log(p2ndnorm  [i]) - np.log(p2ndnorm  [i-1]))/np.log(ratio))
    TInfq.append  (-(np.log(TInfnorm  [i]) - np.log(TInfnorm  [i-1]))/np.log(ratio))
    T2ndq.append  (-(np.log(T2ndnorm  [i]) - np.log(T2ndnorm  [i-1]))/np.log(ratio))
    VxInfq.append (-(np.log(VxInfnorm [i]) - np.log(VxInfnorm [i-1]))/np.log(ratio))
    Vx2ndq.append (-(np.log(Vx2ndnorm [i]) - np.log(Vx2ndnorm [i-1]))/np.log(ratio))
    VyInfq.append (-(np.log(VyInfnorm [i]) - np.log(VyInfnorm [i-1]))/np.log(ratio))
    Vy2ndq.append (-(np.log(Vy2ndnorm [i]) - np.log(Vy2ndnorm [i-1]))/np.log(ratio))

#----------------------------------------------------------
#-                   --------------------------------------
#-Plot infinity norms--------------------------------------
#-                   --------------------------------------
#----------------------------------------------------------
f = plt.figure(figsize=(11.69,8.27))
# Density
plt.subplot(221)
plt.plot(Nx, rhoInfnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Density')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Pressure
plt.subplot(222)
plt.plot(Nx, pInfnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Pressure')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Temperature
plt.subplot(223)
plt.plot(Nx, TInfnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Temperature')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# X-velocity
plt.subplot(224)
plt.plot(Nx, VxInfnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_\infty}$', rotation=0, labelpad=10)
plt.title('X-velocity')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Y-velocity
# plt.subplot(224)
# plt.plot(Nx, VyInfnorm, '-ko')
# plt.xlim(NxLim)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\sqrt{N_c}$')
# plt.ylabel(r'$E_{l_\infty}$', rotation=0, labelpad=10)
# plt.title('Y-velocity')
# plt.grid(b=True, which='major', linestyle='-', linewidth=2)
# plt.grid(b=True, which='minor', linestyle='--')
plt.subplots_adjust(wspace=0.2, hspace=0.4)
# plt.show()
f.savefig('gc_study_infE.pdf', bbox_inches='tight')
#----------------------------------------------------------
#-                    -------------------------------------
#-Plot infinity orders-------------------------------------
#-                    -------------------------------------
#----------------------------------------------------------
g = plt.figure(figsize=(11.69,8.27))
# Density
plt.subplot(221)
plt.plot(Nx[1:nMax], rhoInfq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Density')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Pressure
plt.subplot(222)
plt.plot(Nx[1:nMax], pInfq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Pressure')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Temperature
plt.subplot(223)
plt.plot(Nx[1:nMax], TInfq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_\infty}$', rotation=0, labelpad=10)
plt.title('Temperature')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# X-velocity
plt.subplot(224)
plt.plot(Nx[1:nMax], VxInfq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_\infty}$', rotation=0, labelpad=10)
plt.title('X-velocity')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Y-velocity
# plt.subplot(224)
# plt.plot(Nx[1:nMax], VyInfq, '-ko')
# plt.xlim(NxLim)
# plt.xscale('log')
# plt.xlabel(r'$\sqrt{N_c}$')
# plt.ylabel(r'$q_{l_\infty}$', rotation=0, labelpad=10)
# plt.title('Y-velocity')
# plt.grid(b=True, which='major', linestyle='-', linewidth=2)
# plt.grid(b=True, which='minor', linestyle='--')
plt.subplots_adjust(wspace=0.2, hspace=0.4)
# plt.show()
g.savefig('gc_study_infq.pdf', bbox_inches='tight')
#----------------------------------------------------------
#-                 ----------------------------------------
#-Plot second norms----------------------------------------
#-                 ----------------------------------------
#----------------------------------------------------------
h = plt.figure(figsize=(11.69,8.27))
# Density
plt.subplot(221)
plt.plot(Nx, rho2ndnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_2}$', rotation=0, labelpad=10)
plt.title('Density')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Pressure
plt.subplot(222)
plt.plot(Nx, p2ndnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_2}$', rotation=0, labelpad=10)
plt.title('Pressure')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Temperature
plt.subplot(223)
plt.plot(Nx, T2ndnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_2}$', rotation=0, labelpad=10)
plt.title('Temperature')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# X-velocity
plt.subplot(224)
plt.plot(Nx, Vx2ndnorm, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$E_{l_2}$', rotation=0, labelpad=10)
plt.title('X-velocity')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Y-velocity
# plt.subplot(224)
# plt.plot(Nx, Vy2ndnorm, '-ko')
# plt.xlim(NxLim)
# plt.xscale('log')
# plt.yscale('log')
# plt.xlabel(r'$\sqrt{N_c}$')
# plt.ylabel(r'$E_{l_2}$', rotation=0, labelpad=10)
# plt.title('Y-velocity')
# plt.grid(b=True, which='major', linestyle='-', linewidth=2)
# plt.grid(b=True, which='minor', linestyle='--')
plt.subplots_adjust(wspace=0.2, hspace=0.4)
# plt.show()
h.savefig('gc_study_2ndE.pdf', bbox_inches='tight')
#----------------------------------------------------------
#-                  ---------------------------------------
#-Plot second orders---------------------------------------
#-                  ---------------------------------------
#----------------------------------------------------------
i = plt.figure(figsize=(11.69,8.27))
# Density
plt.subplot(221)
plt.plot(Nx[1:nMax], rho2ndq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_2}$', rotation=0, labelpad=10)
plt.title('Density')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Pressure
plt.subplot(222)
plt.plot(Nx[1:nMax], p2ndq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_2}$', rotation=0, labelpad=10)
plt.title('Pressure')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Temperature
plt.subplot(223)
plt.plot(Nx[1:nMax], T2ndq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_2}$', rotation=0, labelpad=10)
plt.title('Temperature')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# X-velocity
plt.subplot(224)
plt.plot(Nx[1:nMax], Vx2ndq, '-ko')
plt.xlim(NxLim)
plt.xscale('log')
plt.xlabel(r'$\sqrt{N_c}$')
plt.ylabel(r'$q_{l_2}$', rotation=0, labelpad=10)
plt.title('X-velocity')
plt.grid(b=True, which='major', linestyle='-', linewidth=2)
plt.grid(b=True, which='minor', linestyle='--')
# Y-velocity
# plt.subplot(224)
# plt.plot(Nx[1:nMax], Vy2ndq, '-ko')
# plt.xlim(NxLim)
# plt.xscale('log')
# plt.xlabel(r'$\sqrt{N_c}$')
# plt.ylabel(r'$q_{l_2}$', rotation=0, labelpad=10)
# plt.title('Y-velocity')
# plt.grid(b=True, which='major', linestyle='-', linewidth=2)
# plt.grid(b=True, which='minor', linestyle='--')
plt.subplots_adjust(wspace=0.2, hspace=0.4)
# plt.show()
i.savefig('gc_study_2ndq.pdf', bbox_inches='tight')
