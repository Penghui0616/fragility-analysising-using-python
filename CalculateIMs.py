#-*-coding: UTF-8-*-
import numpy as np
import math
import os
import IMs
import xlrd
import xlwt

dt=np.loadtxt("deltaT.txt")
m=len(dt)

cwd=os.getcwd()
if not os.path.exists(cwd+'/IM'):
    os.mkdir('IM')


pga = []
pgv = []
pgd = []
vda = []
aRMS = []
vRMS = []
dRMS = []
ai = []
ic = []
sed = []
cav = []
asi = []
vsi = []
sma = []
smv = []
eda = []
a95 = []
tpre = []
tm = []
td = []
ia = []
fi = []
iv = []
iD = []
sa02 = []
sa10 = []

for i1 in range(m):
    cwd=os.getcwd()
    pathall=os.path.join(cwd,'FiltedAcceleration/',str(i1+1)+".out")
    txtopen=np.loadtxt(pathall)
    imInstance=IMs.IMs(txtopen,dt[i1])

    pga.append(imInstance.PGA())
    pgv.append(imInstance.PGV())
    pgd.append(imInstance.PGD())
    vda.append(imInstance.VmaxDivAmax())
    aRMS.append(imInstance.aRMS())
    vRMS.append(imInstance.vRMS())
    dRMS.append(imInstance.dRMS())
    ai.append(imInstance.AI())
    ic.append(imInstance.Ic())
    sed.append(imInstance.SED())
    cav.append(imInstance.CAV())
    asi.append(imInstance.ASI())
    vsi.append(imInstance.VSI())
    sma.append(imInstance.SMA())
    smv.append(imInstance.SMV())
    eda.append(imInstance.EDA())
    a95.append(imInstance.A95())
    tpre.append(imInstance.Tpre())
    tm.append(imInstance.Tm())
    td_single, ia_single=imInstance.Ia()
    td.append(td_single)
    ia.append(ia_single)
    fi.append(imInstance.FI())
    iv.append(imInstance.Iv())
    iD.append(imInstance.Id())
    T02_single,sa02_single,sv02_single,sd02_single=imInstance.SpPiecewise (1,0.2,0.05)
    T10_single,sa10_single,sv10_single,sd10_single=imInstance.SpPiecewise (1,1,0.05)
    sa02.append(sa02_single)
    sa10.append(sa10_single)

np.savetxt(cwd+'\\IM\\PGA.txt', np.array(pga).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\PGV.txt', np.array(pgv).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\PGD.txt', np.array(pgd).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\VdivA.txt', np.array(vda).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\aRMS.txt', np.array(aRMS).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\vRMS.txt', np.array(vRMS).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\dRMS.txt', np.array(dRMS).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\AI.txt', np.array(ai).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\IC.txt', np.array(ic).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\SED.txt', np.array(sed).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\CAV.txt', np.array(cav).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\ASI.txt', np.array(asi).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\VSI.txt', np.array(vsi).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\SMA.txt', np.array(sma).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\SMV.txt', np.array(smv).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\EDA.txt', np.array(eda).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\A95.txt', np.array(a95).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Tpre.txt', np.array(tpre).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Tm.txt', np.array(tm).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Td.txt', np.array(td).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Ia.txt', np.array(ia).T, fmt='%f')
np.savetxt(cwd+'\\IM\\FI.txt', np.array(fi).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Iv.txt', np.array(iv).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Id.txt', np.array(iD).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Sa02.txt', np.array(sa02).T, fmt='%f') 
np.savetxt(cwd+'\\IM\\Sa10.txt', np.array(sa10).T, fmt='%f')