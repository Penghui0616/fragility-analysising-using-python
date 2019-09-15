#!/usr/bin/env python 
# -*- coding:utf-8 -*-
# @Author   : Penghui Zhang
# @Email    : penghui@tongji.edu.cn
import numpy as np
import math
import os
import regressionpy
import fragilitycurve
import systemfragility
import matplotlib as mpl
import matplotlib.pyplot as plt

plt.rcParams.update({'figure.max_open_warning': 0}) #解除plt打开图片数量限制

cwd=os.getcwd()
if not os.path.exists(cwd+'/regression'):
    os.mkdir('regression')
if not os.path.exists(cwd+'/fagilitycurve'):
    os.mkdir('fagilitycurve')

im = 'PGA'
edps = ['Abutment2_adispa', 'Abutment2_adispp', 'Bearing6_bdisp', 'Pier1_cdrift']
#edps = ['Abutment2_adispa', 'Abutment2_adispp']  
np.random.seed(6) #随机种子

with open(cwd+'/postprocessing/convergence.txt','r') as f:
    convergence =  np.array(f.read().split('\n'))
deleteIndex = np.argwhere(convergence=='No')


gmNumber = int(len(np.delete(convergence, deleteIndex))-1)
edpNumber = int(len(edps))

originalLnEDP = np.zeros((gmNumber,edpNumber))
originalLnEDPbc = np.zeros((gmNumber,edpNumber))

lmCoeffictionList = []
lmInterceptList = []

results = os.listdir('postprocessing')

for i in range(len(edps)):
    for result in results:
        if edps[i] in result:
            componentEDP = result
    lnIM = np.array([math.log(x) for x in np.loadtxt(cwd+'/IM/'+im+'.txt')])
    lnEDP = np.array([math.log(x) for x in np.loadtxt(cwd+'/postprocessing/'+componentEDP)])

    #删除不收敛的地震波
    lnIM = np.delete(lnIM, deleteIndex).copy()
    lnEDP = np.delete(lnEDP, deleteIndex).copy()
    originalLnEDP[:,i] = lnEDP.copy()

    #线性回归（作为比较）
    lmRegression = regressionpy.LinerRegression('ln({})'.format(im),'ln({})'.format(edps[i]))
    lmBeltad = lmRegression.scatterPlot(lnIM, lnEDP)
    lmRegression.residualPlot(lnIM, lnEDP)
    lmRegression.residualHist(lnIM, lnEDP)
    lmCoeffiction = lmRegression.coeffiction
    lmIntercept = lmRegression.intercept
    lmCoeffictionList.append(lmCoeffiction[0])
    lmInterceptList.append(lmIntercept[0])

    #画构件易损性曲线曲线
    immin = math.pow(math.e, np.min(lnIM))
    immax = math.pow(math.e, np.max(lnIM))
    imRange = [0.001, 1.00]
    fragility = fragilitycurve.FragilityCurve(im, edps[i])
    fragility.fragilityCurvePlot (imRange,lmIntercept, lmCoeffiction, lmBeltad)

#画体系易损性曲线
immin = math.pow(math.e, np.min(lnIM))
immax = math.pow(math.e, np.max(lnIM))
imRange = [0.001, 1]
systemFragility = systemfragility.SystemFragilityCurve(im, edps)
systemFragility.systemFragilityCurvePlot (imRange, originalLnEDP, lmCoeffictionList, lmInterceptList)

print(str(gmNumber) + ' ground motions are convergent.')



