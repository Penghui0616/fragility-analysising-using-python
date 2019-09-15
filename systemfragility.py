#!/usr/bin/env python 
# -*- coding:utf-8 -*-
# @Author   : Penghui Zhang
# @Email    : penghui@tongji.edu.cn

import numpy as np
import math
import os
import matplotlib.pyplot as plt

class SystemFragilityCurve:
    def __init__ (self, im, edps):
        self.im = im
        self.edps = edps
        self.n = 100000 #蒙特卡洛抽样次数
        self.dot = 30 #易损性曲线点数

    def lmDemandModel (self, originalLnEDP, lnIM, lmCoeffictionList, lmInterceptList):
        covMtrix = np.cov(originalLnEDP.T)
        lnSdmean = np.array(lmCoeffictionList) * lnIM + np.array(lmInterceptList)
        lnDemand = np.random.multivariate_normal(mean=lnSdmean, cov=covMtrix, size=self.n)
        lmDemand = pow(math.e, lnDemand)
        return lmDemand
    
    def componentCapcityModel (self, edp):
        #对单个构件的能力抽样
        capacity = np.ones((self.n, 4))*(10**8)
        if edp.split('_')[1]=='cdrift':
            sc = np.array([0.005, 0.01, 0.02, 0.025])
            beltac = np.array([0.0025, 0.0025, 0.0046, 0.0046])
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None),\
                                            np.random.lognormal(mean=math.log(sc[2]), sigma=beltac[2], size=None),\
                                            np.random.lognormal(mean=math.log(sc[3]), sigma=beltac[3], size=None)])
                    if randomArray[0]<randomArray[1]<randomArray[2]<randomArray[3]:
                        capacity[i,:] = randomArray.copy()
                        flag = 0

        elif edp.split('_')[1]=='bdisp':
            sc = np.array([0.15, 0.35])
            beltac = np.array([0.35, 0.35])
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None)])
                    if randomArray[0]<randomArray[1]:
                        capacity[i,[0,1]] = randomArray.copy()
                        flag = 0
                

        elif edp.split('_')[1]=='adispa':
            sc = np.array([0.01, 0.038, 0.077])
            beltac = np.array([0.0007, 0.0009, 0.00085])
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None),\
                                            np.random.lognormal(mean=math.log(sc[2]), sigma=beltac[2], size=None)])
                    if randomArray[0]<randomArray[1]<randomArray[2]:
                        capacity[i,[0,1,2]] = randomArray.copy()
                        flag = 0
                

        elif edp.split('_')[1]=='adispp':
            sc = np.array([0.037, 0.147])
            beltac = np.array([0.00046, 0.00046])
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None)])
                    if randomArray[0]<randomArray[1]:
                        capacity[i,[0,1]] = randomArray.copy()
                        flag = 0
        return capacity

    def capcityModel (self):
        #生成3维的能力抽样矩阵
        numComponent = len(self.edps)
        systemCapacity = np.zeros((numComponent, self.n, 4))
        for i in range(numComponent):
            systemCapacity[i,:,:] = self.componentCapcityModel (self.edps[i])
        return systemCapacity

    def mmToInches (self,mm):
        #mm transform to inches
        inches=mm*0.0393700787
        return inches

    def systemFragilityCurvePlot (self, imRange, originalLnEDP, lmCoeffictionList, lmInterceptList):
        #画出体系易损性曲线
        imList = np.linspace(imRange[0], imRange[1], self.dot)
        lnIMList = np.array([math.log(x) for x in imList])
        lmFragilityLS1 = []
        lmFragilityLS2 = []
        lmFragilityLS3 = []
        lmFragilityLS4 = []

        for i in range(self.dot):
            lnIM = lnIMList[i]
            systemCapacity = self.capcityModel()
            lmDemand = self.lmDemandModel (originalLnEDP, lnIM, lmCoeffictionList, lmInterceptList)
            lmFragilityLS1.append(np.sum(np.sum(lmDemand>systemCapacity[:,:,0].T, axis=1)>0)/self.n)
            lmFragilityLS2.append(np.sum(np.sum(lmDemand>systemCapacity[:,:,1].T, axis=1)>0)/self.n)
            lmFragilityLS3.append(np.sum(np.sum(lmDemand>systemCapacity[:,:,2].T, axis=1)>0)/self.n)
            lmFragilityLS4.append(np.sum(np.sum(lmDemand>systemCapacity[:,:,3].T, axis=1)>0)/self.n)


        #画出体系易损性曲线
        width=self.mmToInches(70)
        height=self.mmToInches(50)
        fig = plt.figure(facecolor="white", figsize=(width, height))

        plt.plot(imList,lmFragilityLS1,color='red',linestyle=':',linewidth=1,label='LS1')
        plt.plot(imList,lmFragilityLS2,color='red',linestyle='-.',linewidth=1,label='LS2')
        plt.plot(imList,lmFragilityLS3,color='red',linestyle='--',linewidth=1,label='LS3')
        plt.plot(imList,lmFragilityLS4,color='red',linestyle='-',linewidth=1,label='LS4')


        plt.xlabel(self.im,size=8)
        plt.ylabel('probaility',size=8) 
        plt.xlim(imRange[0], imRange[1])
        plt.ylim(0, 1)
        plt.legend(loc='lower right',frameon=True,edgecolor='black',fontsize=6)

        ax=plt.gca()
        ax.tick_params(direction='in')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

        plt.savefig('fagilitycurve/SystemFragilityCurve.png',dpi = 960, bbox_inches="tight") 
        plt.savefig('fagilitycurve/SystemFragilityCurve.eps',dpi = 960, bbox_inches="tight") 
