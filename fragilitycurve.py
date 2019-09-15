#!/usr/bin/env python 
# -*- coding:utf-8 -*-
# @Author   : Penghui Zhang
# @Email    : penghui@tongji.edu.cn

import matplotlib.pyplot as plt
import numpy as np
import math

class FragilityCurve:
    def __init__ (self, im, edp):
        self.im = im
        self.edp = edp
        self.n = 100000 #蒙特卡洛抽样次数
        self.dot = 30 #易损性曲线点数
 
    def capcityModel (self):
        #能力模型生成100000行的抽样矩阵
        capacity = None
        if self.edp.split('_')[1]=='cdrift':
            sc = np.array([0.005, 0.01, 0.02, 0.025])
            beltac = np.array([0.0025, 0.0025, 0.0046, 0.0046])
            capacity = np.zeros((self.n, 4))
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
                

        elif self.edp.split('_')[1]=='bdisp':
            sc = np.array([0.15, 0.35])
            beltac = np.array([0.35, 0.35])
            capacity = np.zeros((self.n, 2))
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None)])
                    if randomArray[0]<randomArray[1]:
                        capacity[i,:] = randomArray.copy()
                        flag = 0
                

        elif self.edp.split('_')[1]=='adispa':
            sc = np.array([0.01, 0.038, 0.077])
            beltac = np.array([0.0007, 0.0009, 0.00085])
            capacity = np.zeros((self.n, 3))
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None),\
                                            np.random.lognormal(mean=math.log(sc[2]), sigma=beltac[2], size=None)])
                    if randomArray[0]<randomArray[1]<randomArray[2]:
                        capacity[i,:] = randomArray.copy()
                        flag = 0
                

        elif self.edp.split('_')[1]=='adispp':
            sc = np.array([0.037, 0.147])
            beltac = np.array([0.00046, 0.00046])
            capacity = np.zeros((self.n, 2))
            for i in range(self.n):
                flag = 1
                while flag:
                    randomArray = np.array([np.random.lognormal(mean=math.log(sc[0]), sigma=beltac[0], size=None),\
                                            np.random.lognormal(mean=math.log(sc[1]), sigma=beltac[1], size=None)])
                    if randomArray[0]<randomArray[1]:
                        capacity[i,:] = randomArray.copy()
                        flag = 0
        return capacity

    def lmDemandModel (self, lnIM, intercept, coeffiction, beltad):
        lnSd = intercept + coeffiction*lnIM
        lmDemand = np.random.lognormal(mean=lnSd, sigma=beltad, size=(self.n,1))
        return lmDemand

    def mmToInches (self,mm):
        #mm transform to inches
        inches=mm*0.0393700787
        return inches

    def fragilityCurvePlot (self, imRange,lmIntercept, lmCoeffiction, lmBeltad):

        imList = np.linspace(imRange[0], imRange[1], self.dot)
        lnIMList = np.array([math.log(x) for x in imList])

        if self.edp.split('_')[1] == 'cdrift':
            lmfragility = np.zeros((self.dot, 4))
        elif self.edp.split('_')[1] == 'bdisp':
            lmfragility = np.zeros((self.dot, 2))
        elif self.edp.split('_')[1] == 'adispa':
            lmfragility = np.zeros((self.dot, 3)) 
        elif self.edp.split('_')[1] == 'adispp':
            lmfragility = np.zeros((self.dot, 2))

        for i in range(self.dot):
            lnIM = lnIMList[i]
            capacity = self.capcityModel()
            lmdemand = self.lmDemandModel (lnIM, lmIntercept, lmCoeffiction, lmBeltad)
            lmfragility[i,:] = np.sum(lmdemand>capacity, axis=0)/self.n
        print('lmfragility: ',lmfragility[i,:])

        #画出易损性曲线
        width=self.mmToInches(70)
        height=self.mmToInches(50)
        fig = plt.figure(facecolor="white", figsize=(width, height))

        if self.edp.split('_')[1] == 'cdrift':
            plt.plot(imList,lmfragility[:,0],color='red',linestyle=':',linewidth=1,label='LS1')
            plt.plot(imList,lmfragility[:,1],color='red',linestyle='-.',linewidth=1,label='LS2')
            plt.plot(imList,lmfragility[:,2],color='red',linestyle='--',linewidth=1,label='LS3')
            plt.plot(imList,lmfragility[:,3],color='red',linestyle='-',linewidth=1,label='LS4')

        elif self.edp.split('_')[1] == 'bdisp':
            plt.plot(imList,lmfragility[:,0],color='red',linestyle=':',linewidth=1,label='LS1')
            plt.plot(imList,lmfragility[:,1],color='red',linestyle='-.',linewidth=1,label='LS2')

        elif self.edp.split('_')[1] == 'adispa':
            plt.plot(imList,lmfragility[:,0],color='red',linestyle=':',linewidth=1,label='LS1')
            plt.plot(imList,lmfragility[:,1],color='red',linestyle='-.',linewidth=1,label='LS2')
            plt.plot(imList,lmfragility[:,2],color='red',linestyle='--',linewidth=1,label='LS3')

        elif self.edp.split('_')[1] == 'adispp':
            plt.plot(imList,lmfragility[:,0],color='red',linestyle=':',linewidth=1,label='LS1')
            plt.plot(imList,lmfragility[:,1],color='red',linestyle='-.',linewidth=1,label='LS2')


        plt.xlabel(self.im,size=8)
        plt.ylabel('probability',size=8) 
        plt.xlim(imRange[0], imRange[1])
        plt.ylim(0, 1)
        plt.legend(loc='lower right',frameon=True,edgecolor='black',fontsize=6)

        ax=plt.gca()
        ax.tick_params(direction='in')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]

        plt.savefig('fagilitycurve/'+self.edp+" FragilityCurve.png",dpi = 960, bbox_inches="tight") 
        plt.savefig('fagilitycurve/'+self.edp+" FragilityCurve.eps",dpi = 960, bbox_inches="tight") 



