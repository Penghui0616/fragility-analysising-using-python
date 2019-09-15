#!/usr/bin/env python 
# -*- coding:GBK -*-
# @Author   : Penghui Zhang
# @Email    : penghui@tongji.edu.cn
import matplotlib.pyplot as plt
import numpy as np
import math
import seaborn as sns
from scipy import stats
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

class LinerRegression:
    def __init__ (self, xsymbol, ysymbol):
        """初始化Linear Regression模型"""
        self.coeffiction = None
        self.intercept = None
        self.xsymbol = xsymbol
        self.ysymbol = ysymbol
  

    def linearRegression (self, x, y):
        #一般线性回归，data为数据，其中最后一列为y值，x的排序按xo,x1,x2....
        #y=xT*w,其中x为变量列矩阵，w为回归系数列矩阵
        #error=sum((yi-xi*w.T)**2)
        #w=inv(xTx)*xT*y
        #return:w-regression coeffient coeff-correlation coeffient R2-coeffient of determination
        hstackx = np.hstack([np.ones((x.shape[0], 1)), x])
        xMat=np.mat(hstackx)
        yMat=np.mat(y)
        xTx=xMat.T*xMat
        if np.linalg.det(xTx)==0.0:
            print ("Matrix singularity, can not be inverted! ")
            return
        else:
            w=xTx.I*xMat.T*yMat
        yEstimate=xMat*w
        self.intercept = np.array(w)[0].copy()
        self.coeffiction = np.array(w)[1].copy()

        #compute coeffient of determination R2
        r2, sse=self.deterCoeff(yEstimate, y)

        #计算标准差
        beltad = math.sqrt(sse/(x.shape[0]-2))
 
        return  r2, beltad, yEstimate


    def deterCoeff (self,yEstimate,yRecorder):
        #compute coeffient of determination R2
        #input:yEstimate--回归得到的预测值（Nx1）列向量
                ##yRecorder--实验得到的值（Nx1）列向量
        ##output:R2, sse
        ymean=float(sum(yRecorder))/len(yRecorder)
        sst=sum(np.array(np.array(yRecorder)-ymean)**2)
        sse=sum(np.array(np.array(yRecorder)-np.array(yEstimate))**2)
        r2=1-sse/sst
        return r2, sse

    def rank (self, y):
        sorty = np.argsort(y[:,0])
        ranky = np.zeros(y.shape[0])
        for i in range(y.shape[0]):
            ranky[sorty[i]] = i
        return np.array(ranky)

        
    def rankcoeff (self, x, e):
        #计算等级相关系数及p_value
        rankx = self.rank(x)
        ranke = self.rank(abs(e))
        n = e.shape[0]
        rs = 1-np.sum((rankx-ranke)**2)*6/(n*(n**2-1))
        t = math.sqrt(n-2)*rs/math.sqrt(1-rs**2)
        pValue = (1-stats.t.cdf(abs(t), n-2))*2

        return rs, pValue
    

    def mmToInches (self,mm):
        #mm transform to inches
        inches=mm*0.0393700787
        return inches

    def scatterPlot (self, x, y):
        xCol = x.reshape(x.shape[0], 1)
        yCol = y.reshape(y.shape[0], 1)
        r2, beltad, yEstimate = self.linearRegression (xCol, yCol)
        minx, maxx = np.min(x), np.max(x)
        miny, maxy = np.min(y), np.max(y)
        xliner = np.linspace(minx, maxx, 1000)
        yliner = np.array(xliner*self.coeffiction + self.intercept*np.ones(1000))
    
        width=self.mmToInches(60)
        height=self.mmToInches(50)
        fig = plt.figure(facecolor="white", figsize=(width, height))
        plt.scatter(x, y, facecolors='none', edgecolors='b',s=40)
        plt.plot(xliner,yliner,color='r',linestyle='-',linewidth=2)
        plt.xlabel(self.xsymbol,size=8)
        plt.ylabel(self.ysymbol,size=8) 
        plt.xlim(minx-0.1*(maxx-minx), maxx+0.1*(maxx-minx))
        plt.ylim(miny-0.1*(maxy-miny), maxy+0.1*(maxy-miny))
        plt.text(minx,maxy-0.05*(maxy-miny),'R_square = {}'.format(str("%.3f" % r2)),size=7)
        plt.text(minx,maxy-0.15*(maxy-miny),'belta_D = {}'.format(str("%.3f" % beltad)),size=7)

        ax=plt.gca()
        ax.tick_params(direction='in')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.savefig('regression/'+self.ysymbol+" Regression.png",dpi = 960, bbox_inches="tight") 
        plt.savefig('regression/'+self.ysymbol+" Regression.eps",dpi = 960, bbox_inches="tight")
        return beltad

    def residualPlot (self, x, y):
        xCol = x.reshape(x.shape[0], 1)
        yCol = y.reshape(y.shape[0], 1)
        r2, beltad, yEstimate = self.linearRegression (xCol, yCol)
        e = yCol - np.array(yEstimate)
        rs, pValue = self.rankcoeff(xCol, e)

        minx, maxx = np.min(x), np.max(x)
        mine, maxe = np.min(e), np.max(e)
    
        width=self.mmToInches(60)
        height=self.mmToInches(50)
        fig = plt.figure(facecolor="white", figsize=(width, height))
        plt.scatter(x, e, facecolors='none', edgecolors='b',s=40)
        
        plt.xlabel(self.xsymbol,size=8)
        plt.ylabel('e',size=8) 
        plt.xlim(minx-0.1*(maxx-minx), maxx+0.1*(maxx-minx))
        plt.ylim(mine-0.1*(maxe-mine), maxe+0.1*(maxe-mine))
        plt.text(minx-0.05*(maxx-minx),mine+0.04*(maxe-mine),'rs={}'.format(str("%.3f" % rs)),size=7)
        plt.text(minx-0.05*(maxx-minx),mine-0.04*(maxe-mine),'p_value={}'.format(str("%.2f" % (pValue))),size=7)

        ax=plt.gca()
        ax.tick_params(direction='in')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        plt.savefig('regression/'+self.ysymbol+" Residual.png",dpi = 960, bbox_inches="tight") 
        plt.savefig('regression/'+self.ysymbol+" Residual.eps",dpi = 960, bbox_inches="tight")
        
    def residualHist(self, x, y):
        xCol = x.reshape(x.shape[0], 1)
        yCol = y.reshape(y.shape[0], 1)
        r2, beltad, yEstimate = self.linearRegression (xCol, yCol)
        e = yCol - np.array(yEstimate)
        mine, maxe = np.min(e), np.max(e)

        width=self.mmToInches(20)
        height=self.mmToInches(50)
        fig = plt.figure(facecolor="white", figsize=(width, height))
#        plt.hist(e, bins=9, density=True, alpha=0.5, histtype='stepfilled',color='blue')
        sns.distplot(e,color="blue",bins=7,kde=True, kde_kws={"label":"KDE","color":"green"},
                     fit=stats.norm, fit_kws={"label":"Norm","color":"red","linestyle":"--"},vertical=True)

        ax=plt.gca()
        ax.set_ylim(mine-0.1*(maxe-mine), maxe+0.1*(maxe-mine))
        ax.set_yticks([])
        ax.tick_params(which='both',direction='in')
        labels = ax.get_xticklabels() + ax.get_yticklabels()
        [label.set_fontname('Times New Roman') for label in labels]
        ax.legend(loc = 'upper right',frameon=True,edgecolor='black',fontsize=5) 

        plt.savefig('regression/'+self.ysymbol+" Hist.png",dpi = 960, bbox_inches="tight") 
        plt.savefig('regression/'+self.ysymbol+" Hist.eps",dpi = 960, bbox_inches="tight") 
