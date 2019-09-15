#-*-coding: UTF-8-*-
import numpy as np
import math


class IMs():
    def __init__(self,acc,t):
        #acc单列加速度时程(g),t时间间隔(s)
        self.acc=acc
        self.t=t
        self.num=len(self.acc)


    def PGA(self):
        #返回acc的最大峰值加速度(PGA)(g)
        pga=b=np.fabs(self.acc).max()
        return pga


    def AcctoVelocity (self):
        #将加速度(g)转换为速度(cm/s)
        vel=[0]
        acc=self.acc
        for i in range(self.num-1):
            velocity=(acc[i]+acc[i+1])*self.t/2*981+vel[-1]
            vel.append(velocity)
        return vel



    def PGV (self):
        #返回Velocity的最大峰值速度(PGV)(cm/s)
        veloc=self.AcctoVelocity()
        pgv=b=np.fabs(veloc).max()
        return pgv

    def AcctoDisp (self):
        #将加速度(g)转换为位移(cm)
        disp=[0]
        veld=self.veloc=self.AcctoVelocity()
        for i in range(self.num-1):
            displacement=(veld[i]+veld[i+1])*self.t/2+disp[-1]
            disp.append(displacement)
        return disp
    
        
    def PGD (self):
        #返回位移的最大峰值(PGD)(cm)
        dispc=self.AcctoDisp()
        pgd=np.fabs(dispc).max()
        return pgd


    def VmaxDivAmax (self):
        #Vmax/Amax比值(sec)
        vmax=self.PGV()
        amax=self.PGA()*981
        ratio=vmax/amax
        return ratio


    def aRMS (self):
        #aRMS=(integral(a(t)**2,0,ttot)/tot)**0.5(g)
        a2=np.array(self.acc**2)
        ttol=self.num*self.t
        atol=sum(a2)*self.t
        arms=(atol/ttol)**0.5
        return arms


    def vRMS (self):
        ##aRMS=(integral(v(t)**2,0,ttot)/tot)**0.5(cm/s)
        vel=np.array(self.AcctoVelocity())
        vel2=vel**2
        ttol=self.num*self.t
        vtol=sum(vel2)*self.t
        vrms=(vtol/ttol)**0.5
        return vrms


    def dRMS (self):
        ##dRMS=(integral(d(t)**2,0,ttot)/tot)**0.5(cm)
        disp=np.array(self.AcctoDisp())
        disp2=disp**2
        ttol=self.num*self.t
        dtol=sum(disp2)*self.t
        drms=(dtol/ttol)**0.5
        return drms


    def AI (self):
        #Aria Intensity Ia=(integral(a(t)**2,0,ttol)*pi/g)(m/s)
        a2=self.acc*9.81
        a3=np.array(a2)
        a3tol=sum(a3**2)*self.t
        ai=a3tol*math.pi/(2*9.81)
        return ai


    def Ic (self):
        #Characteristic intensity Ic=(aRMS)**(3/2)*(ttot)*0.5(g**1.5*s**0.5)
        ttol=self.num*self.t 
        Ic=(self.aRMS())**(1.5)*(ttol)**0.5
        return Ic


    def SED (self):
        #Specific Energy Density (SED) SED=integral(v(t)**2,0,ttol)(cm2/s)
        vel=self.AcctoVelocity()
        velarr=np.array(vel)
        sed=sum(velarr**2)*self.t
        return sed


    def CAV (self):
        #Cumulative Absolute Velocity (CAV)(cm/s)
        accarr=np.array(self.acc)
        accabs=np.fabs(accarr)*981
        cav=sum(accabs)*self.t
        return cav


    def DisptoVelocity (self,displacement,t):
        #将位移(cm)转换为速度(cm/s)
        #displacement-需要转换为速度的序列(cm)
        #t-时间间隔(s)
        n=len(displacement)
        vel=[0]
        disp=displacement
        for i in range(n-1):
            vell=2*(disp[i+1]-disp[i])/t-vel[-1]
            vel.append(vell)
        return vel


    def VeltoAccele (self,vel,t):
        #将速度(cm/s)转换为加速度(g)
        #vel-需要转换为加速度的序列(cm/s)
        #t-时间间隔(s)
        n=len(vel)
        acc=[0]
        for i in range(n-1):
            accel=2*(vel[i+1]-vel[i])/t-acc[-1]
            acc.append(accel)
        acceleration=np.array(acc)*(0.01/9.81)
        return acceleration
                


    def ASI (self):
        #Acceleration spectrum intensity(ASI) (g.second)
        #ASI=integral(Sa(dratio=0.05,T),0.1,0.5)
        sa=[]
        for i in np.arange(0.1,0.52,0.02):
            a=self.SpPiecewise (1,i,0.05)
            sa.append(a[1])
        asi=sum(sa)*0.02
        return asi


    def VSI (self):
        #Velocity spectrum intensity(VSI) (cm)
        #VSI=integral(Sv(dratio=0.05,T),0.1,2.5)
        sv=[]
        for i in np.arange(0.1,2.52,0.02):
            ainter=self.SpPiecewise (1,i,0.05)
            sv.append(ainter[2])
        vsi=sum(sv)*0.02
        return vsi


    def __SMR (self,response):
        #this parameter gives the sustained maximum acceleration/velocity during three cycles,
        #and is defined as the third highest absolute value of acceleration/velocity in the time-history
        #(note: in order for an absolute value to be considered as a "maximum", it must be larger than values 20 steps before and 20 steps after).
        #Nuttli O.W. [1979] "The relation of sustained maximum ground acceleration and velocity to earthquake intensity and magnitude,
        #" Miscellaneous Paper S-71-1, Report 16, U.S. Army Corps of Engineers, Waterways Experiment Station, Vicksburg, Mississippi.
        a=response
        aBefore=a[0:19]
        aAfter=a[-20:-1]
        aMiddle=a[20:-21]
        b=[]
        c=[]
        for i in range(len(aMiddle)-1):
            if aMiddle[i]*aMiddle[i+1]<0:
                b.append(i)
        for i in np.arange(0,len(b)-1,1):
            c.append(np.fabs(aMiddle[b[i]:b[i+1]]).max())
        c.sort()
        return c[-3]


    def SMA (self):
        ##Sustained maximum acceleration (SMA) (g)
        sma=self.__SMR(self.acc)
        return sma


    def SMV (self):
        ##Sustained maximum velocity (SMV) (cm)
        smv=self.__SMR(self.AcctoVelocity ())
        return smv

    def __Butter (self,n,flowcut,f):
        #Butterworth filter, f is the frequency (Hz)
#       scale=(((f/flowcut)**(2*n))/(1+(f/flowcut)**(2*n)))**0.5
#       return scale
        if f>=flowcut:
            scale=1
        else:
            scale=0
        return scale
        
        


    def lowCutFilter (self,n,flowcut):
        #flowcut filter,i.e. highfrequency pass,n=5, flowcut(Hz)
        npads=math.ceil(1.5*n/flowcut*0.5)
        padlist=[]
        for i in range(int(npads)):
            padlist.append(0)
        acc=padlist+list(self.acc)+padlist

        t=self.t
        fN=1/(2*t)
        xw=np.fft.fft(acc)
        xabs=abs(xw)
        xreal=np.real(xw)
        ximag=np.imag(xw)
        xy= (np.array(xreal)**2+np.array(ximag)**2)**0.5
        xcos=np.array(xreal)/xy
        xsin=np.array(ximag)/xy
        freq = np.fft.fftfreq(len(acc),self.t)

        posdata=[]
        negdata=[]
        for i in range(len(freq)):
            if freq[i]>=0:
                posdata.append(freq[i])
            else:
                negdata.append(np.fabs(freq[i]))

        posscale=[]
        negscale=[]
        for i in range(len(posdata)):
            posscale.append(self.__Butter(n,flowcut,posdata[i]))
        for i in range(len(negdata)):
            negscale.append(self.__Butter(n,flowcut,negdata[i]))


        xabsposdata=[]
        xabsnegdata=[]

        xabsposdata=xabs[:len(posdata)]
        xabsnegdata=xabs[len(posdata):]

        newposdata=[]
        newnegdata=[]
        for i in range(len(posdata)):
            newposdata.append(xabsposdata[i]*posscale[i])
        for i in range(len(negdata)):
            newnegdata.append(xabsnegdata[i]*negscale[i])
        newdata=newposdata+newnegdata

        newxw=[]
        for i in range(len(xreal)):
            newxw.append(complex(newdata[i]*xcos[i],newdata[i]*xsin[i]))
        ixw=np.array(newxw)
        newacc=np.real(np.fft.ifft(ixw))
        finalacc=newacc[int(npads):-int(npads)]
        return finalacc

    def highCutFilter (self,n,flowcut):
        #highcut filter,i.e. lowfrequency pass,N=5,flowcut(Hz)
        npads=math.ceil(1.5*n/flowcut*0.5)
        padlist=[]
        for i in range(int(npads)):
            padlist.append(0)
        acc=padlist+list(self.acc)+padlist

        t=self.t
        fN=1/(2*t)
        xw=np.fft.fft(acc)
        xabs=abs(xw)
        xreal=np.real(xw)
        ximag=np.imag(xw)
        xy= (np.array(xreal)**2+np.array(ximag)**2)**0.5
        xcos=np.array(xreal)/xy
        xsin=np.array(ximag)/xy
        freq = np.fft.fftfreq(len(acc),self.t)

        posdata=[]
        negdata=[]
        for i in range(len(freq)):
            if freq[i]>=0:
                posdata.append(freq[i])
            else:
                negdata.append(np.fabs(freq[i]))

        posscale=[]
        negscale=[]
        for i in range(len(posdata)):
            posscale.append(1-self.__Butter(n,flowcut,posdata[i]))
        for i in range(len(negdata)):
            negscale.append(1-self.__Butter(n,flowcut,negdata[i]))


        xabsposdata=[]
        xabsnegdata=[]

        xabsposdata=xabs[:len(posdata)]
        xabsnegdata=xabs[len(posdata):]

        newposdata=[]
        newnegdata=[]
        for i in range(len(posdata)):
            newposdata.append(xabsposdata[i]*posscale[i])
        for i in range(len(negdata)):
            newnegdata.append(xabsnegdata[i]*negscale[i])
        newdata=newposdata+newnegdata

        newxw=[]
        for i in range(len(xreal)):
            newxw.append(complex(newdata[i]*xcos[i],newdata[i]*xsin[i]))
        ixw=np.array(newxw)
        newacc=np.real(np.fft.ifft(ixw))
        finalacc=newacc[int(npads):-int(npads)]
        return finalacc

    def EDA (self):
        #This parameter corresponds to the peak acceleration value found after (g)
        #lowpass filtering the input time history with a cut-off frequency of 9 Hz [Benjamin and Associates, 1988].
        reacc=self.highCutFilter (5,9)
        eda=np.fabs(reacc).max()
        return eda
        


    def A95 (self):
        #The acceleration level below which 95% of the total Arias intensity is contained.
        #In other words, if the entire accelerogram yields a value of Ia equal to 100, the 
        #A95 parameter is the threshold of acceleration such that integrating all the values 
        #of the accelerogram below it, one gets an Ia=95. (g)
        ai=self.AI()
        threshold=0.95*ai
        acc=np.fabs(self.acc*9.81)
        t=self.t
        acc.sort()
        acclist=[0]
        for i in range(len(acc)):
            acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) 
        for i in range(len(acclist)):
            if acclist[i]>=threshold:
                a95=acc[i-1]/9.81
                break
        return a95


    def Tpre (self):
        #The predominant period (Tp) is the period at which the maximum spectral 
        #acceleration occurs in an acceleration response spectrum calculated at 5% damping.(s)
        spectrum=[]
        T=[]
        for i in np.arange(0.02,10,0.02):
            aspec=self.SpPiecewise  (1,i,0.05)
            spectrum.append(aspec[1])                              
            T.append(aspec[0])
        np.savetxt("spectrum.txt",spectrum,fmt="%f")
        spmax=max(spectrum)
        for i in range(len(spectrum)):
            if spectrum[i]==spmax:
                index=i
                break
        return T[i]


    def Tm (self):
        #According to Rathje et al. [1998] the mean period (Tm) is the best simplified frequency (s)
        #content characterisation parameter, being estimated with the following equation, where 
        #Ci are the Fourier amplitudes, and fi represent the discrete Fourier transform frequencies between 0.25 and 20 Hz.
        acc=self.acc
        xw=np.fft.fft(acc)
        xabs=abs(xw)
        xreal=np.real(xw)
        ximag=np.imag(xw)
        freq = np.fft.fftfreq(len(acc),self.t)
        for i in range(len(freq)):
            if freq[i]>=0.25:
                lowindex=i
                break
        for i in range(len(freq)):
            if freq[i]>=20:
                upperindex=i
                break
        citot=xabs[lowindex:upperindex]
        fitot=freq[lowindex:upperindex]
        cf=sum((np.array(citot)**2)/np.array(fitot))
        tm=cf/sum((np.array(citot)**2))
        return tm


    def DurationBrac(self):
        #Bracketed duration.The total time elapsed between the first and the last excursions of 
        #a specified level of acceleration (usually 0.05g).(Bolt 1969)
        acc=self.acc
        upperindex = 0
        lowindex = 0
        for i in range(len(acc)):
            if np.fabs(acc[i])>=0.05:
                lowindex=i
                break
        list(acc).reverse()
        for j in range(len(acc)):
            if np.fabs(acc[j])>=0.05:
                upperindex=len(acc)-j
                break
        t=(upperindex-lowindex)*self.t
        return t

    def SpPiecewise (self,m,T,dratio):
        #分段精确方法，m-质量(ton),T-周期(second),dratio-阻尼比
        # Sa(g)(绝对加速度=Sarel+Sg) Sv(cm/s) Sd(cm)
        acc=np.array(self.acc)*9.81
        t=self.t
        n=len(acc)
        w=2*math.pi/T
        wd=w*(1-dratio**2)**0.5
        k=m*(w**2)
        c=2*m*w*dratio
        vn=[0]
        vdotn=[0]
        v2dotn=[0]
        for i in range(n-1):
            alpa=-m*(acc[i+1]-acc[i])/t
            A=alpa*c/(k**2)-(-m*acc[i])/k+vn[-1]
            B=(w/wd)*dratio*A-alpa/(k*wd)+vdotn[-1]/wd
            v=math.exp(-dratio*w*t)*(A*math.cos(wd*t)+B*math.sin(wd*t))+(-m*acc[i]+alpa*t)/k-alpa*c/(k**2)
            vn.append(v)
            vdot=-dratio*w*math.exp(-dratio*w*t)*(A*math.cos(wd*t)+B*math.sin(wd*t))+math.exp(-dratio*w*t)\
                *(-A*wd*math.sin(wd*t)+B*wd*math.cos(wd*t))+alpa/k
            vdotn.append(vdot)
            v2dot=(dratio*w)**2*math.exp(-dratio*w*t)*(A*math.cos(wd*t)+B*math.sin(wd*t))-dratio*w\
                *math.exp(-dratio*w*t)*(-A*wd*math.sin(wd*t)+B*wd*math.cos(wd*t))-dratio*w*math.exp(-dratio*w*t)\
                *(-A*wd*math.sin(wd*t)+B*wd*math.cos(wd*t))+math.exp(-dratio*w*t)*(-A*wd**2*math.cos(wd*t)-B*wd**2*math.sin(wd*t))
            v2dotn.append(v2dot)

        disp=np.fabs(vn).max()*100
        vel=np.fabs(vdotn).max()*100
        accel=np.array(v2dotn)+np.array(acc)
        sa=np.fabs(accel).max()/9.81
        sv=vel
        sd=disp
        return T,sa,sv,sd


    def Ia (self):
        #Compound acc.-related IM  (g.s**1/3)
        #Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn30(12):1791–1816
        #Ia=PGA*(td**1/3) td=t2-t1   t1=t(5%AI)  t2=t(95%AI)
        ai=self.AI()
        thresholdt1=0.05*ai
        thresholdt2=0.95*ai
        t=self.t
        acc=self.acc*9.81
        acclist=[0]
        t1=[]
        t2=[]
        for i in range(len(acc)):
            acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) 
        for i1 in range(len(acclist)):
            if acclist[i1]>=thresholdt1:
                t1.append(i1*t)
                break
        for i2 in range(len(acclist)):
            if acclist[i2]>=thresholdt2:
                t2.append(i2*t)
                break
        td=t2[0]-t1[0]
        pga=self.PGA()
        Ia=pga*(td**(1.0/3.0))
        return td, Ia


    def FI (self):
        #Fajfar intensity  (cm/s**0.75)
        #Fajfar P, Vidic T, Fischinger M (1990) A measure of earthquake motion capacity to damage medium-period
        #structures. Soil Dyn Earthq Eng 9(5):236–242
        #FI=PGV*td**0.25   td=t2-t1   t1=t(5%AI)  t2=t(95%AI)
        ai=self.AI()
        thresholdt1=0.05*ai
        thresholdt2=0.95*ai
        t=self.t
        acc=self.acc*9.81
        acclist=[0]
        t1=[]
        t2=[]
        for i in range(len(acc)):
            acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) 
        for i1 in range(len(acclist)):
            if acclist[i1]>=thresholdt1:
                t1.append(i1*t)
                break
        for i2 in range(len(acclist)):
            if acclist[i2]>=thresholdt2:
                t2.append(i2*t)
                break
        td=t2[0]-t1[0]
        pgv=self.PGV()
        fi=pgv*(td**(0.25))
        return fi


    def Iv (self):
        #Compound vel.-related IM (cm**2/3)/s**1/3
        #Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn
        #30(12):1791–1816
        #IvPGV**(2/3)*td(1/3)
        ai=self.AI()
        thresholdt1=0.05*ai
        thresholdt2=0.95*ai
        t=self.t
        acc=self.acc*9.81
        acclist=[0]
        t1=[]
        t2=[]
        for i in range(len(acc)):
            acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) 
        for i1 in range(len(acclist)):
            if acclist[i1]>=thresholdt1:
                t1.append(i1*t)
                break
        for i2 in range(len(acclist)):
            if acclist[i2]>=thresholdt2:
                t2.append(i2*t)
                break
        td=t2[0]-t1[0]
        pgv=self.PGV()
        Iv=pgv**(2.0/3.0)*(td**(1.0/3.0))
        return Iv


    def Id (self):
        #Compound disp.-related IM (cm.s**1/3)
        ##Riddell R, Garcia J (2001) Hysteretic energy spectrum and damage control. Earthq Eng Struct Dyn
        #30(12):1791–1816
        #Id=PGD*td**1/3
        ai=self.AI()
        thresholdt1=0.05*ai
        thresholdt2=0.95*ai
        t=self.t
        acc=self.acc*9.81
        acclist=[0]
        t1=[]
        t2=[]
        for i in range(len(acc)):
            acclist.append((math.pi/(2*9.81))*(acc[i]**2*self.t)+acclist[-1]) 
        for i1 in range(len(acclist)):
            if acclist[i1]>=thresholdt1:
                t1.append(i1*t)
                break
        for i2 in range(len(acclist)):
            if acclist[i2]>=thresholdt2:
                t2.append(i2*t)
                break
        td=t2[0]-t1[0]
        pgd=self.PGD()
        Id=pgd*(td**(1.0/3.0))
        return Id

        

        

        

            

        


            
        
        
        
                

if __name__=='__main__':

    acc=np.loadtxt("FiltedAcceleration/1.out")

    imInstance=IMs(acc,0.02)
    pga=imInstance.PGA()
    print ("PGA=",pga)
    vel=imInstance.AcctoVelocity()
    pgv=imInstance.PGV()
    print ("PGV=",pgv)
    disp=imInstance.AcctoDisp()
    pgd=imInstance.PGD()
    print ("PGD=",pgd)
#   ratio=imInstance.VmaxDivAmax()
#   print (ratio)
#   arem=imInstance.aRMS()
#   print (arem)
#   vrms=imInstance.vRMS()
#   print (vrms)
#   drms=imInstance.dRMS()
#   print (drms)
#   ai=imInstance.AI()
#   print (ai)
#   ic=imInstance.Ic()
#   print (ic)
#   sed=imInstance.SED()
#   print (sed)
#   cav=imInstance.CAV()
#   print (cav)
#   asi=imInstance.ASI()
#   print (asi)
#   vsi=imInstance.VSI()
#   print (vsi)
#   sma=imInstance.SMA()
#   print (sma)
#   smv=imInstance.SMV()
#   print (smv)
#   fft=imInstance.lowCutFilter(5,3)
#   fft=imInstance.highCutFilter(5,9)
#   print (np.fabs(fft).max())
#   a95=imInstance.A95 ()
#   print (a95)
#   Tpredom=imInstance.Tpre()
#   print (Tpredom)
#   tm=imInstance.Tm()
#   print (tm)
#   eda=imInstance.EDA()
#   print (eda)
#   tduration=imInstance.DurationBrac()
#   print tduration
#   duhamel=imInstance.SpPiecewise(1,1,0.05)
#   print (duhamel)
#   ia=imInstance.Ia()
#   print (ia)
#   FI=imInstance.FI()
#   print (FI)
#   iv=imInstance.Iv()
#   print (iv)
#   Id=imInstance.Id()
#   print (Id)
##  a=[]
#   for i3 in np.arange(0.02,4.02,0.02):
#       result=imInstance.SpPiecewise (1,i3,0.05)
#       a.append(result[1])
#   np.savetxt('aaa.txt',np.mat(a).T,fmt="%f")

    



        



    