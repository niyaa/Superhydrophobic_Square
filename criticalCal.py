from glob import glob
from scipy import interpolate
from scipy import  optimize
import numpy as np
import os

def criticalReBeta3(d,prnt=0):

    y=d[:,4];
    Re=list(set(y));
    Re.sort();
    betaCr=[];
    sigmaMax=[];
    crsigmaR=[];
    ReNew=[];
    for i in Re:
        xx=d[d[:,4]==i][:,2];
        yy=d[d[:,4]==i][:,5];
        zz=d[d[:,4]==i][:,6];
        l1, l2 = zip(*sorted(zip(xx, yy)))
        xx=list(l1);yy=list(l2);
        if(prnt==1):print(i);
        xx=np.asarray(xx);yy=np.asarray(yy);
        if len(xx) >=3:
            if(prnt==1):print(xx); print(yy);
            try:
                f=interpolate.interp1d(xx,-yy,kind='quadratic',fill_value="extrapolate");
                f1=interpolate.interp1d(xx,zz,kind='quadratic',fill_value="extrapolate")
            except:
                f=interpolate.interp1d(xx,-yy,kind='quadratic');
                f1=interpolate.interp1d(xx,zz,kind='quadratic')
        else:
            if(prnt==1):print(xx); print(yy);
            f=interpolate.interp1d(xx,-yy,kind='linear',fill_value="extrapolate");
            f1=interpolate.interp1d(xx,zz,kind='linear',fill_value="extrapolate")
        try:
            betaTemp=optimize.fmin(f,round(np.average(xx),1),maxiter=100000);
            betaCr.append(betaTemp[0]);
            sigmaMax.append(-f(betaTemp).flatten()[0]);
            ReNew.append(i);
            crsigmaR.append(f1(betaTemp).flatten()[0]);

        except:
            if(prnt==1):print(str(i)+" No \n");
            continue;

    l3, l2, l1, l0 = zip(*sorted(zip(ReNew,sigmaMax,betaCr, crsigmaR)));
    sigmaMax=list(l2); betaCr=list(l1);ReNew=list(l3); crsigmaR=list(l0);
    if(prnt==1):print('The most important 1. RE 2. sigmaMax 3. betacr 4. crsigmaR \n');
    if(prnt==1):print(ReNew); print(sigmaMax); print(betaCr); print(crsigmaR);
    Re=ReNew;

    f=interpolate.interp1d(Re,sigmaMax,kind='linear');
    f1=interpolate.interp1d(Re,betaCr,kind='linear');
    f2=interpolate.interp1d(Re,crsigmaR,kind='linear');

    s=np.where(np.diff(np.sign(sigmaMax)));
    if(prnt==1):print("the sign change "+str(s))
    recr=optimize.brentq(f,Re[s[0][0]],Re[s[0][0]+1]);
    betacr=f1(recr);
    crsigmaR=f2(recr)
    if(prnt==1):print(Re);print('\n'); print(sigmaMax); print('\n'); print(betaCr);
    print("The Converged values of the critical Re and Beta are " +str(recr) +" ------------ "+str(betacr));
    return recr, betacr.tolist(), crsigmaR.tolist();


path = glob("*/");
print(path)
path.remove("__pycache__/")

ff= []; ReCr = []; BetaCr = []; SigmarCr = [];
cwd = os.getcwd();
for i in path:
    os.chdir(cwd);
    print(i)
    os.chdir(i);
    path1 = glob("*.dat");
    a = np.loadtxt(path1[0], delimiter =",");
    try:
        print(i)
        [recr, betacr, sigmaR] = criticalReBeta3(a, 1);
        ff.append(np.float(i.split("/")[0]));
        ReCr.append(recr);
        BetaCr.append(betacr)
        SigmarCr.append(sigmaR)
    except:
        print(i)
        print(os.getcwd());
        continue

os.chdir(cwd);
print(ff)
print(ReCr);
print(BetaCr);
print(SigmarCr);
c =np.array(SigmarCr)/np.array(BetaCr)
print(c.tolist())
c = np.asarray((ff, ReCr, BetaCr, c)).T
c = c[c[:,0].argsort()]
np.savetxt("data.txt", c);
