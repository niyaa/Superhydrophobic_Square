
# # Stability Driver
# Runs the stability calculation with varing parameters

# In[ ]:

import numpy as np
import sys;
sys.path.append('~/');
import driveClass 
        


# In[163]:
c=driveClass.Case();
#alpha 
alStart = 1;  alEnd = 2; alStep = 0.02;

#S
sStart=0.45; 
#sigmas
siStart= 0; srStart= 0;
#Reynolds
ReStart = 75; ReStep= 5;
#betas
lzStart= 17.4;

clist=[];

#ReStart, lzStart, siStart, srStart=unstableLz(ReStart, lzStart, siStart, srStart, sStart, alStart, clist, 200);





driveClass.alphaContin(alStart, alEnd, alStep, sStart, ReStart, lzStart, siStart, srStart, clist, 40 , ReStep);











