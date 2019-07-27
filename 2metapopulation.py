####uncorrelated configuration model(UCM)

print(" I could have this ")
import math
import random
import numpy as np
import matplotlib.pyplot as plt
import pylab as pl
from collections import Counter


NetNode=1000
NetGama=3
Mmin=2
Mmax=math.floor(math.sqrt(NetNode))

Rbeta=0.1
Rmu=0.1

Dtheta=0
Dalpha=0
Dgama=0



Net=list(range(1,NetNode+1))

####Power-law degree dis
k=[]
while len(k) < NetNode:
        r=random.random()
        powerRandom=0.50+(Mmin-0.50)*np.power(1-r,-1/(NetGama-1) )
        if math.floor(powerRandom<math.sqrt(NetNode)):
                degree=math.floor(powerRandom)
                k.append(float(degree))

if sum(k)%2 !=0:
        k[NetNode-1]=k[NetNode-1]+1

    
####write degreefile
Degfile=open("Degree.txt", "w")
Degfile.write(str(k))
Degfile.close()

NetPara=open("Netparameter.txt", "w")
NetPara.write("Average degree of Network is : ")
NetPara.write( str( np.mean(k) ) )
NetPara.write("\n")
powerk=np.power(k,2)
NetPara.write("Average powder 2 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")
powerk=np.power(k,4)
NetPara.write("Average powder 4 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")
powerk=np.power(k,-1)
NetPara.write("Average powder -1 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

powerk=np.power(k,-2)
NetPara.write("Average powder -2 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

powerk=np.power(k,-4)
NetPara.write("Average powder -4 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

powerk=np.power(k,0.5)
NetPara.write("Average powder 1/2 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

powerk=np.power(k,4)
NetPara.write("Average powder 3 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

powerk=np.power(k,3/2)
NetPara.write("Average powder 3/2 degree of Network is : ")
NetPara.write( str( np.mean(powerk)) )
NetPara.write("\n")

NetPara.close()

data=[]
ou=open("Degree.txt","r+")
for line in ou:
    data.append(line)

####get degree distribution plot
DegPro=Counter(k)
s=dict(DegPro)
x=list(s.keys())
y=list(s.values())

plt.loglog( x, y, marker='o', linestyle='--', color='r')
plt.xlabel('Degree')
plt.ylabel('Degree Distribution')
####plt.show()
plt.savefig("Powerlaw-DegreeDis.png")
plt.close()

####construct UCM
ConnectionMatric=np.zeros((NetNode,NetNode))
Element=[]
for i in range(len(k)):
        for j in range(int(k[i])):
                Element.append(i+1)
        
print("sum K is ", sum(k))
print("sum in the element is ", len(Element))
 

while len(Element)>0:
        index1=random.randint(0,len(Element)-1)
        node1=Element[ index1  ]
        TempleElem=[]
        TempleElem=Element[:]
        TempleElem.pop( index1 )

        index2=random.randint(0,len(TempleElem)-1)
        node2=TempleElem[ index2 ]
        if ConnectionMatric[  node1-1, node2-1 ]==0  and node1!=node2 :
                ConnectionMatric[  node1-1, node2-1 ]=1
                ConnectionMatric[  node2-1, node1-1 ]=1
                TempleElem.pop( index2 )
                Element=[]
                Element=TempleElem[:]        

np.savetxt('NetWork.txt',ConnectionMatric,fmt='%-7.0f')
print("Net UCM is ready !!")

for node in Net:
        if sum(ConnectionMatric[node-1,:])!=k[node-1]:
                print("errorr!!!!")



######initial conditions setting
DensityRho1=np.arange(50,100,50)
DensityRhoSet=DensityRho1/100

for DensityRho in DensityRhoSet:
        Individuals=int( DensityRho*NetNode )
        IniS=random.randint(1,Individuals-1)
        IniI=Individuals-IniS
        IndiMatric=np.zeros((NetNode,2))
        ChoseNode=np.random.permutation(NetNode)
        Sum=IndiMatric.sum(axis=0)
        LeftS=IniS-Sum[0]
        LeftI=IniI-Sum[1]
        ct=0

        while LeftS>0 or LeftI>0:
                s=random.randint(0,LeftS)
                i=random.randint(0,LeftI)
                IndiMatric[ChoseNode[ct]][0]=s
                IndiMatric[ChoseNode[ct]][1]=i
                Sum=IndiMatric.sum(axis=0)
                LeftS=max( IniS-Sum[0],0 )
                LeftI=max( IniI-Sum[1], 0 )
                ct=ct+1

        print("Density Rho is ",DensityRho, "with total individuals is ",     Individuals )        
        print("initial S and initial I is ", sum( IndiMatric))


######partial diffusion coefficient


        PDiffusionMatric=np.zeros((NetNode,NetNode))
        NeigNode=[[]]*NetNode
        for node in Net:
                neig=ConnectionMatric[node-1,:]
                inds=np.nonzero(neig>0)
                neignode=[]
                neigk=[]
        
                for i in range(len(inds[0])):
                        neignode.append(inds[0][i])
                        neigk.append( k[inds[0][i]] )
                NeigNode[node-1]=neignode
                Neigk=np.asarray(neigk)
                tem=k[node-1]*Neigk
                kpow=np.power(tem,Dtheta)
                Pro=kpow/sum(kpow)
                for i in range(len(inds[0])):
                        PDiffusionMatric[node-1, inds[0][i]   ]=Pro[i]
 
        print("PDiffusionMatric is ready ")

######Evoluation process
        Sum=IndiMatric.sum(axis=0)
        RemainS=Sum[0]
        RemainI=Sum[1]
        step=0


        while RemainI>0 and step<100000:
                step=step+1
                print("step is ", step, "  and RemainI is ", RemainI)
####a1step. each infected individual become susceptible with probability mu
####a2 step. each susceptible individual become infrected with probability 1-(1-beta)^nl, nl is the total number of the infected individuals in the node

                for i in Net:
                        ReactionInfectRan=np.random.rand(IndiMatric[i-1][1])
                        ids=np.nonzero(ReactionInfectRan>Rmu)
                        InfectIntoSuspect=len(ids[0])
                        ReactionSuspectRan=np.random.rand(IndiMatric[i-1][0])
                        pro=1-pow((1-Rbeta), IndiMatric[i-1][1])
                        ids=np.nonzero(ReactionSuspectRan>pro)
                        SuspectIntoInfect=len(ids[0])
                        IndiMatric[i-1][1]=IndiMatric[i-1][1]-InfectIntoSuspect+SuspectIntoInfect
                        IndiMatric[i-1][0]=IndiMatric[i-1][0]+InfectIntoSuspect-SuspectIntoInfect

                        Sum=IndiMatric.sum(axis=0)
                        RemainS=Sum[0]
                        RemainI=Sum[1]

####b1 Diffusion, each individual moves into a neighbor with probability A(I) and A(S)
                        DInfectexponent=math.exp(-Dalpha*RemainI/Individuals)
                        DSuspectexponent=math.exp(-Dgama*RemainI/Individuals)
                        for i in Net:
                                if IndiMatric[i-1][1]!=0 :
                                        DiffusionInfectRan=np.random.rand(IndiMatric[i-1][1])
                                        repe=divmod( IndiMatric[i-1][1] , len(NeigNode[i-1])  )
                                        NeigList=NeigNode[i-1]*(repe[0]+1)
                                        NeigSet=NeigList[ 0: int(IndiMatric[i-1][1] ) ]
                                        NeigPro=[]
                                        for j in NeigSet:
                                                diffusionpro=PDiffusionMatric[ i-1, j  ]
                                                NeigPro.append(diffusionpro*DInfectexponent)
                                        DiffusionInfect=np.asarray(NeigPro)
                                        Dif=DiffusionInfect-DiffusionInfectRan
                                        ids=np.nonzero(Dif<0)
                                        if len(ids[0])!=0:
                        ###print("node",i, "transfer Infected", ids[0] )
                        ####print("IndiMatric node before diffusion",IndiMatric[i-1])
                                                if len(ids[0])>IndiMatric[i-1][1] :
                                                        break
                                                IndiMatric[i-1][1]=IndiMatric[i-1][1]-len(ids[0])
                                                for l in range(len(ids[0])):
                                                        IndiMatric[ ids[0][l] ][1]=IndiMatric[ ids[0][l] ][1]+1

                                if IndiMatric[i-1][0]!=0:
                                        DiffusionSuspectRan=np.random.rand(IndiMatric[i-1][0])
                                        repe=divmod( IndiMatric[i-1][0], len(NeigNode[i-1]) )
                                        NeigList=NeigNode[i-1]*(repe[0]+1)
                                        NeigSet=NeigList[  0:int(IndiMatric[i-1][0])  ]
                                        NeigPro=[]
                                        for j in NeigSet:
                                                diffusionpro=PDiffusionMatric[i-1,j]
                                                NeigPro.append(diffusionpro*DSuspectexponent)
                                        DiffusionSuspect=np.asarray(NeigPro)
                                        Dif=DiffusionSuspect-DiffusionSuspectRan
                                        ids=np.nonzero(Dif<0)
                                        if len(ids[0])!=0:
                                ####print("node",i, "transfer Suspected", ids[0] )
                                        ####print("IndiMatric node before diffusion",IndiMatric[i-1])
                                                if len(ids[0])>IndiMatric[i-1][0]:
                                                        break
                                                IndiMatric[i-1][0]=IndiMatric[i-1][0]-len(ids[0])
                                                for l in range(len(ids[0])):
                                                        IndiMatric[ ids[0][l] ][0]=IndiMatric[ ids[0][l] ][0]+1

    
 ###               print("after diffusion the IndiMatric",sum(IndiMatric))
                
             
        ########################write result file####################################
        filename1="1.infectedRation with Dtheta"+str(Dtheta)+" Dalpha "+str(Dalpha)+"Dgama "+str(Dgama)+".txt"
        InfectedRatiofile=open(filename1, "a+")
        result1=[]
        result1.append(DensityRho)
        result1.append(RemainI)
        InfectedRatiofile.write(str(result1) )
        InfectedRatiofile.write("\n")

 

 

