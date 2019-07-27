function Project0904(nn,p)

% clc; clear all;close all;
%Network Deffuant and network MajorityRule By random network with <k>=5, N=10^4
N=10^3;
avgK=5;
% p=0.2; %Interaction Strenth between two network layer
Epsilon=0.1; %For Deffuant model, the difference between one pair *
mu=0.5; %for Deffuant model, convergence parameter
timestep=1e5;
jump=timestep/50;  %for the record of opinion information
InitalMNegative=0.1; % initial percent of Negative opinion


fid0=fopen('trace.txt','wt+');
fid1=fopen('DeffuantOpinionLayer.txt','wt+');
fid2=fopen('MajorityRuleOpinionLayer.txt','wt+');
    
NetD=zeros(N,N);
NetMR=zeros(N,N);
%random networkD and random networkMR
tic
for i=1:N
    i;
     %fprintf(fid0,  ' \n node i   %g\n' ,i );  
     % fprintf(fid0,  'connect to  \t  ' );  
     countD=0;
    while ( sum(NetD(i,:))<avgK && countD<500*N) 
       ConnetN=randi(N,1);
       countD=countD+1;
      %  fprintf(fid0,' chose %1.0f    ', ConnetN );  
       if( NetD(i,ConnetN)==0 && sum(NetD(ConnetN,:))<avgK && i~=ConnetN )
           NetD(i,ConnetN)=1;
           NetD(ConnetN,i)=1;    
       end   
    end
    
    if (countD==500*N)
        break;
    end
end


 if (countD==500*N)
     texte=['p=' num2str(p) '   sampling' num2str(nn) ' , random networkD failed'];
     fprintf(fid0,[texte '\n']);
     disp(texte)
        return
 else
     texte=['NetworkD finished at ' num2str(toc/60) ' mints.'];
disp(texte)
 end
    
%fprintf(' random networkD \n');
tic
for i=1:N
    i;
    countMR=0;
    while ( sum(NetMR(i,:))<avgK && countMR<500*N )
       ConnetN=randi(N,1);
       countMR=countMR+1;
       if( NetMR(i,ConnetN)==0 && sum(NetMR(ConnetN,:))<avgK && i~=ConnetN )
           NetMR(i,ConnetN)=1;
           NetMR(ConnetN,i)=1;    
       end   
    end
     if (countMR==500*N)
        break;
    end
end


 if (countMR==500*N)
         texte=['p=' num2str(p) '   sampling' num2str(nn) ' , random networkMR failed'];
          fprintf(fid0,[texte '\n']);
          disp(texte)
        return
 else
     texte=['NetworkMR finished at ' num2str(toc/60) ' mints.'];
disp(texte)
 end

%fprintf(' random Network MR \n');


%layer connection between two random networks, determinted by p. With
%probability p, each node in one layer can estabilish one link to the other layer.
%LayerConnet is connection between two layers.
RandNet=rand(1,N);
ind=[];
ind=find(RandNet<p);
LayerConnet=zeros(1,N);
LayerConnet(ind)=1;

fprintf(' LayerConnet \n '); 


NetDO=zeros(1,N);
NetMRO=ones(1,N);

%set the initial opinion state of each node
NetDO=-1+2*rand(1,N); % Uniformly distributed pseudorandom number in opinion[-1,1]  for Deffuant Network
% NetMRO1=randi([0 1],1,N);% Uniformly distributed pseudorandom number in opinion[-1,1]  for MajorityRule Network
% ZeroSet=find(NetMRO1==0);
% NetMRO1(ZeroSet)=-1;
% NetMRO=NetMRO1;

InitaNegN=randi(N,N*InitalMNegative,1)' ; % initial percent of Negative opinion
NetMRO(InitaNegN)=-1;


%%the dynamics rule: each layer evloves at same time.
for t=1:timestep

    %%Deffuant layer
    choseN=randi(N);
    NeigNSameLayer=[];
    NeigNDiffLayer=[];
    NeigNSameLayer=find(  NetD(choseN,:)~=0 );
    if ( LayerConnet(choseN)==1 )
        NeigNDiffLayer=choseN;
    end
    NeighN=[ NeigNSameLayer, -NeigNDiffLayer ];
    NNeighSet=[];
    NNeighSet=randperm( length(NeighN) );
    NNeigh=NNeighSet(1);
    %%neighbor node in the same layer
    if (  NNeigh>0  ) %%  
        if ( abs(NetDO(choseN)-NetDO(NNeigh))<Epsilon )
            NetDO(choseN)=NetDO(choseN)+mu*( NetDO(NNeigh)-NetDO(choseN) );
            NetDO(NNeigh)=NetDO(NNeigh)+mu*( NetDO(choseN)-NetDO(NNeigh) );
        end
    end 
    %%neighbor node in the different layer
    if ( NNeigh<0   ) %%
        if (  abs(NetDO(choseN)-NetMRO(-NNeigh))<Epsilon )
             NetDO(choseN)=NetDO(choseN)+mu*( NetMRO(-NNeigh)-NetDO(choseN) );
             NetMRO(-NNeigh)=NetMRO(-NNeigh)+mu*(NetDO(choseN)-NetMRO(-NNeigh));
        end
    end
    
    
     %%Majority Rule layer
    choseN=randi(N);
    NeigNSameLayer=[];
    NeigNDiffLayer=[];
    NeigNSameLayer=find(  NetMR(choseN,:)~=0 );
    if ( LayerConnet(choseN)==1 )
        NeigNDiffLayer=choseN;
    end
    
    
    %%neighbor node in the same layer
    nei=sum( NetMRO(NeigNSameLayer) )+ sum( NetD(NeigNDiffLayer) ) ;
    if ( nei>0 )  
        NetMRO(choseN)=1;
    end
    if (nei<0 )
         NetMRO(choseN)=-1;
    end 
    
    %%output for opinion state
   if ( rem(t,jump)==0 )
       t;
    fprintf(fid1,'  %12.8f  ', NetDO  );
    fprintf(fid1,' \n ' );
    fprintf(fid2,'  %g  ', NetMRO );
    fprintf(fid2,' \n ' );
   end
   
end %%end of evolution one timestep

 fclose(fid0);
 fclose(fid1);
 fclose(fid2);

DeffOpin=load(['DeffuantOpinionLayer.txt']);
[x,y]=size(DeffOpin);
figure(1)
plot(1:x,DeffOpin(:,:) );
xlabel([ 'Evolution Timestep=x*jump(' , num2str(jump), ')'] );
ylabel('Deffuant Layer Opinion');
title(['p=' num2str(p) 'sampling' num2str(nn) 'DeffuantEvolution']);
set(gca, 'YLim',[-1 1]);   
str1=['p=' num2str(p) 'sampling' num2str(nn) 'DeffuantEvolution.jpg'];
saveas(gcf,str1);

MROpin=load(['MajorityRuleOpinionLayer.txt']);
[x,y]=size(MROpin);
figure(2)
trans=MROpin';
plot(1:x,sum( trans(:,:))/N);
xlabel([ 'Evolution Timestep=x*jump(' , num2str(jump), ')'] );
ylabel('MR Layer Opinion');
title(['p=' num2str(p) 'sampling' num2str(nn) 'MRLayerEvolution']);
set(gca, 'YLim',[-1.5 1.5]);
str2=['p=' num2str(p) 'sampling' num2str(nn) 'MRLayerEvolution.jpg'];
saveas(gcf,str2);
