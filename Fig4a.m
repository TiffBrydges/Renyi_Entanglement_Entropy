%% This code generates Fig 4 a)  of the main text (experimental data)


addpath('Subroutines')

type=2;  % Type: 0: all subsystems, 1: all connected, 2: only connected and located left ( [1],[1,2],[1,2,3],... )
    
N=10;  % Number of ions

% Generate (or load) the matrix of coefficients    
filename=strcat('Subroutines/ExpecMatrix_Subs_',int2str(N),'_type_',int2str(type),'.mat');
if not(exist(filename,'file')==2)
        A_Subs=calcExpecMatrix_Subs(N,type);  % Mapping basis of l qubits to all subsystems
        save(filename,'A_Subs');
else
        load(filename);
        disp('Step 1: ExpecMatrixSubs loaded')
end

% Calculate the Purities from the experimental data
CalcErrors=true;
Tlist_clean={'0','1','2','3','4','5'};
TrRho2_clean=zeros(length(Tlist_clean),N);
TrRho2_std_clean=zeros(size(TrRho2_clean));

for t=1:length(Tlist_clean)
    qstates=csvread(strcat('Data_aau4963/10Ions_CleanSystem/MeasuredStates_T_',Tlist_clean{t},'ms.csv'));
    [pur,std]=ExtractPurity_Direct(qstates,A_Subs,N,CalcErrors); 
    
    for l=1:N
        TrRho2_clean(t,l)=pur{l}(1);
        TrRho2_std_clean(t,l)=std{l}(1);
    end
end


Tlist_dis={'01','02','04','06','10','16','20'};
TrRho2_dis=zeros(length(Tlist_dis),N);
TrRho2_std_dis=zeros(size(TrRho2_dis));

for t=1:length(Tlist_dis)
    qstates=csvread(strcat('Data_aau4963/10Ions_withDisorder/MeasuredStates_T_',Tlist_dis{t},'ms.csv'));
    [pur,std]=ExtractPurity_Direct(qstates,A_Subs,N,CalcErrors); 
    
    for l=1:N
        TrRho2_dis(t,l)=pur{l}(1);
        TrRho2_std_dis(t,l)=std{l}(1);
    end
end


figure(20000)
clf
hold on;
errorbar(cellfun(@str2num,Tlist_clean),TrRho2_clean(:,5),TrRho2_std_clean(:,5),'o')
errorbar(cellfun(@str2num,Tlist_dis),TrRho2_dis(:,5),TrRho2_std_dis(:,5),'o')

data=dlmread(strcat('Data_aau4963/10Ions_withDisorder/RenyiTimeEvo_HalfPartition_Numerics.csv'),'\t',1,0);
ax = gca;
ax.ColorOrderIndex = 1;
plot(data(:,1),power(2,-data(:,2)),'--')
ax.ColorOrderIndex = ax.ColorOrderIndex-1;
plot(data(:,1),power(2,-data(:,3)),'-')
plot(data(:,1),power(2,-data(:,4)),'--')
ax.ColorOrderIndex = ax.ColorOrderIndex-1;
plot(data(:,1),power(2,-data(:,5)),'-')
title('Time evolution of purity of half partition')
xlabel('t[ms]')
ylabel('Tr[\rho_{[1\rightarrow 5]}^2]')


figure(20001)
clf
hold on;
errorbar(cellfun(@str2num,Tlist_clean),-log2(TrRho2_clean(:,5)),TrRho2_std_clean(:,5)./(TrRho2_clean(:,5)*log(2)),'o')
errorbar(cellfun(@str2num,Tlist_dis),-log2(TrRho2_dis(:,5)),TrRho2_std_dis(:,5)./(TrRho2_dis(:,5)*log(2)),'o')


data=dlmread(strcat('Data_aau4963/10Ions_withDisorder/RenyiTimeEvo_HalfPartition_Numerics.csv'),'\t',1,0);
ax = gca;
ax.ColorOrderIndex = 1;
plot(data(:,1),data(:,2),'--')
ax.ColorOrderIndex = ax.ColorOrderIndex-1;
plot(data(:,1),data(:,3),'-')
plot(data(:,1),data(:,4),'--')
ax.ColorOrderIndex = ax.ColorOrderIndex-1;
plot(data(:,1),data(:,5),'-')
title('Time evolution of Renyi entropy of half partition')
xlabel('t[ms]')
ylabel('S^{(2)}(\rho_{[1\rightarrow 5]})')


