
%% This code generates Fig 2 a) and b) of the main text (experimental data)

clear variables;
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
Tlist=[0,1,2,3,4,5];

TrRho2=zeros(length(Tlist),N);
TrRho2_std=zeros(size(TrRho2));

for t=1:length(Tlist)
    qstates=csvread(strcat('Data_aau4963/10Ions_CleanSystem/MeasuredStates_T_',int2str(Tlist(t)),'ms.csv'));
    [pur,std]=ExtractPurity_Direct(qstates,A_Subs,N,CalcErrors); 
    
    for l=1:N
        TrRho2(t,l)=pur{l}(1);
        TrRho2_std(t,l)=std{l}(1);
    end
    
end


figure(20000)
clf
hold on;
for t=1:length(Tlist)
        errorbar(1:N,TrRho2(t,:),TrRho2_std(t,:),'o')
end
ax = gca;
ax.ColorOrderIndex = 1;
Tlist={'0','1','2','3','4','5'};
for t=1:length(Tlist)
   data=dlmread(strcat('Data_aau4963/10Ions_CleanSystem/Purity_T_',Tlist{t},'ms.csv'),'\t',1,0);
   plot(1:N,data(:,5))
end
title('Purity of connected subsystems including ion 1')
xlabel('i')
ylabel('Tr[\rho_{[1\rightarrow i]}^2]')
Renyi2=zeros(size(TrRho2));
Renyi2_std=zeros(size(TrRho2));

figure(20001)
clf
hold on;
for t=1:6
        Renyi2(t,:) = -log2(TrRho2(t,:));
        Renyi2_std(t,:) = TrRho2_std(t,:)./(TrRho2(t,:)*log(2));
        errorbar(1:N,Renyi2(t,:),Renyi2_std(t,:),'o')
end
Tlist={'0','1','2','3','4','5'};
ax = gca;
ax.ColorOrderIndex = 1;
for t=1:length(Tlist)
   data=dlmread(strcat('Data_aau4963/10Ions_CleanSystem/RenyiEntropy_T_',Tlist{t},'ms.csv'),'\t',1,0);
   plot(1:N,data(:,5))
end
title('Renyi entropy (connected subsystems including ion 1)')
xlabel('i')
ylabel('S^{(2)}(\rho_{[1\rightarrow i]})')
legend('0ms', '1ms', '2ms', '3ms', '4ms', '5ms')

