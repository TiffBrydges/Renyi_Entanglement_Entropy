
%% This code generates Fig 3 of the main text (experimental data)


addpath('Subroutines')
clear variables

type=2;  % Type: 0: all subsystems, 1: all connected, 2: only connected and located left ( [1],[1,2],[1,2,3],... )
    
N=10;  % Number of ions (measured)

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
Tlist={'00','01','02','03','04','06','10'};

TrRho2=zeros(length(Tlist),N);
TrRho2_std=zeros(size(TrRho2));
TrRho2_num=zeros(size(TrRho2));

for t=1:length(Tlist)
    qstates=csvread(strcat('Data_aau4963/20Ions_CleanSystem/MeasuredStates_T_',Tlist{t},'ms.csv'));
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
        errorbar(1+5:N+5,TrRho2(t,:),TrRho2_std(t,:),'o')
end
title('Purity of connected subsystems including ion 6')
xlabel('i')
ylabel('Tr[\rho_{[6\rightarrow i]}^2]')

xlim([5,16])

figure(20001)
clf
hold on;
for t=1:5
        errorbar(1+5:N+5, -log2(TrRho2(t,:)),TrRho2_std(t,:)./(TrRho2(t,:)*log(2)),'o')
end
errorbar(1+5:8+5, -log2(TrRho2(6,1:8)),TrRho2_std(6,1:8)./(TrRho2(6,1:8)*log(2)),'o')
errorbar(1+5:7+5, -log2(TrRho2(7,1:7)),TrRho2_std(7,1:7)./(TrRho2(7,1:7)*log(2)),'o')
ax = gca;
ax.ColorOrderIndex = 1;
Tlist={'0','1','2','3','4','6','10'};
for t=1:length(Tlist)
   data=dlmread(strcat('Data_aau4963/20Ions_CleanSystem/RenyiEntropy_T_',Tlist{t},'ms.csv'),'\t',1,0);
   plot(1+5:N+5,data(:,5))
end

xlim([5,16])

title('Renyi entropy (connected subsystems including ion 6)')
xlabel('i')
ylabel('S^{(2)}(\rho_{[6\rightarrow i]})')
legend('0ms', '1ms', '2ms', '3ms', '4ms', '6ms', '10ms', 'Location', 'northwest')