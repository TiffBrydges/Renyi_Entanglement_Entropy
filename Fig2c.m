%% This code generates Fig 2 c)  of the main text (experimental data)

clear variables;
addpath('Subroutines')

type=0;  % Type: 0: all subsystems, 1: all connected, 2: only connected and located left ( [1],[1,2],[1,2,3],... )
    
N=10;  % Number of ions

% Generate (or load) the matrix of coefficients    
filename=strcat('Subroutines/ExpecMatrix_Subs_',int2str(N),'_type_',int2str(type),'.mat');
if not(exist(filename,'file')==2)
        A_Subs=calcExpecMatrix_Subs(N,type);  % Mapping basis of l qubits to all subsystems
        save(filename,'A_Subs');
        % This file is ~200 MB! Storing  it will however speed up the
        % calculation.
else
        load(filename);
        disp('Step 1: ExpecMatrixSubs loaded')
end

% Calculate the Purities from the experimental data
CalcErrors=true;

qstates=csvread(strcat('Data_aau4963/10Ions_CleanSystem/MeasuredStates_T_',int2str(5),'ms.csv'));
[pur,std]=ExtractPurity_Direct(qstates,A_Subs,N,CalcErrors); 
    


figure(20000)
clf
hold on;
for l=1:N
        if length(pur{l})==1
            xvalues=l;
        else
            xvalues=l+linspace(-0.3,0.3,length(pur{l}));
        end
        errorbar(xvalues,pur{l},std{l},'o')
end
title(' Purity of all subsystems ')
xlabel('|A|')
ylabel('Tr[\rho_{A}^2]')




figure(20001)
clf
hold on;
for l=1:N
        if length(pur{l})==1
            xvalues=l;
        else
            xvalues=l+linspace(-0.3,0.3,length(pur{l}));
        end
        yvalues=-log2(pur{l});
        yerr=std{l}./(pur{l}*log(2));
        errorbar(xvalues,yvalues,yerr,'o')
end
title('Renyi entropy of all subsystems')
xlabel('|A|')
ylabel('S^{(2)}(\rho_A)')