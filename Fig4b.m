%% This code generates Fig 4 b)  of the main text (experimental data)


addpath('Subroutines')
clear variables;

type=5;  % Type: 0: all subsystems, 1: all connected, 2: only connected and located left ( [1],[1,2],[1,2,3],... )
        %5: Subsystems to calculate mutual information as shown in Fig 4b
        
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



Tlist_dis={'01','02','04','06','10','16','20'};
MutualInfo=zeros(length(Tlist_dis),N);
MutualInfo_std=zeros(size(MutualInfo));


for t=1:length(Tlist_dis)
    qstates=csvread(strcat('Data_aau4963/10Ions_withDisorder/MeasuredStates_T_',Tlist_dis{t},'ms.csv'));
    [NU,NM]=size(qstates);
    [pur,~]=ExtractPurity_Direct_4b(qstates,A_Subs); 
    
    for l=1:3
        
        
        MutualInfo(t,l)=-log2(mean(pur{1}(2*l-1,:)))-log2(mean(pur{1}(2*l,:)))+log2(mean(pur{2}(l,:)));
      
        % Do Jacknife resampling to estimate the standard error
        ls=1:NU;
        jacknifesamples=zeros(NU,1);
   
        for r=1:NU
            jacknifesamples(r,1)=-log2(mean(pur{1}(2*l-1,ls~=r)))-log2(mean(pur{1}(2*l,ls~=r)))+log2(mean(pur{2}(l,ls~=r)));
        end
                
        MutualInfo_std(t,l)=std(jacknifesamples,1)*sqrt(NU-1);
    end
end

data_u=dlmread(strcat('Data_aau4963/10Ions_withDisorder/RenyiMutInfo_Numerics_Unitary.csv'),'\t',1,0);
data_d=dlmread(strcat('Data_aau4963/10Ions_withDisorder/RenyiMutInfo_Numerics_inclDec.csv'),'\t',1,0);



figure(20001)
clf
ax = gca;
hold on;
for l=1:3
    errorbar(cellfun(@str2num,Tlist_dis),MutualInfo(:,l),MutualInfo_std(:,l),'o')
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    plot(data_u(:,1),data_u(:,l+1),'--');
    ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
    plot(data_d(:,1),data_d(:,l+1),'-');
end


title('Mutual information')
xlabel('t[ms]')
ylabel('I^{(2)}(\rho_A : \rho_B)')

