function [A_Subs,Subsets]=calcExpecMatrix_Subs(N,type)

    filename=strcat('Subroutines/ExpecMatrix_',int2str(N),'.mat');
    if not(exist(filename,'file')==2)
        A=calcExpecMatrix(N);  % Mapping basis of l qubits to all subsystems
        save(filename,'A');
        disp('Step 1a: ExpecMatrix calculated')
    else
        load(filename);
        disp('Step 1a: ExpecMatrix loaded')
    end
    
    filename=strcat('Subroutines/BasisMap_',int2str(N),'_',int2str(type),'.mat');
    if not(exist(filename,'file')==2)
        [BasisMap,Subsets]=TransformBasis_type(N,type);  % Mapping basis of l qubits to all subsystems
        save(filename,'BasisMap','Subsets','type');
        disp('Step 1b: Basis Map calculated')
    else
        load(filename);
        disp('Step 1b: Basis Map loaded')
    end
    clear('filename')
    
        number_of_subsystem_sizes=length(Subsets);
    
        A_Subs=cell(number_of_subsystem_sizes,2);
    
        for k=1:number_of_subsystem_sizes
            
            l=length(Subsets{k}(1,:));

            [number_of_subsets,~]=size(Subsets{k});
            Al=int16(A{l});

            A_Subs{k,1}=zeros(number_of_subsets,2^N,2^N,'int16');

            for s=1: number_of_subsets
               basis=BasisMap{k}(s,:);
               A_Subs{k,1}(s,:,:)=Al(basis,basis);
            end
            
            A_Subs{k,2}=l;

        end
        
     disp('Step 1c: ExpecMatrixSubs calculated')
