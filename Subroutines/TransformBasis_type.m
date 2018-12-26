function [BasisMap_ret,Subsets_ret]=TransformBasis_type(N,type)
% Returns the BasisMap (a cell array of dimension Nx1) to map the basis of a system with N spins
% to all its possible subsystems with number of spins l=1...N
% BasisMap{l}(s,i) maps the i-th basis state of the global system (N spins)
% to the basis state of the subsystem s consisting of l spins.

% Subsets is a cell array of dimension Nx1, Subsets{l} contains all
% subsystems consisting of l spins.

% type
% 0: all subsystems
% 1: all connected subsystems
% 2: all connected subsystems containing ion 1, [1],[1,2],...
% 3: subsystems of the form [N/2],[N/2+1], [N/2,N/2+1],  [N/2-1],[N/2+2],[N/2-1,N/2+2] ,..
% 4: subsystems of the form [N/2,N/2-1],[N/2+1,N/2+2],
% [N/2,N/2-1,N/2+1,N/2+2],  [N/2-1,N/2-2],[N/2+2,N/2+3],
% [N/2-1,N/2-2,N/2+2,N/2+3], ...
% 5: subsystems of the form [N/2,N/2-1,N/2-2],[N/2+1,N/2+2,N/2+3],
% [N/2,N/2-1,N/2-2,N/2+1,N/2+2,N/2+3],  [N/2-1,N/2-2,N/2-3],[N/2+2,N/2+3,N/2+4],
% [N/2-1,N/2-2,N/2-3,N/2+2,N/2+3,N/2+4]

    Nh=2^N;

    I=cell(Nh,1);
    for i=1:Nh
        I{i}=dec2bin(i-1,N);
    end

    Subsets=cell(N,1);
    BasisMap=cell(N,1);
    for l=1:N
       Subsets{l}=nchoosek(1:N,l);
       BasisMap{l}=zeros(nchoosek(N,l),Nh);
       for s=1:nchoosek(N,l)             
            for i=1:Nh
                    BasisMap{l}(s,i)=bin2dec(I{i}(Subsets{l}(s,:)))+1;
            end
        end
    end
    
    if type ==1
        
        SubsetsConnected=cell(N,1);
        BasisMapConnected=cell(N,1);
        for l=1:N
            SubsetsConnected{l,1}=zeros(N-(l-1),l);
            BasisMapConnected{l,1}=zeros(N-(l-1),Nh);
            index=1;

            for i=1:N-(l-1)
                SubsetsConnected{l,1}(i,:)=Subsets{l,1}(index,:);
                BasisMapConnected{l,1}(i,:)=BasisMap{l,1}(index,:);
                index=index+nchoosek(N-i,l-1);
            end


        end
        
        
        BasisMap_ret=BasisMapConnected;
        Subsets_ret=SubsetsConnected;
        
    elseif type ==2
        
        
        SubsetsConnected=cell(N,1);
        BasisMapConnected=cell(N,1);
        
        for l=1:N
            SubsetsConnected{l,1}=zeros(1,l);
            BasisMapConnected{l,1}=zeros(1,Nh);
            
            index=1;
            SubsetsConnected{l,1}(1,:)=Subsets{l,1}(index,:);
            BasisMapConnected{l,1}(1,:)=BasisMap{l,1}(index,:);


        end
        
        BasisMap_ret=BasisMapConnected;
        Subsets_ret=SubsetsConnected;
        
     elseif type==3
        
        Subsets_ret=cell(2,1);
        BasisMap_ret=cell(2,1);
        
        for l=1:1
            Subsets_ret{1,1}=zeros((N-2),l);
            BasisMap_ret{1,1}=zeros((N-2),2^N);
            
            for i=1:N/2
                
                curr_subset=[N/2-(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
           
                Subsets_ret{1,1}(2*i-1,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i-1,:)=BasisMap{l,1}(index,:);
                
                
                curr_subset=[N/2+1+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
                
                
                Subsets_ret{1,1}(2*i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i,:)=BasisMap{l,1}(index,:);
                
            end


        end
        
        for l=2:2
            
            
            Subsets_ret{2,1}=zeros((N-2)/2,l);
            BasisMap_ret{2,1}=zeros((N-2)/2,2^N);
            
            for i=1:N/2
                curr_subset=[N/2-(i-1),N/2+1+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
                
                Subsets_ret{2,1}(i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{2,1}(i,:)=BasisMap{l,1}(index,:);
               
            end
                

        end    
        
        
    elseif type ==4
        
        Subsets_ret=cell(2,1);
        BasisMap_ret=cell(2,1);
        
        
        for l=2:2
            
            Subsets_ret{1,1}=zeros((N-2),l);
            BasisMap_ret{1,1}=zeros((N-2),2^N);
            
            for i=1:(N-2)/2
                
                curr_subset=[N/2-1-(i-1),N/2-(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
  
                Subsets_ret{1,1}(2*i-1,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i-1,:)=BasisMap{l,1}(index,:);
                
                
                curr_subset=[N/2+1+(i-1),N/2+2+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
                
                
                Subsets_ret{1,1}(2*i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i,:)=BasisMap{l,1}(index,:);
                
            end

        end
                
        for l=4:4

            Subsets_ret{2,1}=zeros((N-2)/2,l);
            BasisMap_ret{2,1}=zeros((N-2)/2,2^N);
            
            for i=1:(N-2)/2
                curr_subset=[N/2-1-(i-1),N/2-(i-1),N/2+1+(i-1),N/2+2+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
                
                Subsets_ret{2,1}(i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{2,1}(i,:)=BasisMap{l,1}(index,:);
                
                
            end
                
        end
        
    elseif type==5
        
        
        for l=3:3
            Subsets_ret{1,1}=zeros((N-4),l);
            BasisMap_ret{1,1}=zeros((N-4),2^N);
            
            for i=1:(N-4)/2
                
                curr_subset=[N/2-(i-1)-2,N/2-(i-1)-1,N/2-(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
           
                Subsets_ret{1,1}(2*i-1,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i-1,:)=BasisMap{l,1}(index,:);
                
                curr_subset=[N/2+1+(i-1),N/2+2+(i-1),N/2+3+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
              
                
                Subsets_ret{1,1}(2*i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{1,1}(2*i,:)=BasisMap{l,1}(index,:);
                
            end


        end
        
        for l=6:6

            Subsets_ret{2,1}=zeros((N-4)/2,l);
            BasisMap_ret{2,1}=zeros((N-4)/2,2^N);            
            for i=1:(N-4)/2
                curr_subset=[N/2-(i-1)-2,N/2-(i-1)-1,N/2-(i-1),N/2+1+(i-1),N/2+2+(i-1),N/2+3+(i-1)];
                
                index=ismember(Subsets{l,1}(:,1:l),curr_subset,'rows');
                index=find(index);
                
                Subsets_ret{2,1}(i,:)=Subsets{l,1}(index,:);
                BasisMap_ret{2,1}(i,:)=BasisMap{l,1}(index,:);
                
            end
                

        end
        
        
    else
        BasisMap_ret=BasisMap;
        Subsets_ret=Subsets;
    end
        
end