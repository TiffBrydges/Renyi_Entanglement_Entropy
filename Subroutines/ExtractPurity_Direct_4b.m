function [TrRho2_av,TrRho2_av_std]=ExtractPurity_Direct_4b(qstates,A_Subs)
        
        
        [NU,NM]=size(qstates);
        
        
        
        number_of_subsystem_sizes=length(A_Subs);

        TrRho2_av=cell(number_of_subsystem_sizes,1);
        TrRho2_av_std=cell(number_of_subsystem_sizes,1);
        
        for l=1:number_of_subsystem_sizes
            
            sub_size=A_Subs{l,2};
            
            [numberofsubs_l,~,~]=size( A_Subs{l,1});
            TrRho2_av{l}=zeros(numberofsubs_l,NU);
            TrRho2_av_std{l}=zeros(numberofsubs_l,NU);
            
            for s=1:numberofsubs_l
                
                A=squeeze(A_Subs{l,1}(s,:,:));
                
                purity=zeros(NU,1);
                
                for r=1:NU
                    B=A(qstates(r,:)+1,qstates(r,:)+1);
                    purity(r)=sum(B(:));
                end
                
                TrRho2_av{l}(s,:)=purity/(NM*(NM-1)) - 2^sub_size/(NM-1);
                    
            end
            
        end
end