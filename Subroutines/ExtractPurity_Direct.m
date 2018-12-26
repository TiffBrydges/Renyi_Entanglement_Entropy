function [TrRho2_av,TrRho2_av_std]=ExtractPurity_Direct(qstates,A_Subs,N,CalcErrors)
        
        
        [NU,NM]=size(qstates);
        
        number_of_subsystem_sizes=length(A_Subs);

        TrRho2_av=cell(number_of_subsystem_sizes,1);
        TrRho2_av_std=cell(number_of_subsystem_sizes,1);
        
        for l=1:number_of_subsystem_sizes
            
            sub_size=A_Subs{l,2};
            
            [numberofsubs_l,~,~]=size( A_Subs{l,1});
            TrRho2_av{l}=zeros(numberofsubs_l,1);
            TrRho2_av_std{l}=zeros(numberofsubs_l,1);
            
            for s=1:numberofsubs_l
                
                A=squeeze(A_Subs{l,1}(s,:,:));
                
                purity=zeros(NU,1);
                
                for r=1:NU
                    B=A(qstates(r,:)+1,qstates(r,:)+1);
                    purity(r)=sum(B(:));
                end
                
                if not(CalcErrors)

                    TrRho2_av{l}(s)=mean(purity)/(NM*(NM-1)) - 2^sub_size/(NM-1);
                    
                else
                    % Jackknife estimation of the standard error, note that
                    % since purity estimation is linear this coincides with
                    % the usual standard error of the mean std(purity)/(NM*(NM-1)*sqrt(NU)).
                   
                    jacknifesamples=zeros(NU,1);
                    meanpurity=mean(purity);
   
                    for r=1:NU
                        
                        jacknifesamples(r,1)=(meanpurity - purity(r)/NU)*NU/(NU-1) ;
                    end
 
                    TrRho2_av{l}(s)=mean(jacknifesamples)/(NM*(NM-1)) - 2^sub_size/(NM-1);
                    TrRho2_av_std{l}(s)=std(jacknifesamples,1)/(NM*(NM-1))*sqrt(NU-1);
                end
                
            end
            
        end
end