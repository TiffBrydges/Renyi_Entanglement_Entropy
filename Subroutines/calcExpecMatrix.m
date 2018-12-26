function [A]=calcExpecMatrix(N)

    alpha = sqrt(2-sqrt(3))*(sqrt(3)+2);
    beta = -1./(sqrt(3)+2);

    A=cell(N,1);

    for l=1:N

            A{l} = zeros(2^l,2^l);


            for k=1:2^l
                    for k0=1:2^l


                        dk=pdist2(dec2bin(k-1,l)-'0',dec2bin(k0-1,l)-'0','Minkowski',1);

                        A{l}(k,k0) =A{l}(k,k0)+ alpha^(l)*beta^(dk);

                    end
            end

            A{l}  =  A{l}*transpose(A{l})/2^l ;
    end

end