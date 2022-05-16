function [ID,AcIndex] = Detect_MaxID(Y,N,K,index_allz,p1)

                [~,Bindex] = sort(Y);
                AcIndex = sort(Bindex(N-K+1:N),'descend'); 
                AcIndex = AcIndex';
                ID = -1;
                for ii=1:2.^p1
                    if(sum(index_allz(ii,:)==AcIndex)==K)
                        ID = ii-1;   
                    end
                end