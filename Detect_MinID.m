function [ID,AcIndex] = Detect_MinID(Y,K,index_allz,p1)

                [~,Bindex] = sort(Y);
                AcIndex = sort(Bindex(1:K),'descend'); 
                AcIndex = AcIndex';
                ID=-1;
                for ii=1:2^p1
                    if(sum(index_allz(ii,:)==AcIndex)==K)
                        ID = ii-1; 
                    end
                end  
                
