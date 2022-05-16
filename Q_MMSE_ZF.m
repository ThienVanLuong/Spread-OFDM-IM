function [Q] = Q_MMSE_ZF(N,hh,ZM)      

 % a=avSNR*diag(h(:,jj))*sc;
               Q=zeros(N,1);
                if(ZM==1)
                   % h2=inv(a'*a)*a';
                   Q=hh.^-1;
                else
                    %h2=inv(a'*a+eye(N))*a';
                    for t=1:N
                        Q(t)=conj(hh(t))./(abs(hh(t)).^2+1);
                    end
                end