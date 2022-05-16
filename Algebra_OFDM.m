function [G] = Algebra_OFDM(N,Tn)

G=zeros(N,N);
for i=1:N
    for k=1:N
        G(i,k)=Tn(i).^(k-1);
    end
end
fa = sqrt(trace(G'*G)./N);
G=G./fa;

