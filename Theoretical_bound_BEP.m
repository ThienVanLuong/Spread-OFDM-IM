% Thien Van Luong 10/2017 - Queen's University Belfast, UK. Email:
% thienctnb@gmail.com.
% I am now with Phenikaa University, Vietnam. 
% My homepage: www.tvluong.wordpress.com
% The code is based on our my paper:
% [2] T. V. Luong and Y. Ko, "Spread OFDM-IM With Precoding Matrix and Low-Complexity Detection Designs," 
% IEEE Transactions on Vehicular Technology, vol. 67, no. 12, pp. 11619-11626, Dec. 2018.

%% ============================ Theoretical bound on BER of spread OFDM-IM ================
%
clear

M=4;
N=4;
K=1;

SScode = 5; % 1= Zadoff, 2=Walsh, 3=DFT, 4=Algebra_OFDM, 5 = rotated WH
roZC = 0; % rotation angle of ZC spreading matrix

if(M==4)
    ro=1;
else
    ro=0;
end

%Mary=2; % 1 PSK, 2 QAM
if(M>8)
    Mary=2;
else
    Mary=1;
end
if(M==8)
    QAM = (5*M-4)./6;
else
    QAM = (2/3)*(M-1);
end

tic
%% ======================= Misc Parameters ================================
EbN0dB = 0:5:25;
EbN0 = 10.^(EbN0dB/10);

% Es/N0 parameter
PwrSC = N/K; % Average Tx power per active sub-carrier
bps = log2(M); % bits per symbol 
EsN0dB = EbN0dB; % + 10*log10(bps);%+10*log10(1/PwrSC);
EsN0 = 10.^(EsN0dB/10);
c = 2^floor(log2(nchoosek(N,K))); % Effective Carrier Combinations
p1 = floor(log2(nchoosek(N,K)));  % index bit length per cluster
p2 = K*bps; % information bit length per cluster
p2_re=bps;
p=p1+p2;
SE=p/N;
sigma = sqrt(1./EsN0);
T=c*M.^K;

%% Spreading matrix generator
% Algebraic Constructions
if(N==4) % for N = 3,4,5,6,8
    Tn=[exp(-1i*pi./8) exp(-1i*5*pi./8) exp(-1i*9*pi./8) exp(-1i*13*pi./8)];
elseif(N==5)
    Tn=[2.^(1/5)*exp(-1i*pi./20) 2.^(1/5)*exp(-1i*9*pi./20) 2.^(1/5)*exp(-1i*17*pi./20) 2.^(1/5)*exp(-1i*25*pi./20) 2.^(1/5)*exp(-1i*33*pi./20)];
elseif(N==6)
    Tn=[exp(-1i*2*pi./7) exp(-1i*4*pi./7) exp(-1i*6*pi./7) exp(-1i*8*pi./7) exp(-1i*10*pi./7) exp(-1i*12*pi./7)];
elseif(N==3)
    Tn=[2.^(1/3)*exp(-1i*pi./12) 2.^(1/3)*exp(-1i*9*pi./12) 2.^(1/5)*exp(-1i*17*pi./12)];
elseif(N==2)
    Tn=[exp(-1i*pi./4) exp(-1i*5*pi./4)];
else
    Tn=[exp(-1i*pi./16) exp(-1i*5*pi./16) exp(-1i*9*pi./16) exp(-1i*13*pi./16) exp(-1i*17*pi./16) exp(-1i*21*pi./16) exp(-1i*25*pi./16) exp(-1i*29*pi./16)];
end

Dr=zeros(1,N);
if(1==1)
    for k=1:N
        Dr(k)=exp(1i*2*pi*(k-1)./N./M);
    end
end
sc = zeros(N,N); 
        if(SScode==1)
            u_Zad = 1;
            Zad = zeros(N,1);
            for kk=1:N
                if(mod(N,2)==0)
                    Zad(kk) = exp(-1j*pi*u_Zad*kk.^2./N);
                else
                    Zad(kk) = exp(-1j*pi*u_Zad*kk*(kk+1)./N);
                end
            end          
            for kk=1:N
                if(roZC==0)
                    sc(:,kk)=circshift(Zad,kk);
                else
                    sc(:,kk)=Dr(kk)*circshift(Zad,kk);
                end
            end
        elseif(SScode==2)
            sc = hadamard(N);
        elseif(SScode==3)
            sc = dftmtx(N);
        elseif(SScode==4)
            sc = sqrt(N)*Algebra_OFDM(N,Tn);
        else
            sc=hadamard(N)*diag(Dr);
        end
        
sc = sc./sqrt(N); 

        if(K==2&&N==4)
            index_all = [1 0;2 0;3 1;3 2];
            %index_all = Combin_Md(N,K);
        else
            index_all = Combin_Md(N,K);
        end

%% ==================== Loop for SNR =========================
BER_theo = zeros(1,size(sigma,2)); % index symbol error

dd=[0:T-1]';
bit2=de2bi(dd,p);
X=zeros(N,T);
Z=zeros(N,T);
tic
for s1 = 1:size(sigma,2) 
    fprintf('== EbN0(dB) is %g == \n',EbN0dB(s1))
    avSNR=sqrt(EsN0(s1)); 
        info_bit = bit2(:,p1+1:end);           
        sym=[];
        x=1;
        for i=1:K
            y=bps*i;
            info_bit_i= info_bit(:,x:y);
            x=y+1;
            info_dec_i = bi2de(info_bit_i);
           % sym_i = sym_test(info_dec_i+1);
           if(Mary==1)
                sym_i = pskmod(info_dec_i,M,ro*pi./M,'gray');
           else
                sym_i = qammod(info_dec_i,M,0,'gray');
           end
            sym(:,i)=sym_i;
        end        
        % index bits (p1)
        index_bit = bit2(:,1:p1);
        % index symbol ( bit to decimal ), select indices from combinatorial method
        index_sym = BitoDe(index_bit);
        
        % Set average symbol power to 1
        if(Mary==1)
            sym_norm = sym.*(1./abs(sym));
        else
            sym_norm = sym.*(1./abs(QAM));
        end
       
        % Power reallocation
        sym_tx = sym_norm.*sqrt(PwrSC);
        % transmitted OFDM symbols
        tx_sym1 = zeros(N,T);
        for kk = 1:T            
            kk_index = index_sym(kk)+1;            
            indices = index_all(kk_index,:)+1;
            tx_sym1(indices,kk) = sym_tx(kk,:);
        end   
       % X=avSNR*tx_sym1;
        tx_sym=zeros(size(tx_sym1));
        for k=1:T
            tx_sym(:,k)=sc*tx_sym1(:,k);
        end
        Z=avSNR*tx_sym;
        
        sum=0;
         for i=1:T
             for k=1:T
                 Dis=abs(Z(:,i)-Z(:,k)).^2;
                 w=length(find(bit2(i,:)~=bit2(k,:)));
                 F=@(x)(sin(x).^2./(sin(x).^2+Dis(1)./4))*(sin(x).^2./(sin(x).^2+Dis(2)./4))*(sin(x).^2./(sin(x).^2+Dis(3)./4))*(sin(x).^2./(sin(x).^2+Dis(4)./4));
                 ite=quadv(F,1.0e-12,pi/2,1.0e-12)./pi;
                 sum=sum+ite*w;
             end
         end
         BER_theo(s1)=sum./T./p;              
end  
figure(2)
semilogy(EbN0dB,BER_theo,'k >:','LineWidth',1.5,'MarkerSize',10)
hold  on
toc