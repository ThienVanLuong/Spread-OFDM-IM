
% Thien Van Luong 10/2017 - Queen's University Belfast, UK. Email:
% thienctnb@gmail.com.
% I am now with Phenikaa University, Vietnam. 
% My homepage: www.tvluong.wordpress.com
% The code is based on our my paper:
% [2] T. V. Luong and Y. Ko, "Spread OFDM-IM With Precoding Matrix and Low-Complexity Detection Designs," 
% IEEE Transactions on Vehicular Technology, vol. 67, no. 12, pp. 11619-11626, Dec. 2018.

%% ======================================= Main spread IM program ========================
%
clear
%% Important Paremeters
M=4;
N=4;
K=1;
Ang=M*N;

SScode = 5; % 1= Zadoff, 2=Walsh, 3=DFT, 4=Algebra_OFDM, 5 = rotated WH
Detect_method = 1; %1=ML. 2=MMSE+LLR, 3=New MMSE-LLR, 4=MMSE+GD; 5=LowML ; 6 = MRC

ZM=1; % 1 = ZF, 2 = MMSE: Q-matrix

ch = 2; %1 is LowML/ IP-MMSE, 2 is LowML-new/ EIP-MMSE -> further reduce complexity
roZC = 1; % rotation angle of ZC spreading matrix

%% 
D=1; % New MMSE-LLR
var = 0.05;
mmse = 1;
CSI=1;

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
%% ======================= Spread-IM Parameters ================================
iter = 2;  % Iterations
nSymPerFrame = 1e4; % Number of symbol per frame(1 OFDM symbol)
EbN0dB = 0:5:25;
%EbN0dB = [0 5 10 15 17];
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
%PwrSC=N./(p1+p2);
%% ==================== Loop for SNR =========================
PEP = zeros(1,size(sigma,2)); % index symbol error
OFDM_SER = zeros(1,size(sigma,2)); % ofdm symbol error
Total_SER = zeros(1,size(sigma,2));
BER=zeros(1,size(sigma,2));
BER1=zeros(1,size(sigma,2));
BER2=zeros(1,size(sigma,2));

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
        Dr(k)=exp(1i*2*pi*(k-1)./Ang);
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
        
%% Mary ref_sym
sym_test=zeros(M,1);
        for qq=1:M
            if(Mary==1)
                sym_test(qq)=pskmod(qq-1,M,ro*pi./M,'gray');
            else
                sym_test(qq)=qammod(qq-1,M,0,'gray');
            end
        end  
        
ref_sym = sym_test; 
        if(Mary==1)
            ref_symmm = ref_sym.*(1./abs(ref_sym)); % PSK
        else
            ref_symmm = ref_sym.*(1/sqrt(QAM)); % QAM
        end   
%% ML Mary detection
if(Detect_method==1)
    A=(1:M)';
    for x=1:K-1
        E=cell(1,M);
        for i=1:M
            C=ones(length(A),1)*i;
            E{i}=[C A];
        end
        combs=cell2mat(E');
        A=combs;
    end

ref_sym_ML = sym_test(A);
    %ref_sym = sym_test.';
    if(Mary==1)
        ref_symmm_ML = ref_sym_ML.*(1./abs(ref_sym_ML));
    else
        ref_symmm_ML = ref_sym_ML.*(1/sqrt(QAM));
    end
end

        if(K==2&&N==4)
            index_all = [1 0;2 0;3 1;3 2];
            %index_all = Combin_Md(N,K);
        else
            index_all = Combin_Md(N,K);
        end
index_allz=index_all+1;
%% =============== SNR For Loop ==============================
for s1 = 1:size(sigma,2)    
    fprintf('== EbN0(dB) is %g == \n',EbN0dB(s1))
    %% ==================== Loop for iteration =======================
    symerr_mcik = zeros(1,iter);
    symerr_ofdm = zeros(1,iter);
    symerr_iter= zeros(1,iter);  
    BER_iter= zeros(1,iter); 
    BER_iter_1= zeros(1,iter);
    BER_iter_2= zeros(1,iter);
    for s2 = 1:iter        
        fprintf('== EbN0(dB) is %g and iteration is %g == \n',EbN0dB(s1),s2)        
        %% ===================== Bit generator =========================
        bit = randi([0 1],1,(p1+p2)*nSymPerFrame);
        bit2 = reshape(bit.',p1+p2,nSymPerFrame).';        
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
        
        index_bit = bit2(:,1:p1);
        index_sym = BitoDe(index_bit);
        
        if(Mary==1)
            sym_norm = sym.*(1./abs(sym));
        else
            sym_norm = sym.*(1./sqrt(QAM));
        end
       
        sym_tx = sym_norm.*sqrt(PwrSC);
        
        % transmitted OFDM symbols
        tx_sym1 = zeros(N,nSymPerFrame);
        for kk = 1:nSymPerFrame            
            kk_index = index_sym(kk)+1;            
            indices = index_all(kk_index,:)+1;
            tx_sym1(indices,kk) = sym_tx(kk,:);
        end   
        
        tx_sym=zeros(size(tx_sym1));
        for k=1:nSymPerFrame
            tx_sym(:,k)=sc*tx_sym1(:,k);
        end       
        
 %% Imperfect CSI and Received signals                      
        if(CSI==1)
            eps=0;
        elseif(CSI==2)
            eps=var;
        else
            eps=1./(1+mmse*EsN0(s1));
        end
        
        noise = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)));        
        h = 1/sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym)))*sqrt(1-eps);                     
        e=sqrt(eps)./sqrt(2)*(randn(size(tx_sym))+1i*randn(size(tx_sym))); 
        h1=h+e; 
        y = sqrt(EsN0(s1))*h1.*tx_sym+noise;         
        avSNR=sqrt(EsN0(s1)); 
        
        %% ================== ML / LLR / Greedy detect ====================
        index_sym_de = zeros(1,nSymPerFrame);
        indices_de = zeros(nSymPerFrame,K);
        re_sym = zeros(nSymPerFrame,K);
      
        for jj=1:nSymPerFrame
%% ML detector
            if (Detect_method == 1)   
                [BB,MM] = ML_Detector_Spread_IM(avSNR,M,K,p1,PwrSC,index_all,y,h,N,jj,sc,ref_sym_ML,ref_symmm_ML);
                index_sym_de(jj) = BB-1;
                re_sym(jj,:) = MM;
%% ZF/MMSE Detector + LLR
            elseif (Detect_method == 2)
                hh=avSNR*h(:,jj);
                Q=Q_MMSE_ZF(N,hh,ZM);
                h2=sc'*diag(Q);
                y2=h2*y(:,jj);               
                lamda = zeros(N,1);
                sym_m=zeros(1,N);
                for i = 1:N  
                        A=abs(y2(i)-sqrt(PwrSC).*ref_symmm).^2;
                        [minA,I]=min(A);
                        lamda(i)=abs(y2(i)).^2-minA;
                        sym_m(i)=ref_sym(I);
                end
                [index_sym_de(jj),indices_de(jj,:)] = Detect_MaxID(lamda,N,K,index_allz,p1);          
                re_sym(jj,:) = sym_m(indices_de(jj,:));  
                
%% ======================== ZF/MMSE Detector + LowC ML ---- New MMSE-LLR========== Don't care
            elseif(Detect_method == 3) 
                hh=avSNR*h(:,jj);
                H=diag(h(:,jj))*sc;
                Y=y(:,jj);
                Q=Q_MMSE_ZF(N,hh,ZM);
                
                u=zeros(N,1);
                for i1=1:N
                    hk=(conj(sc(:,i1)).*Q).';
                    z=hk*Y;                 
                    u(i1)=abs(z);
                end                
                wid=zeros(1,2^p1); %% LLR for each index symbol of total 2^p1
                for k=1:2^p1
                    wid(k)=sum(u(index_allz(k,:)));
                end
                [Aw, Bw]=sort(wid);                

                w=zeros(D,1);
                sym_w=zeros(D,K);
                    for i=1:D
                        a=avSNR*H(:,index_allz(Bw(2^p1-i+1),:));
                        if(ch==1)
                            h2=inv(a'*a+eye(K))*a';
                        else
                            h2=sc(:,index_allz(Bw(2^p1-i+1),:))'*diag(Q);
                        end
                        y2=h2*Y;       
                        sym_m=Mary_Decision(y2,K,ref_sym,QAM,PwrSC,M);
                        sym_w(i,:)=sym_m;
                        w(i)=norm(Y-sqrt(PwrSC)*a*sym_m.');
                    end
                    [~,I]=min(w);
                    index_sym_de(jj) = Bw(2^p1-I+1)-1;  %2^p1-I+1
                    re_sym(jj,:) = sym_w(I,:);                    
                    
%%================= ZF/MMSE Detector + GD ==============
            elseif(Detect_method == 4) 
                hh=avSNR*h(:,jj);
                Q=Q_MMSE_ZF(N,hh,ZM);            
                h2=sc'*diag(Q);               
                y2=h2*y(:,jj);
                Y = abs(y2).^2;
                [index_sym_de(jj),indices_de(jj,:)] = Detect_MaxID(Y,N,K,index_allz,p1);
                
                for i1=1:K
                        idx = indices_de(jj,i1);
                        dis = zeros(1,M);
                        for k1=1:M
                            dis(k1)=abs(y2(idx)-sqrt(PwrSC)*ref_symmm(k1)).^2;
                        end
                        [~,I] = min(dis);
                        re_sym(jj,i1) =ref_sym(I);                    
                end   
%% =========================== MMSE - LowML   =========================
            elseif(Detect_method == 5)  
                H=diag(h(:,jj))*sc;
                Y=y(:,jj);
                w=zeros(2^p1,1);                              
                if(K==9)    % 1 or 9          
%% %%%%%%%%%%%%%%%%%%% OS-MMSE for K=1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
                    for i=1:2^p1
                        hk=avSNR*H(:,i);
                        z=hk'*Y;%./(hk'*hk);
                        w(i)=abs(z).^2;
                    end              
                    [~,I]=max(w);
                    index_sym_de(jj) = I-1;                   
                    y2=H(:,I)'*Y; %./(H(:,I)'*H(:,I)+1);
                    re_sym(jj,:)=Mary_Decision(y2,K,ref_sym,QAM,PwrSC,M);
                else
%% %%%%%%%%%%%%%%%%%%% MMSE - LowML for any K %%%%%%%%%%%%%%%%%%%%%%%%%%
                   if(ch==2)
                       hh=avSNR*h(:,jj);
                       Q=Q_MMSE_ZF(N,hh,ZM);
                   end
                    sym_w=zeros(2^p1,K);
                    for i=1:2^p1
                        a=avSNR*H(:,index_allz(i,:));
                        if(ch==1)                           
                            if(K>1)
                                h2=inv(a'*a+eye(K))*a';
                            else
                                h2=a'./(a'*a+1);
                            end
                        else
                            h2=sc(:,index_allz(i,:))'*diag(Q);
                        end
                        y2=h2*Y;   
                        sym_m=Mary_Decision(y2,K,ref_sym,QAM,PwrSC,M);
                        sym_w(i,:)=sym_m;
                        if(Mary==1)
                            w(i)=norm(Y-sqrt(PwrSC)*a*sym_m.');
                        else
                            w(i)=norm(Y-sqrt(PwrSC)*a*sym_m.'./sqrt(QAM));
                        end
                    end
                    [~,I]=min(w);
                    index_sym_de(jj) = I-1;  
                    re_sym(jj,:) = sym_w(I,:);
                end
 %% ===================== % MRC , near ML ===================== Don't care
            else
                H=diag(h(:,jj))*sc;
                Y=y(:,jj);     
                w=zeros(N,1);               
                hh=avSNR*h(:,jj);              
                Q=zeros(N,1);
                for t=1:N
                    Q(t)=conj(hh(t))./(abs(hh(t)).^2+1);
                end
                %% Near ML detection, based on MCIK
                sym_w= zeros(1,N);     
                for i1=1:N
                   % z=zeros(K,1);
                    %hk=H(:,i1);
                    hk=(conj(sc(:,i1)).*Q).';
                    z=hk*Y;%./norm(hk).^2;
                    h3=sqrt(PwrSC)*avSNR*hk*H(:,i1);
                    z=z./h3;
                    v=abs(z-ref_sym).^2; 
                    [~,I]=min(v);
                    sym_w(i1)=ref_sym(I);
                    s=ref_sym(I);                   
                    w(i1)=(abs(s*h3)).^2-2*real(conj(z)*h3*s);
                end
%                
                [index_sym_de(jj),AcIndex] = Detect_MinID(w,K,index_allz,p1);
                re_sym(jj,:)=sym_w(AcIndex);
%                 
            %% OS-MMSE
%                 u=zeros(1,2^p1);
%                 for i1=1:2^p1
%                     id=index_allz(i1,:);
%                     z=zeros(1,K);
%                     for i2=1:K
%                          hk=(conj(sc(:,id(i2))).*Q).';
%                          z(i2)=hk*Y./norm(hk).^2;                        
%                     end
%                     u(i1)=sum(abs(z).^2);
%                 end
%                 [~,I]=max(u);
%                 index_sym_de(jj) = I-1;
%                 AcIndex=index_allz(I,:);

%% MRC
                for i1=1:N
                    hk=(conj(sc(:,i1)).*Q).';
                    z=hk*Y;                 
                    w(i1)=abs(z);
                end
                [index_sym_de(jj),AcIndex] = Detect_MaxID(w,N,K,index_allz,p1);
                                               
%% M-ary detection
                    a=avSNR*H(:,AcIndex);
                    h2=inv(a'*a+eye(K))*a';
                   % h2=sc(:,index_allz(i,:))'*diag(Q);                   
                    y2=h2*Y;
                    re_sym(jj,:)=Mary_Decision(y2,K,ref_sym,QAM,PwrSC,M);
            end
        end
        %% =================error rate computation====================
        % ofdm symbol error
        ofdm_symerr = sum(sum(sym~=re_sym));
        % index symbol error
        ind_symerr = sum(index_sym~=index_sym_de);
        % index symbol to bit, index bit error
        index_bit_de = DetoBit(index_sym_de,p1);
        index_bit_err=sum(sum(index_bit~=index_bit_de));
        
        % QAM symbol to bit
        if(Mary==1)
            info_de_re=pskdemod(re_sym,M,ro*pi./M,'gray');
        else
            info_de_re=qamdemod(re_sym,M,0,'gray');
        end
        info_bit_re= zeros(nSymPerFrame,K*bps);
        for kk=1:K
            info_bit_re(:,(kk-1)*bps+1:kk*bps)=de2bi(info_de_re(:,kk),bps);
        end
        info_bit_err=sum(sum(info_bit~=info_bit_re));
        
        %% ===========symbol & bit error rate  1 iteration==========        
        % MCIK sym error
        symerr_mcik(s2) = ind_symerr/nSymPerFrame;
        % OFDM sym error
        symerr_ofdm(s2) = ofdm_symerr/(K*nSymPerFrame);
        % symbol error rate
        symerr_iter(s2) = (ind_symerr+ofdm_symerr)/(nSymPerFrame+K*nSymPerFrame);
        
        %%% Bit error rate BER
        BER_iter(s2)=(info_bit_err+index_bit_err)./((p1+p2)*nSymPerFrame);
        BER_iter_1(s2) = index_bit_err./p1./nSymPerFrame;
        BER_iter_2(s2) = info_bit_err./p2./nSymPerFrame;
    end    
    
    %% =============average bit error rate================
    PEP(s1) = sum(symerr_mcik)/iter;
    OFDM_SER(s1) = sum(symerr_ofdm)/iter;
    Total_SER(s1) = sum(symerr_iter)/iter; 
    BER(s1)= sum(BER_iter)./iter;
    BER1(s1)= sum(BER_iter_1)./iter;
    BER2(s1)= sum(BER_iter_2)./iter;
end

%% display plot figure
% print choosen Detector and Spreading matrix
if(Detect_method==2)
    Detector = 'MMSE-LLR';
elseif(Detect_method==5)
    if(ch==1)
        Detector = 'LowML';
    else
        Detector = 'LowML-new';
    end
elseif(Detect_method==1)
    Detector='ML';
else
    Detector = 'Unknown';
end
if(SScode==1)
    if(roZC==1)
        Smatrix = 'roZC';
    else
        Smatrix='ZC'
    end
elseif(SScode==2)
    Smatrix = 'WH';
elseif(SScode==3)
    Smatrix='DFT';
elseif(SScode==4)
    Smatrix='AL';
elseif(SScode==5)
    Smatrix='roWH';
else
    Smatrix='Unknown';
end
fprintf('Detector: %s, Smatrix: %s, M = %d, N = %d, K = %d \n',Detector, Smatrix,M,N,K)

SM=OFDM_SER;

figure (3)
if(Detect_method==1)
    semilogy(EbN0dB,BER,'r o-','LineWidth',1.5,'MarkerSize',10)
elseif(Detect_method==2)
    if(ZM==2)
        semilogy(EbN0dB,BER,'b >-','LineWidth',1.5,'MarkerSize',10)
    else
        semilogy(EbN0dB,BER,'b <:','LineWidth',1.5,'MarkerSize',10)
    end
elseif(Detect_method==4)
    if(ZM==2)
        semilogy(EbN0dB,BER,'k *-','LineWidth',1.5,'MarkerSize',10)
    else
        semilogy(EbN0dB,BER,'k +:','LineWidth',1.5,'MarkerSize',10)
    end
else
    if(ZM==2)
        semilogy(EbN0dB,BER,'r *-','LineWidth',1.5,'MarkerSize',10)
    else
        semilogy(EbN0dB,BER,'r o:','LineWidth',1.5,'MarkerSize',10)
    end
end
hold on
BER
axis([0 25 10^-5 10^0])
grid on
hold on
title('')
xlabel('SNR (dB)')
ylabel('BEP')
toc