% ============Maximum Likelihood Detector
function [BB,MM] = ML_Detector_Spread_IM(avSNR,M,K,p1,PwrSC,index_all,y,h,N,jj,sc,ref_sym_ML,ref_symmm_ML)

dis_all = zeros(2^p1,M^K); % all possible realizations, index & symbol
%dis_m = zeros(2^p1,1); %% unused
%Heq=diag(h(:,jj))*sc;
for mm = 1:M^K    
    sym_norm = ref_symmm_ML(mm,:);
    sym_m = sqrt(PwrSC)*sym_norm;
    dis_m = zeros(2^p1,1);
    for bb = 1:2^p1
        pp = index_all(bb,:)+1;
        sym_b = zeros(N,1);
        sym_b(pp) = sym_m;
        sym_b=sc*sym_b;
        tmp = norm(y(:,jj)-avSNR*h(:,jj).*sym_b).^2;
       % tmp = norm(y(:,jj)-avSNR*Heq*sym_b).^2;
        dis_m(bb) = tmp;
    end
    dis_all(:,mm) = dis_m;
end
% find minimun row (symbol)
[tmp0 I] = min(dis_all);
% find minimum column (index)
[tmp1 min_col] = min(tmp0);
% index demodulate
BB = I(min_col);
MM = ref_sym_ML(min_col,:);