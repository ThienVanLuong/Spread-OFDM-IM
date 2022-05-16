function [sym_m] = Mary_Decision(y2,K,ref_sym,QAM,PwrSC,M)

sym_m=zeros(1,K);
Y2=zeros(M,1);
for k=1:K
                            if(M==4)
                                b=real(y2(k));
                                f=imag(y2(k));
                                if ((b>0)&&(f>0))
                                    sym_m(k)=ref_sym(1);
                                elseif ((b<0)&&(f>0))
                                    sym_m(k)=ref_sym(2);
                                elseif ((b>0)&&(f<0))
                                    sym_m(k)=ref_sym(3);
                                else
                                    sym_m(k)=ref_sym(4);
                                end    
                            elseif(M==2)
                                 b=real(y2(k));
                                 if(b>0)
                                      sym_m(k)=ref_sym(1);
                                 else
                                      sym_m(k)=ref_sym(2);
                                 end
                            elseif(M==8) % 8 PSK
                                pha=angle(y2(k)); 
                                if((pha>-pi/8)&&(pha<pi/8))
                                     sym_m(k)=ref_sym(1);
                                elseif((pha>pi/8)&&(pha<3*pi/8))
                                    sym_m(k)=ref_sym(2);
                                elseif((pha>3*pi/8)&&(pha<5*pi/8))
                                    sym_m(k)=ref_sym(4);
                                elseif((pha>5*pi/8)&&(pha<7*pi/8))
                                    sym_m(k)=ref_sym(3);
                                elseif((pha>-7*pi/8)&&(pha<-5*pi/8))
                                    sym_m(k)=ref_sym(8);
                                elseif((pha>-5*pi/8)&&(pha<-3*pi/8))
                                    sym_m(k)=ref_sym(6);
                                elseif((pha>-3*pi/8)&&(pha<-pi/8))
                                    sym_m(k)=ref_sym(5);
                                else
                                    sym_m(k)=ref_sym(7);
                                end  
                            elseif(M==16) % 16 QAM     
                                y3=y2(k)*sqrt(QAM)./sqrt(PwrSC);
                                b=real(y3);
                                f=imag(y3);
                                if(b>2)
                                    re=3;
                                elseif(b<-2)
                                    re=-3;
                                elseif((b<2)&&(b>0))
                                    re=1;
                                else
                                    re=-1;
                                end
                                if(f>2)
                                    im=3;
                                elseif(f<-2)
                                    im=-3;
                                elseif((f<2)&&(f>0))
                                    im=1;
                                else
                                    im=-1;
                                end
                                sym_m(k)=re+1i*im;    
                            else
                              %  DisY = zeros(M,1);
                                y3=y2(k)*sqrt(QAM)./sqrt(PwrSC);
                                for m=1:M
                                    Y2(m)=y3;
                                    %DisY(m)=abs(y2(k)-ref_sym(m));
                                end
                               % [~,Imin] = min(DisY);
                                [~,Imin] = min(abs(Y2-ref_sym));
                                sym_m(k)=ref_sym(Imin);                                
                            end
end