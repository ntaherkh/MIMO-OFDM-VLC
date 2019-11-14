clc
clear
close all

SNR= [0 10 15 20 25 30 35 40];
loop=10;
N=128;
N_t=130;
sishma_orginal=0.276;
M_qam=16;
M_PPM=4;
N_sub=N/4;
I_H=1;
% lambda_h=2;
% lambda_l=2;

% I_OFF_org=(lambda_l*sishma_tuned_org);
% I_ON_org=(I_H*0.5)+(lambda_l*sishma_tuned_org);
% I_ON_org=(I_H)-(lambda_h*sishma_tuned_org);
% delta_last=0.15*floor((min((1/((3*lambda_l/delta_org)+(lambda_h/delta_org))),(1/((3*lambda_h/delta_org)+(lambda_l/delta_org)))))/0.15);
tau=M_PPM/N_sub;
Count_location=nchoosek([1:N_sub],M_PPM); %  can we optimize N_sub and M_PPM ??
b_MPPM=floor(log2(length(Count_location))); %% bit in MPPM
M_location=Count_location(1:2^b_MPPM,:);
b_dco=(N/2)*log2(M_qam); %%  bit in dco

phi_plot=[1.8:0.4:3.4];
etta_plot=[0.8:0.2:1.2];
delta_plot=[ 10 12 14  16 18] ;

for delta_inx=1:length(delta_plot)
    delta_org=delta_plot(delta_inx);
    sishma_tuned_org=(1/delta_org)*I_H;
    phi_index=1;
    for phi=1.8:0.4:3.4
        phi
        sishma_tuned=sishma_tuned_org;
        normalize_coff=sishma_tuned/sishma_orginal;
        sishma=sishma_tuned;
        sishma_tuned_vect(phi_index)=sishma_tuned;
        I_OFF=(phi*sishma_tuned);
        I_ON=(I_H)-(phi*sishma_tuned);
        interval=(I_ON-I_OFF)/2;
        interval_vect(phi_index)=interval;
        etta_index=1;
        lambda_l=phi;
        lambda_h=phi;
        for etta=1:1:1
            %         etta
            I_m=I_OFF+(interval*etta);
            %     I_L=(1-delta)*I_H;
            p_1=(sishma^2)*(((lambda_l^2)*qfunc(lambda_l))+((lambda_h^2)*qfunc(lambda_h)));
            p_2=(sishma^2)*(1/sqrt(2*pi))*(-((exp(-(lambda_l^2)/2))*lambda_l)-((exp(-(lambda_h^2)/2))*lambda_h));
            p_3=(sishma^2)*(1-qfunc(lambda_l)-qfunc(lambda_h));
            mean_clipped =(sishma)*(((1/sqrt(2*pi))*((exp(-(lambda_l^2)/2))-(exp(-(lambda_h^2)/2))))-lambda_l*(qfunc(lambda_l))+lambda_h*qfunc(lambda_h));
            sq_clipped=p_1+p_2+p_3;
            sishma_clipped=sqrt(sq_clipped);
            sishma_clipped_vect(phi_index)=sishma_clipped;
            mean_clipped_vect(phi_index)=mean_clipped ;
            PP_index=0;
            PP_vect=0;
            for a=1:size(SNR,2)
                sumerr=0;
                bc_1=b_MPPM*(N/N_sub);
                bc_2=b_dco;
                snr_p=SNR(a);
                snr_li=sqrt(10^(snr_p/10));
                E_n=SNR(a)-0;
                E_n_lin=10^(E_n/10);
                E_MPPM_DCO=sishma_clipped^2+ ((tau*(I_ON^2))+((1-tau)*(I_OFF^2)));
                Power_VECTOR(phi_index)=E_MPPM_DCO;
                %             stadeviation_n_1=sqrt((1/(log2(M_qam)))*(E_MPPM_DCO/E_n_lin));
                stadeviation_n_1=sqrt(E_MPPM_DCO/E_n_lin);
                I_mppm=I_OFF*ones(1,N_t);
                %      I_mppm(1)=0;
                %      I_mppm(N_t/2)=0;
                out_en=zeros(1,N_t);
                xe_dco=zeros(1,N_t);
                xe_mppm_dis=zeros(1,N_t);
                xe_mppm=zeros(1,N_t);
                thr_l=sishma*lambda_l;
                thr_h=sishma*lambda_h;
                cnt=0;
                for jj=1:loop
                    I_mppm=I_OFF*ones(1,N_t);
                    %      I_mppm(1)=0;
                    %      I_mppm(N_t/2)=0;
                    out_en=zeros(1,N_t);
                    xe_dco=zeros(1,N_t);
                    xe_mppm_dis=zeros(1,N_t);
                    xe_mppm=zeros(1,N_t);
                    b_1=randi([0 1],bc_1,1);
                    b_2=randi([0 1],bc_2,1);
                    %   MPPM modulation
                    smppm_1=bi2de(b_1(1:b_MPPM)')+1;
                    loc_ON_1=M_location(smppm_1,:)+1;
                    %          1
                    smppm_2=bi2de(b_1(b_MPPM+1:2*b_MPPM)')+1;
                    loc_ON_2=M_location(smppm_2,:)+33;
                    %          2
                    smppm_3=bi2de(b_1(2*b_MPPM+1:3*b_MPPM)')+1;
                    loc_ON_3=M_location(smppm_3,:)+66;
                    %          3
                    smppm_4=bi2de(b_1(3*b_MPPM+1:4*b_MPPM)')+1;
                    loc_ON_4=M_location(smppm_4,:)+98;
                    %          4
                    ON_vect_Mppm=[loc_ON_1 loc_ON_2 loc_ON_3 loc_ON_4];
                    %  DCO modulation
                    dco_data=bi2de(reshape(b_2',log2(M_qam),N/2)');
                    in=qammod(dco_data,M_qam);
                    out_en(2:N_t/2) =in;
                    out_en(1)=0;
                    out_en(N_t/2+1)=0;
                    I_mppm(ON_vect_Mppm)=I_ON;
                    for i=N_t/2+1:N_t
                        out_en(i)=conj(out_en(N_t+2-i));
                    end
                    g=normalize_coff*(ifft(out_en));
                    g_c=zeros(1,size(g,2));
                    g_c=g;
                    for iw=1:size(g,2)
                        if(g(iw)<=-thr_l)
                            g_c(iw)=-thr_l;
                            cnt=cnt+1;
                        end
                    end
                    for iw=1:size(g,2)
                        if(g(iw)>=thr_h)
                            g_c(iw)=thr_h;
                            cnt=cnt+1;
                        end
                    end
                    I_dco_t=g_c;
                    I_t=I_dco_t+I_mppm;
                    PP_index=PP_index+1;
                    PP_vect(PP_index)=max(I_t);
                    noise_code=normrnd(0,stadeviation_n_1,[1 N_t]);
                    %          noise_code=0;
                    re_out=I_t+noise_code;
                    %          re_out(1)=0;
                    %          re_out(N_t)=0;
                    xe_dco(1)=re_out(1)-I_OFF;
                    for j=2:N_t
                        if(re_out(j)>= I_m)
                            xe_dco(j)=re_out(j)-I_ON;
                            xe_mppm(j)=I_ON;
                            xe_mppm_dis(j)=abs(I_ON-re_out(j));
                        else
                            xe_dco(j)=re_out(j)-I_OFF;
                            xe_mppm(j)=I_OFF;
                            xe_mppm_dis(j)=abs(I_ON-re_out(j));
                            %                   can we detect in a better way? like by comparing to I_OFF
                        end
                    end
                    Ie_t=xe_dco+xe_mppm;
                    %          MPPM Estimation
                    P_1=xe_mppm_dis(2:N_sub+1);
                    P_2=xe_mppm_dis(N_sub+2:N_sub*2+1);
                    P_3=xe_mppm_dis(2*N_sub+3:3*N_sub+2);
                    P_4=xe_mppm_dis(3*N_sub+3:4*N_sub+2);
                    L1=sort(P_1);
                    L2=sort(P_2);
                    L3=sort(P_3);
                    L4=sort(P_4);
                    %                     [B_1,I_1]=mink(P_1,M_PPM);
                    %                     [B_2,I_2]=mink(P_2,M_PPM);
                    %                     [B_3,I_3]=mink(P_3,M_PPM);
                    %                     [B_4,I_4]=mink(P_4,M_PPM);
                    %
                    vec1=find(P_1==L1(1) | P_1==L1(2) | P_1==L1(3)| P_1==L1(4));
                    %                     %          if( (vec1(1)>= 15) | ((vec1(1)== 14)&(vec1(2)> 21)))
                    %                     %              vec1(1)=13;
                    %                     %          end
                    vec2=find(P_2==L2(1) | P_2==L2(2) | P_2==L2(3)| P_2==L2(4));
                    %                     %          if( (vec2(1)>= 15) | ((vec2(1)== 14)&(vec2(2)> 21)))
                    %                     %              vec2(1)=13;
                    % %                     %          end
                    vec3=find(P_3==L3(1) | P_3==L3(2) | P_3==L3(3)| P_3==L3(4));
                    %                     %          if( (vec3(1)>= 15) | ((vec3(1)== 14)&(vec3(2)> 21)))
                    % %                     %              vec3(1)=13;
                    %                     %          end
                    vec4=find(P_4==L4(1) | P_4==L4(2) | P_4==L4(3)| P_4==L4(4));
                    %                     %          if( (vec4(1)>= 15) | ((vec4(1)== 14)&(vec4(2)> 21)))
                    %                     %              vec4(1)=13;
                    %          end
                    % %                     ext_bit_1=de2bi((find(Count_location(:,1)==I_1(1) & Count_location(:,2)==I_1(2)& Count_location(:,3)==vec1(3)& Count_location(:,4)==vec1(4)& Count_location(:,5)==vec1(5)& Count_location(:,6)==vec1(6)& Count_location(:,7)==vec1(7)& Count_location(:,8)==vec1(8)))-1,b_MPPM+1);
                    %                     ext_bit_1=find(I_1==Count_location);
                    %                     ext_bit_2=find(I_2==Count_location);
                    %                     ext_bit_3=find(I_3==Count_location);
                    %                     ext_bit_4=find(I_4==Count_location);
                    ext_bit_1=de2bi((find(Count_location(:,1)==vec1(1) & Count_location(:,2)==vec1(2)& Count_location(:,3)==vec1(3)& Count_location(:,4)==vec1(4)))-1,b_MPPM+1);
                    ext_bit_2=de2bi((find(Count_location(:,1)==vec2(1) & Count_location(:,2)==vec2(2)& Count_location(:,3)==vec2(3)& Count_location(:,4)==vec2(4)))-1,b_MPPM+1);
                    ext_bit_3=de2bi((find(Count_location(:,1)==vec3(1) & Count_location(:,2)==vec3(2)& Count_location(:,3)==vec3(3)& Count_location(:,4)==vec3(4)))-1,b_MPPM+1);
                    ext_bit_4=de2bi((find(Count_location(:,1)==vec4(1) & Count_location(:,2)==vec4(2)& Count_location(:,3)==vec4(3)& Count_location(:,4)==vec4(4)))-1,b_MPPM+1);
                    ext_bit_total=[ext_bit_1(1:b_MPPM)  ext_bit_2(1:b_MPPM)  ext_bit_3(1:b_MPPM)  ext_bit_4(1:b_MPPM) ];
                    %    DCO Estimation
                    out_fft_dco=fft(xe_dco/normalize_coff);
                    out_dco=out_fft_dco(2:N_t/2);
                    f_out_dco=qamdemod(out_dco,M_qam);
                    bit_out_dco=de2bi(f_out_dco);
                    reshaped_bit_dco=reshape(bit_out_dco',N/2*log2(M_qam),1);
                    [ber_2]=biterr(b_2,reshaped_bit_dco);
                    [ber_1]=biterr(b_1,ext_bit_total');
                    sumerr=sumerr + ber_1 +ber_2;
                end
                
                cnt_vect(delta_inx,phi_index,a)=cnt;
                err_rate(delta_inx,phi_index,a)=sumerr/(loop*(bc_1+bc_2));
                err_mppm(delta_inx,phi_index,a)=ber_1;
                err_dco(delta_inx,phi_index,a)=ber_2;
            end
%             figure
%             mean_max_vect(phi_index)=mean(PP_vect);
%             magnit_max_vect(phi_index)=max(PP_vect);
%             plot(Ie_t,'r');hold on;plot(I_t,'b');hold on;plot(I_OFF*ones(1,N_t),'g');hold on; plot(I_ON*ones(1,N_t),'g');plot(I_m*ones(1,N_t));
            PAPR_vect=(PP_vect.^2)/E_MPPM_DCO;
            PAPR_vect_dB=10*log10(PAPR_vect);
            length_PAPR_vect=length(PAPR_vect_dB);
            in_dB=0;
            for dB=0:0.5:15
                in_dB=in_dB+1;
                PAPR_CCDF(in_dB,phi_index,delta_inx)=(length(find(PAPR_vect_dB > dB)))/length_PAPR_vect;
            end
        end
        phi_index=phi_index+1;
    end
end

xtime_c=SNR;

for i=1:length(phi_plot)
figure
semilogy(SNR,squeeze(err_rate(1,i,:)),'b--');hold on;semilogy(SNR,squeeze(err_rate(2,i,:)),'r--');hold on;semilogy(SNR,squeeze(err_rate(3,i,:)),'k--')
hold on;semilogy(SNR,squeeze(err_rate(4,i,:)),'g--');hold on;semilogy(SNR,squeeze(err_rate(5,i,:)),'c--')
title([ 'phi' '=' num2str(phi_plot(i)) ])
end

for j=1:length(delta_plot)
figure
semilogy(SNR,squeeze(err_rate(j,1,:)),'b--');hold on;semilogy(SNR,squeeze(err_rate(j,2,:)),'r--');hold on;semilogy(SNR,squeeze(err_rate(j,3,:)),'k--')
hold on;semilogy(SNR,squeeze(err_rate(j,4,:)),'g--');hold on;semilogy(SNR,squeeze(err_rate(j,5,:)),'c--')
title([ 'delta' '=' num2str(delta_plot(j))])
end
dB=[0:0.5:15];
for j=1:length(phi_plot)
figure
semilogy(dB,squeeze(PAPR_CCDF(:,j,1)),'b--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,2)),'r--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,3)),'k--')
hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,4)),'g--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,5)),'c--')
title([ 'phi' '=' num2str(phi_plot(j))])
end

for j=1:length(delta_plot)
figure
semilogy(dB,squeeze(PAPR_CCDF(:,1,j)),'b--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,2,j)),'r--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,3,j)),'k--')
hold on;semilogy(dB,squeeze(PAPR_CCDF(:,4,j)),'g--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,5,j)),'c--');hold on;
title([ 'delta' '=' num2str(delta_plot(j))  ])
end

