clc
clear
close all

SNR= [0 10  20  30  40 45];
loop=100
N=128;
N_t=130;

sishma_original_2=0.2412; % 170
% 162 0.5062  64 qam
% 0.2486; % N=178 16QAM
sishma_original_1=  0.2104  % 130

% 130 0.4528 64 qam
% 0.1956; % N=130  16QAM
M_qam=16;
M_PPM=12;
N_sub=N/8;
I_H=1;
tau=M_PPM/N_sub;
Count_location=nchoosek([1:N_sub],M_PPM); %  can we optimize N_sub and M_PPM ??
b_MPPM=floor(log2(length(Count_location))); %% bit in MPPM
M_location=Count_location(1:2^b_MPPM,:);
DCO_bit_extra=(N/N_sub)*b_MPPM;
b_dco=(N/2)*log2(M_qam); %%  bit in dco
b_dco_2=(b_dco+DCO_bit_extra);
N_2=floor(b_dco_2/(log2(M_qam)));
N_half_2=N_2+1;
phi_plot=[2:0.4:3.2];
etta_plot=[0.8:0.2:1.2];
delta_plot=[ 8 12 14  16 22] ;

for delta_inx=1:length(delta_plot)
    delta_org=delta_plot(delta_inx);
    sishma_tuned_2=(1/delta_org)*I_H;
    sishma_tuned_1=sishma_tuned_2*(sishma_original_1/sishma_original_2);
    phi_index=1;
    for phi=2:0.4:3.2
        phi
        normalize_coff_1=sishma_tuned_1/sishma_original_1;
        normalize_coff_2=sishma_tuned_2/sishma_original_2;
        sishma_1=sishma_tuned_1;
        sishma_2=sishma_tuned_2;
        sishma_tuned_vect(phi_index)=sishma_tuned_1;
        I_OFF=(phi*sishma_tuned_1);
        I_ON=(I_H)-(phi*sishma_tuned_1);
        interval=(I_ON-I_OFF)/2;
        interval_vect(phi_index)=interval;
        etta_index=1;
        lambda_l=phi;
        lambda_h=phi;
        I_BIAS=(I_H*0.5);
        alpha=(0.5*I_H)/(lambda_h*sishma_2);
        for etta=1:1:1
            %         etta
            I_m=I_OFF+(interval*etta);
            %     I_L=(1-delta)*I_H;
            p_1=(sishma_1^2)*(((lambda_l^2)*qfunc(lambda_l))+((lambda_h^2)*qfunc(lambda_h)));
            p_2=(sishma_1^2)*(1/sqrt(2*pi))*(-((exp(-(lambda_l^2)/2))*lambda_l)-((exp(-(lambda_h^2)/2))*lambda_h));
            p_3=(sishma_1^2)*(1-qfunc(lambda_l)-qfunc(lambda_h));
            mean_clipped =(sishma_1)*(((1/sqrt(2*pi))*((exp(-(lambda_l^2)/2))-(exp(-(lambda_h^2)/2))))-lambda_l*(qfunc(lambda_l))+lambda_h*qfunc(lambda_h));
            sq_clipped=p_1+p_2+p_3;
            sishma_clipped=sqrt(sq_clipped);
            sishma_clipped_vect(phi_index)=sishma_clipped;
            % %
            p_1=(sishma_2^2)*(((lambda_l^2)*qfunc(lambda_l))+((lambda_h^2)*qfunc(lambda_h)));
            p_2=(sishma_2^2)*(1/sqrt(2*pi))*(-((exp(-(lambda_l^2)/2))*lambda_l)-((exp(-(lambda_h^2)/2))*lambda_h));
            p_3=(sishma_2^2)*(1-qfunc(lambda_l)-qfunc(lambda_h));
            mean_clipped =(sishma_2)*(((1/sqrt(2*pi))*((exp(-(lambda_l^2)/2))-(exp(-(lambda_h^2)/2))))-lambda_l*(qfunc(lambda_l))+lambda_h*qfunc(lambda_h));
            sq_clipped=p_1+p_2+p_3;
            sishma_clipped_2=sqrt(sq_clipped);
            sishma_clipped_vect_2(phi_index)=sishma_clipped_2;
            mean_clipped_vect_2(phi_index)=mean_clipped ;
            % %
            PP_index=0;
            PP_vect=0;
            for a=1:size(SNR,2)
                sumerr=0;
                sumerr_2=0;
                bc_1=b_MPPM*(N/N_sub);
                bc_2=b_dco;
                bc_3=N_2*log2(M_qam);
                snr_p=SNR(a);
                snr_li=sqrt(10^(snr_p/10));
                E_n=SNR(a)-0;
                E_n_lin=10^(E_n/10);
                E_MPPM_DCO=sishma_clipped^2+ ((tau*(I_ON^2))+((1-tau)*(I_OFF^2)));
                E_DCO_2=(I_BIAS^2)+(alpha^2)*(sishma_clipped_2^2);
                Power_VECTOR(phi_index)=E_MPPM_DCO;
                %             stadeviation_n_1=sqrt((1/(log2(M_qam)))*(E_MPPM_DCO/E_n_lin));
                stadeviation_n_1=sqrt(E_MPPM_DCO/E_n_lin);
                stadeviation_n_2=sqrt(E_DCO_2/E_n_lin);
                I_mppm=I_OFF*ones(1,N_t);
                %      I_mppm(1)=0;
                %      I_mppm(N_t/2)=0;
                out_en=zeros(1,N_t);
                xe_dco=zeros(1,N_t);
                out_en_2=zeros(1,2*N_half_2);
                xe_dco_2=zeros(1,2*N_half_2);
                xe_mppm_dis=zeros(1,N_t);
                xe_mppm=zeros(1,N_t);
                thr_l_1=sishma_1*lambda_l;
                thr_h_1=sishma_1*lambda_h;
                thr_l_2=sishma_2*lambda_l;
                thr_h_2=sishma_2*lambda_h;
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
                    b_3=randi([0 1],bc_3,1);
                    %   MPPM modulation
                    smppm_1=bi2de(b_1(1:b_MPPM)')+1;
                    loc_ON_1=M_location(smppm_1,:)+1;
                    %          1
                    smppm_2=bi2de(b_1(b_MPPM+1:2*b_MPPM)')+1;
                    loc_ON_2=M_location(smppm_2,:)+17;
                    %          2
                    smppm_3=bi2de(b_1(2*b_MPPM+1:3*b_MPPM)')+1;
                    loc_ON_3=M_location(smppm_3,:)+33;
                    %          3
                    smppm_4=bi2de(b_1(3*b_MPPM+1:4*b_MPPM)')+1;
                    loc_ON_4=M_location(smppm_4,:)+49;
                    
                    smppm_5=bi2de(b_1(4*b_MPPM+1:5*b_MPPM)')+1;
                    loc_ON_5=M_location(smppm_5,:)+66;
                    
                    smppm_6=bi2de(b_1(5*b_MPPM+1:6*b_MPPM)')+1;
                    loc_ON_6=M_location(smppm_6,:)+82;
                    
                    smppm_7=bi2de(b_1(6*b_MPPM+1:7*b_MPPM)')+1;
                    loc_ON_7=M_location(smppm_7,:)+98;
                    
                    smppm_8=bi2de(b_1(7*b_MPPM+1:8*b_MPPM)')+1;
                    loc_ON_8=M_location(smppm_8,:)+114;
                    %          4
                    ON_vect_Mppm=[loc_ON_1 loc_ON_2 loc_ON_3 loc_ON_4 loc_ON_5 loc_ON_6 loc_ON_7 loc_ON_8];
                    %  DCO modulation
                    dco_data=bi2de(reshape(b_2',log2(M_qam),N/2)');
                    dco_data_2=bi2de(reshape(b_3',log2(M_qam),N_2)');
                    in=qammod(dco_data,M_qam);
                    in_2=qammod(dco_data_2,M_qam);
                    out_en(2:N_t/2) =in;
                    out_en_2(2:N_half_2)=in_2;
                    out_en(1)=0;
                    out_en_2(1)=0;
                    out_en(N_t/2+1)=0;
                    out_en_2(N_half_2)=0;
                    I_mppm(ON_vect_Mppm)=I_ON;
                    for i=N_t/2+1:N_t
                        out_en(i)=conj(out_en(N_t+2-i));
                    end
                    for i=N_half_2:2*N_half_2
                        out_en_2(i)=conj(out_en_2(2*N_half_2+2-i));
                    end
                    g=normalize_coff_1*(ifft(out_en));
                    g_2=normalize_coff_2*(ifft(out_en_2));
                    g_c=zeros(1,size(g,2));
                    g_c=g;
                    g_c_2=zeros(1,size(g_2,2));
                    g_c_2=g_2;
                    for iw=1:size(g,2)
                        if(g(iw)<=-thr_l_1)
                            g_c(iw)=-thr_l_1;
                            cnt=cnt+1;
                        end
                        if(g(iw)>=thr_h_1)
                            g_c(iw)=thr_h_1;
                            cnt=cnt+1;
                        end
                    end
                    for iw=1:size(g_2,2)
                        if(g_2(iw)<=-thr_l_2)
                            g_c_2(iw)=-thr_l_2;
                            cnt=cnt+1;
                        end
                        if(g_2(iw)>=thr_h_2)
                            g_c_2(iw)=thr_h_2;
                            cnt=cnt+1;
                        end
                    end
                    I_dco_t=g_c;
                    I_dco_t_2=g_c_2;
                    I_t=I_dco_t+I_mppm;
                    PP_index=PP_index+1;
                    PP_vect(PP_index)=max(I_t);
                    PP_vect_2(PP_index)=max(alpha*I_dco_t_2+I_BIAS);
                    noise_code=normrnd(0,stadeviation_n_1,[1 N_t]);
                    noise_code_2=normrnd(0,stadeviation_n_2,[1 2*N_half_2]);
                    %          noise_code=0;
                    re_out=I_t+noise_code;
                    re_out_2=I_dco_t_2+(noise_code_2/alpha);
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
                        P_3=xe_mppm_dis(N_sub*2+2:3*N_sub+1);
                        P_4=xe_mppm_dis(3*N_sub+2:4*N_sub+1);
                        P_5=xe_mppm_dis(4*N_sub+3:5*N_sub+2);
                        P_6=xe_mppm_dis(5*N_sub+3:6*N_sub+2);
                        P_7=xe_mppm_dis(6*N_sub+3:7*N_sub+2);
                        P_8=xe_mppm_dis(7*N_sub+3:8*N_sub+2);
                        [V1,I1]=mink(P_1,M_PPM);
                        [V2,I2]=mink(P_2,M_PPM);
                        [V3,I3]=mink(P_3,M_PPM);
                        [V4,I4]=mink(P_4,M_PPM);
                        [V5,I5]=mink(P_5,M_PPM);
                        [V6,I6]=mink(P_6,M_PPM);
                        [V7,I7]=mink(P_7,M_PPM);
                        [V8,I8]=mink(P_8,M_PPM);
                        vec1=sort(I1);
                        vec2=sort(I2);
                        vec3=sort(I3);
                        vec4=sort(I4);
                        vec5=sort(I5);
                        vec6=sort(I6);
                        vec7=sort(I7);
                        vec8=sort(I8);
                        d1=ismember(Count_location,vec1,'row');
                        ext_bit_1=de2bi((find(d1~=0))-1,b_MPPM+1);
                        d2=ismember(Count_location,vec2,'row');
                        ext_bit_2=de2bi((find(d2~=0))-1,b_MPPM+1);
                        d3=ismember(Count_location,vec3,'row');
                        ext_bit_3=de2bi((find(d3~=0))-1,b_MPPM+1);
                        d4=ismember(Count_location,vec4,'row');
                        ext_bit_4=de2bi((find(d4~=0))-1,b_MPPM+1);
                        d5=ismember(Count_location,vec5,'row');
                        ext_bit_5=de2bi((find(d5~=0))-1,b_MPPM+1);
                        d6=ismember(Count_location,vec6,'row');
                        ext_bit_6=de2bi((find(d6~=0))-1,b_MPPM+1);
                        d7=ismember(Count_location,vec7,'row');
                        ext_bit_7=de2bi((find(d7~=0))-1,b_MPPM+1);
                        d8=ismember(Count_location,vec8,'row');
                        ext_bit_8=de2bi((find(d8~=0))-1,b_MPPM+1);
                        ext_bit_total=[ext_bit_1(1:b_MPPM)  ext_bit_2(1:b_MPPM)  ext_bit_3(1:b_MPPM)  ext_bit_4(1:b_MPPM)  ext_bit_5(1:b_MPPM) ext_bit_6(1:b_MPPM)   ext_bit_7(1:b_MPPM)  ext_bit_8(1:b_MPPM)];
                        %    DCO Estimation
                        out_fft_dco=fft(xe_dco/normalize_coff_1);
                        out_fft_dco_2=fft(re_out_2/normalize_coff_2);
                        
                        out_dco=out_fft_dco(2:N_t/2);
                        out_dco_2=out_fft_dco_2(2:N_half_2);
                        f_out_dco=qamdemod(out_dco,M_qam);
                        f_out_dco_2=qamdemod(out_dco_2,M_qam);
                        bit_out_dco=de2bi(f_out_dco);
                        bit_out_dco_2=de2bi(f_out_dco_2);
                        reshaped_bit_dco=reshape(bit_out_dco',N/2*log2(M_qam),1);
                        reshaped_bit_dco_2=reshape(bit_out_dco_2',N_2*log2(M_qam),1);
                        [ber_3]=biterr(b_3,reshaped_bit_dco_2);
                        [ber_2]=biterr(b_2,reshaped_bit_dco);
                        [ber_1]=biterr(b_1,ext_bit_total');
                        sumerr=sumerr + ber_1 +ber_2;
                        sumerr_2=sumerr_2 + ber_3;
                end
                
                cnt_vect(delta_inx,phi_index,a)=cnt;
                err_rate(delta_inx,phi_index,a)=sumerr/(loop*(bc_1+bc_2));
                err_rate_2(delta_inx,phi_index,a)=sumerr_2/(loop*(bc_3));
                err_mppm(delta_inx,phi_index,a)=ber_1;
                err_dco(delta_inx,phi_index,a)=ber_2;
            end
            %             figure
            %             mean_max_vect(phi_index)=mean(PP_vect);
            %             magnit_max_vect(phi_index)=max(PP_vect);
            %             plot(Ie_t,'r');hold on;plot(I_t,'b');hold on;plot(I_OFF*ones(1,N_t),'g');hold on; plot(I_ON*ones(1,N_t),'g');plot(I_m*ones(1,N_t));
            PAPR_vect=(PP_vect.^2)/E_MPPM_DCO;
            PAPR_vect_2=(PP_vect_2.^2)/E_DCO_2;
            PAPR_vect_dB=10*log10(PAPR_vect);
            PAPR_vect_dB_2=10*log10(PAPR_vect_2);
            length_PAPR_vect=length(PAPR_vect_dB);
            length_PAPR_vect_2=length(PAPR_vect_dB_2);
            in_dB=0;
            for dB=0:0.5:15
                in_dB=in_dB+1;
                PAPR_CCDF(in_dB,phi_index,delta_inx)=(length(find(PAPR_vect_dB > dB)))/length_PAPR_vect;
                PAPR_CCDF_2(in_dB,phi_index,delta_inx)=(length(find(PAPR_vect_dB_2 > dB)))/length_PAPR_vect_2;
            end
        end
        phi_index=phi_index+1;
    end
end

xtime_c=SNR;

for i=1:length(phi_plot)
figure
semilogy(SNR,squeeze(err_rate(1,i,:)),'b--');hold on;semilogy(SNR,squeeze(err_rate(2,i,:)),'r--');hold on;semilogy(SNR,squeeze(err_rate(3,i,:)),'k--')
hold on;semilogy(SNR,squeeze(err_rate(4,i,:)),'g--')
% ;hold on;semilogy(SNR,squeeze(err_rate(5,i,:)),'c--')
hold on;semilogy(SNR,squeeze(err_rate_2(1,i,:)),'bv-');hold on;semilogy(SNR,squeeze(err_rate_2(2,i,:)),'rv-');hold on;semilogy(SNR,squeeze(err_rate_2(3,i,:)),'kv-')
hold on;semilogy(SNR,squeeze(err_rate_2(4,i,:)),'gv-')
% ;hold on;semilogy(SNR,squeeze(err_rate_2(5,i,:)),'cv-')
title([ 'phi' '=' num2str(phi_plot(i)) ])
end

for j=1:length(delta_plot)
figure
semilogy(SNR,squeeze(err_rate(j,1,:)),'b--');hold on;semilogy(SNR,squeeze(err_rate(j,2,:)),'r--');hold on;semilogy(SNR,squeeze(err_rate(j,3,:)),'k--')
hold on;semilogy(SNR,squeeze(err_rate(j,4,:)),'g--');
% hold on;semilogy(SNR,squeeze(err_rate(j,5,:)),'c--')
hold on;semilogy(SNR,squeeze(err_rate_2(j,1,:)),'bv-');hold on;semilogy(SNR,squeeze(err_rate_2(j,2,:)),'rv-');hold on;semilogy(SNR,squeeze(err_rate_2(j,3,:)),'kv-')
hold on;semilogy(SNR,squeeze(err_rate_2(j,4,:)),'gv-');
% hold on;semilogy(SNR,squeeze(err_rate_2(j,5,:)),'cv-')
title([ 'delta' '=' num2str(delta_plot(j))])
end
dB=[0:0.5:15];
for j=1:length(phi_plot)
figure
semilogy(dB,squeeze(PAPR_CCDF(:,j,1)),'b--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,2)),'r--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,3)),'k--')
hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,4)),'g--');
% hold on;semilogy(dB,squeeze(PAPR_CCDF(:,j,5)),'c--')
hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,j,1)),'bv-');hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,j,2)),'rv-');hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,j,3)),'kv-')
hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,j,4)),'gv-');
% hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,j,5)),'cv-')
title([ 'phi' '=' num2str(phi_plot(j))])
end

for j=1:length(delta_plot)
figure
semilogy(dB,squeeze(PAPR_CCDF(:,1,j)),'b--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,2,j)),'r--');hold on;semilogy(dB,squeeze(PAPR_CCDF(:,3,j)),'k--')
hold on;semilogy(dB,squeeze(PAPR_CCDF(:,4,j)),'g--')
% ;hold on;semilogy(dB,squeeze(PAPR_CCDF(:,5,j)),'c--');hold on;
hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,1,j)),'bv-');hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,2,j)),'rv-');hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,3,j)),'kv-')
hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,4,j)),'gv-')
% ;hold on;semilogy(dB,squeeze(PAPR_CCDF_2(:,5,j)),'cv-');hold on;
title([ 'delta' '=' num2str(delta_plot(j))  ])
end

for j=1:length(SNR)
figure
semilogy(phi_plot,squeeze(err_rate_2(1,:,j)),'bv-');hold on;semilogy(phi_plot,squeeze(err_rate_2(2,:,j)),'rv-');hold on;semilogy(phi_plot,squeeze(err_rate_2(3,:,j)),'kv-')
hold on;semilogy(phi_plot,squeeze(err_rate_2(4,:,j)),'gv-')
% ;hold on;semilogy(phi_plot,squeeze(err_rate_2(5,:,j)),'cv-')
title([ 'SNR' '=' num2str(SNR(j))  ])
end
