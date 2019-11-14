clear
clc
loop=20000;

%% QAM Modulation order 
loog_M=6;
%  4-QAM , N=130, HAZV=0,SISHMA=0.1231
%% Full subcarrier size
N=62;
hazv=2;
 0.1969
 0.1965
 0.1961
for jj=1:loop
jj
b_2=randi([0 1],loog_M*(N/2-1),1);
c_2=reshape(b_2,(N/2-1),loog_M);
d_2=bi2de(c_2);
in_2=qammod(d_2,2^loog_M);
out_en=zeros(1,N);
% out_en(2:64-numPuncs)=in;
out_en_2(2:N/2-hazv) =in_2(1:N/2-1-hazv);
out_en_2(1)=0;
out_en_2(N/2+1)=0;
for i=N/2+1:N
    out_en_2(i)=conj(out_en_2(N+2-i));
end
g_2=ifft(out_en_2);
l_1=(jj-1)*N+1
l_2=N*jj;
g_vect(l_1:l_2)=g_2;
end
meaned_vect=g_vect-mean(g_vect);
% mean(g_vect);
varr= mean(meaned_vect.^2);
sigma=sqrt(varr)
max(g_vect)