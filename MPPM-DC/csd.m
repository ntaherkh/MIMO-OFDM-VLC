d=4;
N=2^d;
W=exp(-2*i*(pi/N));
a=[0:N-1];
J_K=a'*a;
F_d=W.^(J_K)
% a= [ 2 3 ; 3 4 ];
% A=[a a*2 ; a.^2  a-2];
b=reshape(F_d,2,2,2,2,2,2,2,2);
c=permute(b,[1,5,2,6,3,7,4,8]);
d=reshape(c,4,4,4,4);
E=tt_matrix(b);
F=tt_tensor(d);
E.core'
F.core'


