clear
d=3;
N=2^d;
W=exp(-2*i*(pi/N));
a=[0:N-1];
J_K=a'*a;
F_d=W.^(J_K);
n_k=4;
% a= [ 2 3 ; 3 4 ];
% A=[a a*2 ; a.^2  a-2];
b=reshape(F_d,2,2,2,2,2,2);
% C_d=tt_matrix(b);
% full(C_d)
c=permute(b,[1,4,2,5,3,6]);
d=reshape(c,4,4,4);
z=reshape(d,size(b));
zb=permute(z,[1,3,5,2,4,6]);
zf=reshape(zb,size(F_d))

% E=tt_matrix(b);
E=tt_tensor(d);

l1=n_k*E.r(2);
l2=n_k*E.r(2)*E.r(3)+l1;
l3=n_k*E.r(3)*E.r(4)+l2;



vect1=E.core(1:l1);
H1=reshape(vect1,E.r(2),n_k);
% H1=permute(H1,[2,1]);

vect2=E.core(l1+1:l2);
H2=reshape(vect2,E.r(2),n_k,E.r(3));
H2=permute(H2,[2,1,3]);

vect3=E.core(l2+1:l3);
H3=reshape(vect3,E.r(3),n_k,E.r(4));
H3=conj(H3');

for i=1:4
    for j=1:4
        for k=1:4
            h2=squeeze(H2(j,:,:));
            EST_2(i,j,k)=H1(i,:)*h2*(conj(H3(k,:)'));
            
        end
    end
end
Z=reshape(EST_2,size(b));
ZB=permute(Z,[1,3,5,2,4,6]);
MAT_2=reshape(ZB,size(F_d))
% EST_2
% SS1=permute(EST_2,[3,2,1])
% MAT_2=reshape(SS1,8,8);
F_d
MAT_2