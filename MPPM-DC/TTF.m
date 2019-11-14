clear
d=4;
N=2^d;
W=exp(-2*i*(pi/N));
a=[0:N-1];
J_K=a'*a;
F_d=W.^(J_K);
n_k=4;
% a= [ 2 3 ; 3 4 ];
% A=[a a*2 ; a.^2  a-2];
b=reshape(F_d,2,2,2,2,2,2,2,2);
% C_d=tt_matrix(b);
% full(C_d)
c=permute(b,[1,5,2,6,3,7,4,8]);
d=reshape(c,4,4,4,4);
% E=tt_matrix(b);
E=tt_tensor(d);

l1=n_k*E.r(2);
l2=n_k*E.r(2)*E.r(3)+l1;
l3=n_k*E.r(3)*E.r(4)+l2;
l4=n_k*E.r(4)*E.r(5)+l3;


vect1=E.core(1:l1);
H1=reshape(vect1,E.r(2),n_k);
% H1=permute(H1,[2,1]);

vect2=E.core(l1+1:l2);
H2=reshape(vect2,E.r(2),n_k,E.r(3));
H2=permute(H2,[2,1,3]);

vect3=E.core(l2+1:l3);
H3=reshape(vect3,E.r(3),n_k,E.r(4));
H3=permute(H3,[2,1,3]);

vect4=E.core(l3+1:l4);
H4=reshape(vect4,E.r(4),n_k);
H4=conj(H4');

for i=1:4
    for j=1:4
        for k=1:4
            for m=1:4
                h2=squeeze(H2(j,:,:));
                h3=squeeze(H3(k,:,:));
                EST_2(i,j,k,m)=H1(i,:)*h2*h3*(conj(H4(m,:)'));
            end
        end
    end
end
Z=reshape(EST_2,size(b));
ZB=permute(Z,[1,3,5,7,2,4,6,8]);
MAT_2=reshape(ZB,size(F_d))
% EST_2
% SS1=permute(EST_2,[2,1,3,4])
% MAT_2=reshape(SS1,16,16);
F_d
MAT_2

sum(sum(abs(F_d-MAT_2)))