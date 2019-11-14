clear
a= [ 2 3 4 5; 3 -3 5 3; 3 4 21 2; 2 3 0 5];
A=[a a*2 ; a.^2  a-2];
b=reshape(A,4,4,4);
n_k=4;
E=tt_tensor(b)

b1=reshape(b,4,16);
[U1,S1,V1]=svd(b1);


M1=S1(1:4,1:4)*V1(:,1:4)';
b2=reshape(M1,16, 4)
[u2,S2,V2]=svd(b2);
U22=u2(:,1:4);
U23=reshape(U22,4,4,4);
U2=permute(U23,[2,1,3]);
M2=S2(1:4,1:4)*V2';


G1=U1(:,1:4);
G2=U2;
G3=M2'; 

for i=1:4
    for j=1:4
        for k=1:4
         g2=squeeze(G2(j,:,:));
         EST(i,j,k)=G1(i,:)*g2*G3(k,:)';
        end
    end
end

MAT=reshape(EST,8,8);
         
vect1=E.core(1:n_k*E.r(2));
H1=reshape(vect1,n_k,E.r(2));

vect2=E.core(n_k*E.r(2)+1:n_k*E.r(2)*E.r(3)+n_k*E.r(2));
H2=reshape(vect2,E.r(2),n_k,E.r(3));
H2=permute(H2,[2,1,3]);
vect3=E.core(end-n_k*E.r(3)+1:end);
H3=reshape(vect3,E.r(3),n_k);
H3=H3';

for i=1:4
    for j=1:4
        for k=1:4
         h2=squeeze(H2(j,:,:));
         EST_2(i,j,k)=H1(i,:)*h2*H3(k,:)';
        end
    end
end

MAT_2=reshape(EST_2,8,8)
MAT
A
