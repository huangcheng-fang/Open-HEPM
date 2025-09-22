%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%Inner product of two rank-4 tensors
%Cijmn=Aijkl*Bklmn ---> summation convention for repeated subscripts
function C=Inner_R4R4(A,B)
Asize=size(A,[1,2,3,4]);
Bsize=size(B,[1,2,3,4]);
C=zeros([Asize([1,2]),Bsize([3,4])]);
% syms C [Asize(1),Asize(2),Bsize(3),Bsize(4)]
for i=1:1:Asize(1)
    for j=1:1:Asize(2)
        for m=1:1:Bsize(3)
            for n=1:1:Bsize(4)
                a=A(i,j,:,:);
                b=B(:,:,m,n);
                C(i,j,m,n)=a(:).'*b(:);
            end
        end
    end
end
end