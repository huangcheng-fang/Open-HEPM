%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%Dyadic product of two rank-2 tensors
%Cijkl=Aij*Bkl
function C=Dyadic_R2R2(A,B,permution)
Asize=size(A,[1,2]);
Bsize=size(B,[1,2,3,4]);
C=zeros([Asize,Bsize]);
% syms C [Asize,Bsize]
%==========================================================================
if isequal(A,eye(3))
    for i=1:1:Asize(1)
        C(i,i,:,:,:,:)=B;
    end
else
    for i=1:1:Asize(1)
        for j=1:1:Asize(2)
            C(i,j,:,:,:,:)=pagemtimes(A(i,j,:,:),B);
        end
    end
end
%==========================================================================
if ~isempty(permution)
    permution=[permution,5,6];
    C=permute(C,permution);
end
end