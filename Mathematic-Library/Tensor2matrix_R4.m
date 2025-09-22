%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function B=Tensor2matrix_R4(A)
Asize=size(A,[1,2,3,4,5,6]);
B=zeros(Asize(1)*Asize(2),Asize(3)*Asize(4),Asize(5),Asize(6));
% syms B [Asize(1)*Asize(2),Asize(3)*Asize(4)]
for i=1:1:3
    for j=1:1:3
        for m=1:1:3
            for n=1:1:3
                B((i-1)*3+j,(m-1)*3+n,:,:)=A(i,j,m,n,:,:);
            end
        end
    end
end
end