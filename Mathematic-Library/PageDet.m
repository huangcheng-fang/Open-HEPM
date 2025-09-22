%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function B=PageDet(A)
asize=size(A);
num=numel(A)/(asize(1)*asize(2));
B=zeros(1,1,num);
for ii=1:1:num
    B(1,1,ii)=det(A(:,:,ii));
end
asize(1:2)=[1,1];
B=reshape(B,asize);
end