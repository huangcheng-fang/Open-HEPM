%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function dN=Quadrilateral_interpolation_derivative(ip)
% auxX=[1,-1,-1,1];
% auxY=[1,1,-1,-1];
% dN=zeros(2,4,size(ip,1));
% dN(2,:,:)=(0.25*(auxY+ip(:,1)*(auxY.*auxX)))';
% dN(1,:,:)=(0.25*(auxX+ip(:,2)*(auxY.*auxX)))';
%等效如下
x=ip(:,1)';y=ip(:,2)';
dN=zeros(2,4,size(ip,1));
dN(2,:,:)=repmat([0.25;0.25;-0.25;-0.25],1,numel(x))+[0.25;-0.25;0.25;-0.25]*x;
dN(1,:,:)=repmat([0.25;-0.25;-0.25;0.25],1,numel(x))+[0.25;-0.25;0.25;-0.25]*y;
end