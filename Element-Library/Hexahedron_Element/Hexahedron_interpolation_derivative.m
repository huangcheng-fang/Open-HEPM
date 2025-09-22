%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function dN=Hexahedron_interpolation_derivative(ip)
X=ip(:,1);Y=ip(:,2);Z=ip(:,3);
auxX=[-1,-1,-1,-1,1,1,1,1];
auxY=[1,-1,-1,1,1,-1,-1,1];
auxZ=[1,1,-1,-1,1,1,-1,-1];
dN=zeros(3,8,size(ip,1));
dN(3,:,:)=(0.125*(1+X*auxX).*(repmat(auxZ,numel(Y),1)+Y*(auxY.*auxZ)))';
dN(2,:,:)=(0.125*(1+X*auxX).*(repmat(auxY,numel(Z),1)+Z*(auxZ.*auxY)))';
dN(1,:,:)=(0.125*(1+Y*auxY).*(repmat(auxX,numel(Z),1)+Z*(auxZ.*auxX)))';
end