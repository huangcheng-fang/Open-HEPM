%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%kappa: swell/recompression index
%niu: poisson_ratio
%e: Void ratio
%p: Mean stress
%============================Secant stiffness==============================
function [D_elastic,D_elastic0,e]=Form_D_porous_elastic(kappa,niu,e,p,dstrainv)
dstrainv(dstrainv==0)=1e-6;
dp=p.*(exp(1/kappa*(1+e).*(1-exp(-dstrainv))))-p;
K=dp./dstrainv;
G=3*(1-2*niu)/2/(1+niu)*K;
G2=2*G;H=K-G2/3;F=G2+H;
O=zeros(1,1,size(e,3),size(e,4));
D_elastic0=[F, H, H, O, O, O;
           H, F, H, O, O, O;
           H, H, F, O, O, O;
           O, O, O, G, O, O;
           O, O, O, O, G, O;
           O, O, O, O, O, G;];
e=max(exp(-dstrainv).*(1+e)-1,0);
K=1/kappa*(1+e).*(p+dp);
G=3*(1-2*niu)/2/(1+niu)*K;
G2=2*G;H=K-G2/3;F=G2+H;
O=zeros(1,1,size(e,3),size(e,4));
D_elastic=[F, H, H, O, O, O;
           H, F, H, O, O, O;
           H, H, F, O, O, O;
           O, O, O, G, O, O;
           O, O, O, O, G, O;
           O, O, O, O, O, G;];
% loc=dp<0;
% D_elastic(:,:,loc)=D_elastic0(:,:,loc);
D_elastic=D_elastic0;
end
%============================Tangent stiffness==============================
% function [D_elastic,e]=Form_D_porous_elastic(kappa,niu,e,p,dstrainv)
% K=1/kappa*(1+e).*p;
% G=3*(1-2*niu)/2/(1+niu)*K;
% G2=2*G;H=K-G2/3;F=G2+H;
% O=zeros(1,1,size(e,3),size(e,4));
% D_elastic=[F, H, H, O, O, O;
%            H, F, H, O, O, O;
%            H, H, F, O, O, O;
%            O, O, O, G, O, O;
%            O, O, O, O, G, O;
%            O, O, O, O, O, G;];
% e=max(exp(-dstrainv).*(1+e)-1,0);
% end