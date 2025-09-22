%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [dstrain,spin]=Tetrahedron_large_strain(Cnodes,CNodes_stepU)
[ipE,wtE]=Gauss_Integral(1,1);
ipE=(ipE+1)/2;wtE=wtE/2;
IPNUM=1;
ENUM=size(Cnodes,4);
%--------------------------------------------------------------
dstrain=zeros(6,1,IPNUM,ENUM);
spin=zeros(3,3,IPNUM,ENUM);
for ei=1:1:ENUM
    %----------------------------------------------------------
    element_u=CNodes_stepU(:,:,:,ei);
    element_coor=Cnodes(:,:,:,ei)-element_u;
    %==========================================================
    Jac=element_coor(2:end,:)-repmat(element_coor(1,:),3,1);
    dJac=element_u(2:end,:)-repmat(element_u(1,:),3,1);
    %==========================================================
    dUdX=0;
    for kk=1:1:numel(ipE)
        dUdX=dUdX+(Jac+ipE(kk)*dJac)\dJac*wtE(kk);
    end
    dEstrain=[dUdX(1,1);dUdX(2,2);dUdX(3,3);dUdX(2,1)+dUdX(1,2);dUdX(3,2)+dUdX(2,3);dUdX(1,3)+dUdX(3,1)];
    %==========================================================
    dUdX=(Jac+0.5*dJac)\dJac;
    dw=[dUdX(2,1)-dUdX(1,2);dUdX(3,2)-dUdX(2,3);dUdX(1,3)-dUdX(3,1)]*0.5;
    Espin=[0,dw(1),-dw(3);-dw(1),0,dw(2);dw(3),-dw(2),0];
    % Espin=eye(3)+(eye(3)-0.5*Espin)\Espin; 
    %====================Output================================
    dstrain(:,:,:,ei)=dEstrain;
    %spin(:,:,:,ei)=Espin;
    spin(:,:,:,ei)=dUdX.';
end
end


% w=[dUdX(3,2)-dUdX(2,3),dUdX(1,3)-dUdX(3,1),dUdX(2,1)-dUdX(1,2)]*0.5;
% normw=norm(w);theta=normw;w=w/norm(w);
% if normw<1e-12;w=[1,0,0];end
% Espin=[w(1)^2*(1-cos(theta))+cos(theta),w(1)*w(2)*(1-cos(theta))+w(3)*sin(theta),w(1)*w(3)*(1-cos(theta))-w(2)*sin(theta);
%        w(1)*w(2)*(1-cos(theta))-w(3)*sin(theta),w(2)*w(2)*(1-cos(theta))+cos(theta),w(2)*w(3)*(1-cos(theta))+w(1)*sin(theta)
%        w(1)*w(3)*(1-cos(theta))+w(2)*sin(theta),w(2)*w(3)*(1-cos(theta))-w(1)*sin(theta),w(3)*w(3)*(1-cos(theta))+cos(theta)];