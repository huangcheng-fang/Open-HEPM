%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [dstrain,spin]=Triangular_large_strain(elements,nodes,dU)
nodeU=reshape(dU,2,[]).';
nodes=nodes-nodeU;
%--------------------------------------------------------------
dstrain=zeros(4,size(elements,1));
spin=zeros(1,size(elements,1));
for ei=1:1:size(elements,1)
    %----------------------------------------------------------
    nid=elements(ei,1:3);
    element_coor=nodes(nid,1:2);
    element_u=nodeU(nid,1:2);
    %==========================================================
    Jac=element_coor(2:end,:)-repmat(element_coor(1,:),2,1);
    dJac=element_u(2:end,:)-repmat(element_u(1,:),2,1);
    Nint=min(ceil(max(sum(abs(dJac),2)./sum(abs(Jac),2))/0.0001+1e-32),1000);
    ddJac=dJac/Nint;
    ddUdx=[-1,1,0;-1,0,1]*element_u/Nint;
    dEstrain=zeros(4,1);
    Espin=0;
    for nti=1:1:Nint
        Jac=Jac+ddJac;
        ddUdX=Jac\ddUdx;
        Espin=Espin+(ddUdX(1,2)-ddUdX(2,1));
        dEstrain=dEstrain+[ddUdX(1,1)+0.5*(ddUdX(1,1)^2+ddUdX(1,2)^2);
            ddUdX(2,2)+0.5*(ddUdX(2,1)^2+ddUdX(2,2)^2);
            0;
            (ddUdX(2,1)+ddUdX(1,2)+ddUdX(1,1)*ddUdX(2,1)+ddUdX(1,2)*ddUdX(2,2));];
    end
    dstrain(:,ei)=dEstrain;
    spin(ei)=Espin*0.5;
end
end