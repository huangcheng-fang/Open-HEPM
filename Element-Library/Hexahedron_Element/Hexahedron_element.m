%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Hex=Hexahedron_element(elements,nodes,IntegralOder)
[ip,wt]=Gauss_Integral(3,IntegralOder);
N_matrix=Hexahedron_interpolation(ip);
dN=Hexahedron_interpolation_derivative(ip);
IPNUM=IntegralOder^3;
%--------------------------------------------------------------------------
Hex.R_matrix=zeros(0,0,1,size(elements,1));
Hex.B_matrix=zeros(6,24,IPNUM,size(elements,1));
Hex.DetJac=zeros(1,1,IPNUM,size(elements,1));
Hex.Nodes=zeros(8,3,1,size(elements,1));
Hex.Elements=permute(elements,[4,2,3,1]);
Hex.N_matrix=permute(N_matrix,[3,2,1]);
Bs_element=zeros(6,24,IPNUM);detJacs=zeros(1,1,IPNUM);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:8);
    element_nodes=nodes(nid,:);
    %----------------------------------------------------------------------
    for pii=1:1:IPNUM
        [B_element,detJac]=Form_B_Hex(element_nodes,dN(:,:,pii));
        detJacs(pii)=wt(pii)*detJac;
        Bs_element(:,:,pii)=B_element;
    end
    %----------------------------------------------------------------------
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    Hex.DetJac(:,:,:,ei)=detJacs;
    Hex.B_matrix(:,:,:,ei)=Bs_element;
    Hex.Nodes(:,:,1,ei)=element_nodes;
end
end

function [B_element,detJac]=Form_B_Hex(element_nodes,dNi)
J=dNi*element_nodes;
LB=J\dNi;detJac=det(J);
B_element=[LB(1,1),0,0,LB(1,2),0,0,LB(1,3),0,0,LB(1,4),0,0,LB(1,5),0,0,LB(1,6),0,0,LB(1,7),0,0,LB(1,8),0,0;
           0,LB(2,1),0,0,LB(2,2),0,0,LB(2,3),0,0,LB(2,4),0,0,LB(2,5),0,0,LB(2,6),0,0,LB(2,7),0,0,LB(2,8),0;
           0,0,LB(3,1),0,0,LB(3,2),0,0,LB(3,3),0,0,LB(3,4),0,0,LB(3,5),0,0,LB(3,6),0,0,LB(3,7),0,0,LB(3,8);
           LB(2,1),LB(1,1),0,LB(2,2),LB(1,2),0,LB(2,3),LB(1,3),0,LB(2,4),LB(1,4),0,LB(2,5),LB(1,5),0,LB(2,6),LB(1,6),0,LB(2,7),LB(1,7),0,LB(2,8),LB(1,8),0;
           0,LB(3,1),LB(2,1),0,LB(3,2),LB(2,2),0,LB(3,3),LB(2,3),0,LB(3,4),LB(2,4),0,LB(3,5),LB(2,5),0,LB(3,6),LB(2,6),0,LB(3,7),LB(2,7),0,LB(3,8),LB(2,8);
           LB(3,1),0,LB(1,1),LB(3,2),0,LB(1,2),LB(3,3),0,LB(1,3),LB(3,4),0,LB(1,4),LB(3,5),0,LB(1,5),LB(3,6),0,LB(1,6),LB(3,7),0,LB(1,7),LB(3,8),0,LB(1,8);];
end