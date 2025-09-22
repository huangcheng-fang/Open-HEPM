%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Tet=Tetrahedron_quadratic_element(elements,nodes,IntegralOrder)
[ip,wt]=Hammer_integral(3,IntegralOrder);
N_matrix=Tetrahedron_quadratic_interpolation(ip);
dN=Tetrahedron_quadratic_interpolation_derivative(ip);
IPNUM=numel(wt);wt=permute(wt,[3,2,1]);
%--------------------------------------------------------------------------
Tet.R_matrix=zeros(0,0,IPNUM,size(elements,1));
Tet.B_matrix=zeros(6,30,IPNUM,size(elements,1));
Tet.DetJac=zeros(1,1,IPNUM,size(elements,1));
Tet.Nodes=zeros(10,3,1,size(elements,1));
Tet.Elements=permute(elements,[4,2,3,1]);
Tet.N_matrix=permute(N_matrix,[3,2,1]);
Bs_element=zeros(6,30,IPNUM);detJacs=zeros(1,1,IPNUM);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:10);
    element_nodes=nodes(nid,:);
    %----------------------------------------------------------------------
    for pii=1:1:IPNUM
        [B_element,detJac]=Form_B_Tet_quad(element_nodes,dN(:,:,pii));
        detJacs(pii)=wt(pii)*detJac;
        Bs_element(:,:,pii)=B_element;
    end
    %----------------------------------------------------------------------
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    Tet.DetJac(:,:,:,ei)=detJacs;
    Tet.B_matrix(:,:,:,ei)=Bs_element;
    Tet.Nodes(:,:,1,ei)=element_nodes;
end
end

function [B_element,detJac]=Form_B_Tet_quad(element_nodes,dNi)
J=dNi*element_nodes;
LB=J\dNi;detJac=det(J);
B_element=zeros(6,3*size(element_nodes,1));
B_element(1,1:3:end)=LB(1,:);
B_element(2,2:3:end)=LB(2,:);
B_element(3,3:3:end)=LB(3,:);
B_element(4,1:3:end)=LB(2,:);
B_element(4,2:3:end)=LB(1,:);
B_element(5,2:3:end)=LB(3,:);
B_element(5,3:3:end)=LB(2,:);
B_element(6,1:3:end)=LB(3,:);
B_element(6,3:3:end)=LB(1,:);
end