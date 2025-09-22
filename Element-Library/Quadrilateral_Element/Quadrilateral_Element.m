%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Quad=Quadrilateral_Element(elements,nodes,IntegralOder)
[ip,wt]=Gauss_Integral(2,IntegralOder);
dN=Quadrilateral_interpolation_derivative(ip);
IPNUM=size(ip,1);Enum=size(elements,1);
%--------------------------------------------------------------------------
Quad.type=324;
Quad.R_matrix=zeros(8,12,1,Enum);
Quad.B_matrix=zeros(4,12,IPNUM,Enum);
Quad.FB_matrix=zeros(2,4,IPNUM,Enum);
Quad.detJac=zeros(1,1,IPNUM,Enum);
Quad.dof=zeros(12,Enum);
Bs_element=zeros(3,8,IPNUM);detJacs=zeros(1,1,IPNUM);
for ei=1:1:Enum
    %----------------------------------------------------------------------
    nid=elements(ei,1:4);
    dof0=nid*3;dof=[dof0-2;dof0-1;dof0];
    element_nodes3D=nodes(nid,1:3);
    %----------------------------------------------------------------------
    Normal_vector=Get_normal_vector(element_nodes3D);
    R=Normal_to_tangential(Normal_vector/norm(Normal_vector));
    element_nodes=element_nodes3D*R.';
    %----------------------------------------------------------------------
    for pii=1:1:IPNUM
        [B_element,detJac]=Form_B_Quad(element_nodes,dN(:,:,pii));
        detJacs(1,1,pii)=wt(pii)*detJac;
        Bs_element(:,:,pii)=B_element;
    end
    %----------------------------------------------------------------------
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    blank=zeros(2,3);Rs=[R,blank,blank,blank;blank,R,blank,blank;blank,blank,R,blank;blank,blank,blank,R];
    Quad.detJac(:,:,:,ei)=detJacs;
    Quad.R_matrix(:,:,1,ei)=Rs;
    Quad.B_matrix([1,2,4],:,:,ei)=pagemtimes(Bs_element,Rs);
    Quad.FB_matrix([1,2],:,:,ei)=[Bs_element(1,1:2:end,:);Bs_element(2,2:2:end,:)];
    Quad.dof(:,ei)=dof(:);
end
end

function [B_element,detJ]=Form_B_Quad(element_nodes,dNi)
J=dNi*element_nodes;
detJ=det(J);
LB=J\dNi;
B_element=[LB(1,1),0,LB(1,2),0,LB(1,3),0,LB(1,4),0;
           0,LB(2,1),0,LB(2,2),0,LB(2,3),0,LB(2,4);
           LB(2,1),LB(1,1),LB(2,2),LB(1,2),LB(2,3),LB(1,3),LB(2,4),LB(1,4);];
end
