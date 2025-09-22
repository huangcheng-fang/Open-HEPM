%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Link=Link_Element_312(elements,nodes,IntegralOder)
[ip,wt]=Gauss_Integral(1,IntegralOder);
dN=Line_interpolation_derivative(ip);
IPNUM=IntegralOder;
%--------------------------------------------------------------------------
Link.type=312;
Link.R_matrix=zeros(2,6,IPNUM,size(elements,1));
Link.B_matrix=zeros(3,6,IPNUM,size(elements,1));
Link.FB_matrix=zeros(1,2,IPNUM,size(elements,1));
Link.detJac=zeros(1,1,IPNUM,size(elements,1));
Link.dof=zeros(6,size(elements,1));
Bs_element=zeros(1,2,IPNUM);detJacs=zeros(1,1,IPNUM);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:2);
    dof0=nid*3;dof=[dof0-2;dof0-1;dof0];
    %----------------------------------------------------------------------
    R=nodes(nid(2),:)-nodes(nid(1),:);
    L=norm(R);R=R/L;
    element_nodes=[0;L];
    %----------------------------------------------------------------------
    for pii=1:1:IPNUM
        [B_element,detJac]=Form_B_Link(element_nodes,dN(:,:,pii));
        detJacs(pii)=wt(pii)*detJac;
        Bs_element(:,:,pii)=B_element;
    end
    %----------------------------------------------------------------------
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    blank=zeros(1,3);Rs=[R,blank;blank,R];
    Link.detJac(:,:,:,ei)=detJacs;
    Link.B_matrix(1,:,:,ei)=Bs_element*Rs;
    Link.FB_matrix(:,:,:,ei)=Bs_element;
    Link.dof(:,ei)=dof(:);
end
end

function [B_element,detJac]=Form_B_Link(element_nodes,dNi)
J=dNi*element_nodes;
B_element=J\dNi;detJac=det(J);
end