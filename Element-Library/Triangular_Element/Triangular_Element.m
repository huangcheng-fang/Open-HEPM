%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Tri=Triangular_Element(elements,nodes)
Tri.type=223;
Tri.B_matrix=zeros(3,6,1,size(elements,1));
Tri.detJac=zeros(1,1,1,size(elements,1));
Tri.dof=zeros(6,size(elements,1));
B_element=zeros(3,6);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:3);
    dof=[nid*2-1;nid*2];
    element_nodes=nodes(nid,1:2);
    %----------------------------------------------------------------------
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),2,1);
    dN=Jac\[-1,1,0;-1,0,1];
    B_element([1,7,13;5,11,17])=dN;
    B_element([6,12,18;3,9,15])=dN;
    %----------------------------------------------------------------------
    detJacs=det(Jac)/2;
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    Tri.detJac(1,1,1,ei)=detJacs;
    Tri.B_matrix(:,:,1,ei)=B_element;
    Tri.dof(:,ei)=dof(:);
end
end