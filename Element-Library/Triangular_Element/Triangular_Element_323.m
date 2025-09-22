%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Tri3D=Triangular_Element_323(elements,nodes)
Tri3D.type=323;
Tri3D.R_matrix=zeros(6,9,1,size(elements,1));
Tri3D.B_matrix=zeros(4,9,1,size(elements,1));
Tri3D.FB_matrix=zeros(2,3,1,size(elements,1));
Tri3D.detJac=zeros(1,1,1,size(elements,1));
Tri3D.dof=zeros(9,size(elements,1));
B_element=zeros(3,6);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:3);
    dof0=nid*3;dof=[dof0-2;dof0-1;dof0];
    element_nodes3D=nodes(nid,1:3);
    %----------------------------------------------------------------------
    Normal_vector=cross(element_nodes3D(2,:)-element_nodes3D(1,:),element_nodes3D(3,:)-element_nodes3D(1,:));
    R=Normal_to_tangential(Normal_vector/norm(Normal_vector));
    element_nodes=element_nodes3D*R.';
    %----------------------------------------------------------------------
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),2,1);
    dN=Jac\[-1,1,0;-1,0,1];
    B_element([1,7,13;5,11,17])=dN;
    B_element([6,12,18;3,9,15])=dN;
    %----------------------------------------------------------------------
    blank=zeros(2,3);Rs=[R,blank,blank;blank,R,blank;blank,blank,R];
    Tri3D.detJac(1,1,1,ei)=det(Jac)/2;
    if Tri3D.detJac(1,1,1,ei)<0;error('negative area');end
    Tri3D.R_matrix(:,:,1,ei)=Rs;
    Tri3D.B_matrix([1,2,4],:,1,ei)=B_element*Rs;
    Tri3D.FB_matrix([1,2],:,1,ei)=[B_element(1,1:2:end);B_element(2,2:2:end)];
    Tri3D.dof(:,ei)=dof(:);
end
end