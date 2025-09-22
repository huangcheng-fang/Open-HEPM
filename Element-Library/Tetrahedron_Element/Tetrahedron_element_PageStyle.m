%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Tet=Tetrahedron_element_PageStyle(elements,nodes)
Tet.B_matrix=zeros(6,12,size(elements,1));
Tet.detJac=zeros(size(elements,1),1);
Tet.dof=zeros(size(elements,1),12);
%--------------------------------------------------------------------------
element_nodes=pagetranspose(reshape(nodes(elements(:,1:4).',:).',3,4,[]));
Jac=element_nodes(2:end,:,:)-repmat(element_nodes(1,:,:),3,1,1);
[invJ,detJ]=Page_inv_3D(Jac);
dN=pagemtimes(invJ,[-1,1,0,0;-1,0,1,0;-1,0,0,1]);
% dN=pagemldivide(Jac,[-1,1,0,0;-1,0,1,0;-1,0,0,1]);
ndof=elements(:,1:4)*3;
%--------------------------------------------------------------------------
Tet.B_matrix(1,1:3:end,:)=dN(1,:,:);
Tet.B_matrix(2,2:3:end,:)=dN(2,:,:);
Tet.B_matrix(3,3:3:end,:)=dN(3,:,:);

Tet.B_matrix(4,1:3:end,:)=dN(2,:,:);
Tet.B_matrix(4,2:3:end,:)=dN(1,:,:);

Tet.B_matrix(5,2:3:end,:)=dN(3,:,:);
Tet.B_matrix(5,3:3:end,:)=dN(2,:,:);

Tet.B_matrix(6,1:3:end,:)=dN(3,:,:);
Tet.B_matrix(6,3:3:end,:)=dN(1,:,:);
%--------------------------------------------------------------------------
Tet.detJac=detJ/6;
Tet.dof(:,[1:3:end,2:3:end,3:3:end])=[ndof-2,ndof-1,ndof];
end