%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Tet=Tetrahedron_element(elements,nodes,IntegralOrder)
[ip,wt]=Hammer_integral(3,IntegralOrder);
N_matrix=Tetrahedron_interpolation(ip);
IPNUM=numel(wt);wt=permute(wt,[3,2,1]);
%--------------------------------------------------------------------------
Tet.R_matrix=zeros(0,0,IPNUM,size(elements,1));
Tet.B_matrix=zeros(6,12,IPNUM,size(elements,1));
Tet.DetJac=zeros(1,1,IPNUM,size(elements,1));
Tet.Nodes=zeros(4,3,1,size(elements,1));
Tet.Elements=permute(elements,[4,2,3,1]);
Tet.N_matrix=permute(N_matrix,[3,2,1]);
B_element=zeros(6,12);
for ei=1:1:size(elements,1)
    %----------------------------------------------------------------------
    nid=elements(ei,1:4);
    element_nodes=nodes(nid,:);
    %----------------------------------------------------------------------
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),3,1);
    dN=Jac\[-1,1,0,0;-1,0,1,0;-1,0,0,1];
    B_element([1,19,37,55;8,26,44,62;15,33,51,69])=dN;
    B_element([10,28,46,64;4,22,40,58;11,29,47,65])=dN;
    B_element([18,36,54,72;17,35,53,71;6,24,42,60])=dN;
    %----------------------------------------------------------------------
    detJacs=det(Jac);
    if any(detJacs<0);error('negative area');end
    %----------------------------------------------------------------------
    Tet.DetJac(1,1,IPNUM,ei)=detJacs*wt;
    Tet.B_matrix(:,:,IPNUM,ei)=repmat(B_element,1,1,IPNUM);
    Tet.Nodes(:,:,1,ei)=element_nodes;
end
end