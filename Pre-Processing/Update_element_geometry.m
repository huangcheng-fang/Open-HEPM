%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Mesh=Update_element_geometry(Mesh)
nodes=Mesh.nodes;
elements=Mesh.elements;
enum=size(elements,1);
dim=size(nodes,2);
ecenter=zeros(enum,dim);
e_n_num=zeros(enum,1);
%==========================================================================
for j=1:1:size(elements,2)
    nid=elements(:,j);
    eid=find(nid~=0&~isnan(nid));
    for k=1:1:dim
        ecenter(:,k)=ecenter(:,k)+accumarray(eid,nodes(nid(eid),k),[enum,1]);
    end
    e_n_num=e_n_num+accumarray(eid,1,[enum,1]);
end
ecenter=e_n_num.\ecenter;
%==============================Output======================================
Mesh.ecenter=ecenter;
end