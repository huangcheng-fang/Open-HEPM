%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Mesh=Update_surface_geometry(Mesh)
nodes=Mesh.nodes;
surfaces=Mesh.surfaces;
dim=size(nodes,2);
snum=size(surfaces,1);
s_n_num=zeros(snum,1);
scenter=zeros(snum,dim);
%==========================================================================
for j=1:1:size(surfaces,2)
    nid=surfaces(:,j);
    sid=find(nid~=0&~isnan(nid));
    for k=1:1:dim
        scenter(:,k)=scenter(:,k)+accumarray(sid,nodes(nid(sid),k),[snum,1]);
    end
    s_n_num=s_n_num+accumarray(sid,1,[snum,1]);
end
scenter=s_n_num.\scenter;
%===============================Output=====================================
Mesh.scenter=scenter;
end