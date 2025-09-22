%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function MID=Get_particle_material_ID(Material,Element_set,Mesh)
elements=Mesh.elements;
Nnum=size(Mesh.nodes,1);
Nactivation=Mesh.nactivation;
MID=zeros(Nnum,1);
for mi=1:1:numel(Material)
    setn=Material(mi).set(:);
    GN=Merge_cell(Element_set,setn);
    nid=reshape(elements(GN,:),[],1);
    nid(isnan(nid))=[];
    MID(nid,1)=mi;
end
if ismember(0,MID(Nactivation))
    error('Some nodes are not assigned material properties')
end
end