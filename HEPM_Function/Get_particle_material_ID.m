%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function PMID=Get_particle_material_ID(Material,Particle_set,Pnum)
PMID=zeros(Pnum,1);
for mi=1:1:numel(Material)
    setn=Material(mi).particle_set(:);
    GN=Merge_cell(Particle_set,setn);
    PMID(GN,1)=mi;
end
if ismember(0,PMID)
    error('Some particles are not assigned material properties')
end
end