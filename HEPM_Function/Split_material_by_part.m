%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [newMaterial,Element_set,Particle_set,CPart_list]=Split_material_by_part(Material,Element_set,Particle_set,Einpart,Pinpart)
coder.varsize('eset_id')
eset_id=coder.nullcopy(cell(0,1));
coder.varsize('pset_id')
pset_id=coder.nullcopy(cell(0,1));
coder.varsize('m_id')
m_id=coder.nullcopy(cell(0,1));
coder.varsize('part_list')
CPart_list=coder.nullcopy(cell(0,1));
for mi=1:1:numel(Material)
    esetn=Material(mi).set(:);
    eid=Merge_cell(Element_set,esetn);
    psetn=Material(mi).particle_set(:);
    pid=Merge_cell(Particle_set,psetn);
    einpart=Einpart(eid);
    pinpart=Pinpart(pid);
    part_list1=unique(einpart);
    part_list2=unique(pinpart);
    part_list=unique([part_list1;part_list2]);
    for i=1:1:numel(part_list)
        Element_set{end+1,1}=eid(einpart==part_list(i));
        Particle_set{end+1,1}=pid(pinpart==part_list(i));
        eset_id{end+1,1}=size(Element_set,1);
        pset_id{end+1,1}=size(Particle_set,1);
        m_id{end+1,1}=mi;
        CPart_list{end+1,1}=part_list(i);
    end
end
%==========================================================================
Materialcell=struct2cell(Material);
m_id=Montage_cell(m_id,1);
newMaterial=struct('set',eset_id,'particle_set',pset_id,'elasticity',Materialcell(3,m_id).','section',Materialcell(4,m_id).','plasticity',Materialcell(5,m_id).');
CPart_list=Montage_cell(CPart_list,1);
end