%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function MID=Get_element_material_ID(Material,Element_set,Enum)
MID=zeros(Enum,1);
for mi=1:1:numel(Material)
    setn=Material(mi).set(:);
    GN=Merge_cell(Element_set,setn);
    MID(GN,1)=mi;
end
if ismember(0,MID)
    error('Some elements are not assigned material properties')
end
end