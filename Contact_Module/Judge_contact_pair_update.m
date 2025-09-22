%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Judge_contact_pair_update(Contact_pair,Contact_condition,positions,tol)
coor=positions.';coor=coor(:);
for ci=1:1:numel(Contact_pair)
    if ~Contact_condition(ci).large_sliding
        Contact_pair(ci).update=false;
        continue
    end
    contact_matrix=Contact_pair(ci).contact_matrix;
    R_matrix=Contact_pair(ci).R_matrix;
    area_vector=Contact_pair(ci).area_vector;
    gap=R_matrix*contact_matrix*coor;
    gapabs=abs(gap(1:3:end))+abs(gap(2:3:end));
    resd=norm(gapabs)/norm(sqrt(area_vector));
    Contact_pair(ci).update=resd>tol;
end
end