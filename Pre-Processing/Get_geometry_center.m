%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function center=Get_geometry_center(elements,nodes)
center=zeros(size(elements,1),size(nodes,2));
num=zeros(size(elements,1),1);
for i=1:1:size(elements,2)
    nid=elements(:,i);
    loc=~isnan(nid);
    center(loc,:)=center(loc,:)+nodes(nid(loc),:);
    num(loc)=num(loc)+1;
end
center=center./num;
end