%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function quality=Tet_element_quality(tet_elements,nodes)
quality=zeros(size(tet_elements,1),1);
invW=inv([1,0.5,0.5;0,sqrt(3)/2,sqrt(3)/6;0,0,sqrt(2/3)]);
for ei=1:1:size(tet_elements,1)
    nid=tet_elements(ei,1:4);
    element_nodes=nodes(nid,1:3);
    Jac=element_nodes(2:end,:)-repmat(element_nodes(1,:),3,1);
    a(1)=norm(cross(Jac(1,:),Jac(2,:)));
    a(2)=norm(cross(Jac(2,:),Jac(3,:)));
    a(3)=norm(cross(Jac(3,:),Jac(1,:)));
    h=det(Jac)/max(a);
    l=max(vecnorm(Jac'));
    quality(ei)=h/l;
    % S=Jac*invW;
    % quality(ei)=3*det(S)^(2/3)/sum(S.^2,'all');
end
end