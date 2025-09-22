%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Mesh=Creat_new_mesh()
Mesh.nodes=zeros(0,3);
Mesh.positions=zeros(0,3);
Mesh.nactivation=true(0,1);
Mesh.ninpart=zeros(0,1);
Mesh.elements=zeros(0,4);
Mesh.etype=zeros(0,1);
Mesh.einpart=zeros(0,1);
Mesh.ecenter=zeros(0,3);
Mesh.eactivation=true(0,1);
Mesh.surfaces=zeros(0,3);
Mesh.stype=zeros(0,1);
Mesh.sinpart=zeros(0,1);
Mesh.sinfacet=zeros(0,1);
Mesh.scenter=zeros(0,3);
Mesh.particles=zeros(0,3);
Mesh.pinpart=zeros(0,1);
end