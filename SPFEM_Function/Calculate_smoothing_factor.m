%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Fai,node_volume]=Calculate_smoothing_factor(Mesh,CEdstrain,CEid)
nodes=Mesh.nodes;
elements=Mesh.elements;
Nactivation=Mesh.nactivation;
node_num=size(nodes,1);
element_num=size(elements,1);
Evolumes=zeros(element_num,1);
for ei=1:1:element_num
    nid=elements(ei,1:4);
    Evolumes(ei)=det(nodes(nid(2:end),:)-nodes(nid(1),:));
end
Evolumes=Evolumes/6;
%--------------------------------------------------------------------------
Ecenter=Get_geometry_center(elements,nodes);
Distance=zeros(element_num,size(elements,2));
for i=1:1:size(elements,2)
    Distance(:,i)=1./vecnorm(Ecenter-nodes(elements(:,i),:),2,2);
end
Distance=sum(Distance,2).\Distance*4;
%--------------------------------------------------------------------------
eis=(1:1:size(elements,1))';
fai=sparse(node_num,element_num);
for i=1:1:size(elements,2)
fai=fai+sparse(elements(:,i),eis,Evolumes.*Distance(:,i),node_num,element_num);
end
node_volume=sum(fai,2);
node_volume(~Nactivation)=1e-12;
fai=diag(node_volume)\fai;
%==========================================================================
Nconnection=Get_node_connection(elements);
normal_area=cell(node_num,1);
normal_dof=cell(node_num,1);
num=0;
for ni=1:1:node_num
    if Nactivation(ni)
        eid=Nconnection(ni,1:Nconnection(ni,end))';
    else
         eid=zeros(1,0);
    end
    [nsurface,inelememt]=Get_surface_temp(elements(eid,:));
    nseid=eid(inelememt);
    L1=nodes(nsurface(:,2),:)-nodes(nsurface(:,1),:);
    L2=nodes(nsurface(:,3),:)-nodes(nsurface(:,1),:);
    outwardsXarea=cross(L1,L2,2)/2;
    outwardsXarea(abs(outwardsXarea)<1e-12)=0;
    normal_area{ni}=outwardsXarea;
    normal_dof{ni}=nseid;
    num=num+numel(nseid);
end
%==========================================================================
dfai_area=zeros(num,3);
dof4dfai=zeros(num,2);
flag=1;
for ni=1:1:node_num
    flagend=flag+numel(normal_dof{ni})-1;
    dfai_area(flag:flagend,:)=normal_area{ni};
    dof4dfai(flag:flagend,:)=[repmat(ni,numel(normal_dof{ni}),1),normal_dof{ni}];
    flag=flagend+1;
end
dfaidx_area=sparse(dof4dfai(:,1),dof4dfai(:,2),dfai_area(:,1),node_num,element_num);
dfaidy_area=sparse(dof4dfai(:,1),dof4dfai(:,2),dfai_area(:,2),node_num,element_num);
dfaidz_area=sparse(dof4dfai(:,1),dof4dfai(:,2),dfai_area(:,3),node_num,element_num);
%==========================================================================
r=(0.75/pi).^(1/3)*(node_volume/4).^(1/3);
dfaidx=diag(node_volume.\r)*dfaidx_area;
dfaidy=diag(node_volume.\r)*dfaidy_area;
dfaidz=diag(node_volume.\r)*dfaidz_area;
node_volume=node_volume/4;
node_volume(~Nactivation)=0;
%==========================================================================
Estrain_step=zeros(element_num,6);
for ci=1:1:numel(CEdstrain)
Estrain_step(CEid{ci},:)=permute(CEdstrain{ci},[4,1,3,2]);
end
strain=fai'*(node_volume.*(fai*Estrain_step));
dstrain=dfaidx'*(node_volume/3.*(dfaidx*Estrain_step))...
       +dfaidy'*(node_volume/3.*(dfaidy*Estrain_step))...
       +dfaidz'*(node_volume/3.*(dfaidz*Estrain_step));
% strain=Evolumes.\strain;
% dstrain=Evolumes.\dstrain;
D=Form_D_elastic(1,0.3);
A=sum(Evolumes.*sum((D*Estrain_step')'.*Estrain_step,2),1);
B=sum(strain.*(D*Estrain_step')','all');
C=sum(dstrain.*(D*Estrain_step')','all');
alpha=sqrt((A-B)/C/2);
alpha=sqrt(2);
if isnan(alpha);alpha=1;end
%================================Output====================================
Fai{1}=fai+alpha*dfaidx;
Fai{2}=fai-alpha*dfaidx;
Fai{3}=fai+alpha*dfaidy;
Fai{4}=fai-alpha*dfaidy;
Fai{5}=fai+alpha*dfaidz;
Fai{6}=fai-alpha*dfaidz;
end