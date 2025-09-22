%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Import_FE_mesh(Model,nodeCell,elementCell,etypeCell,partIdCell)
coder.varsize('eactivationCell')
eactivationCell=coder.nullcopy(cell(numel(elementCell),1));
for pii=1:numel(eactivationCell)
    eactivationCell{pii}=true(size(elementCell{pii},1),1);
end
%==========================================================================
Mesh=Model.Mesh;
Mesh.elements(Mesh.elements==0)=nan;
loc=ismember(Mesh.einpart,Montage_cell(partIdCell,1));loc=~loc;
nodes0=Mesh.nodes;
elements0=Mesh.elements(loc,:);
Einpart0=Mesh.einpart(loc,:);
Etype0=Mesh.etype(loc,:);
Eactivation0=Mesh.eactivation(loc,:);
partid0=unique(Einpart0);
flag=numel(elementCell);
for pii=1:1:numel(partid0)
    loc=Einpart0==partid0(pii);
    temp_element=elements0(loc,:);
    [nid,~,temp_element]=unique(temp_element);
    nid(isnan(nid))=[];
    temp_element(temp_element>numel(nid))=nan;
    elementCell{flag+1}=reshape(temp_element,[],size(elements0,2));
    nodeCell{flag+1}=nodes0(nid,:);
    etypeCell{flag+1}=Etype0(loc);
    eactivationCell{flag+1}=Eactivation0(loc);
    partIdCell{flag+1}=partid0(pii);
end
[~,loc]=sort(Montage_cell(partIdCell,1));
nodeCell=nodeCell(loc);
elementCell=elementCell(loc);
etypeCell=etypeCell(loc);
eactivationCell=eactivationCell(loc);
partIdCell=partIdCell(loc);
%==========================================================================
coder.varsize('ninpartCell')
ninpartCell=coder.nullcopy(cell(numel(nodeCell),1));
coder.varsize('einpartCell')
einpartCell=coder.nullcopy(cell(numel(elementCell),1));
flag=0;
for pii=1:numel(elementCell)
    pid=partIdCell{pii};
    nnum=size(nodeCell{pii},1);
    enum=size(elementCell{pii},1);
    elementCell{pii}=elementCell{pii}+flag;
    flag=flag+nnum;
    ninpartCell{pii}=repmat(pid,nnum,1);
    einpartCell{pii}=repmat(pid,enum,1);
end
%==========================================================================
Mesh=Creat_new_mesh();
Mesh.nodes=Montage_cell(nodeCell,1);
Mesh.positions=Mesh.nodes;
Mesh.ninpart=Montage_cell(ninpartCell,1);
Mesh.nactivation=false(size(Mesh.nodes,1),1);
Mesh.elements=Montage_cell(elementCell,1);
Mesh.elements(Mesh.elements==0)=nan;
Mesh.etype=Montage_cell(etypeCell,1);
Mesh.einpart=Montage_cell(einpartCell,1);
Mesh.ecenter=Get_geometry_center(Mesh.elements,Mesh.nodes);
Mesh.eactivation=Montage_cell(eactivationCell,1);
nid=Mesh.elements(Mesh.eactivation,:);nid(isnan(nid))=[];
Mesh.nactivation(nid(:))=true;
Mesh=Get_surface(Mesh);
Mesh.particles=Model.Mesh.particles;
Mesh.pinpart=Model.Mesh.pinpart;
%==========================================================================
Model.Mesh=Mesh;
end