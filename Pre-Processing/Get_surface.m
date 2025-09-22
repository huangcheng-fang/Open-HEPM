%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Mesh=Get_surface(Mesh)
Eactivation=Mesh.eactivation;
elements=Mesh.elements(Eactivation,:);
Etype=Mesh.etype(Eactivation,:);
Einpart=Mesh.einpart(Eactivation,:);
facetnid=save_infacet(Mesh.surfaces,Mesh.sinfacet);
%------------------------------------------------------------------------
fsize=0;flag=0;
faces=zeros(size(elements,1)*8,3);finfo=zeros(size(elements,1)*8,2);
loc=find(Etype==223);
if ~isempty(loc)
fsize=2;
face=elements(loc,[1 2]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[2+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[2,3]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[2+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[3,1]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[2+loc*0,loc];
flag=flag+size(face,1);
end

loc=find(Etype==334);
if ~isempty(loc)
fsize=3;
face=elements(loc,[2,1,3]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[3+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[4,1,2]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[3+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[3,1,4]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[3+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[2,3,4]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[3+loc*0,loc];
flag=flag+size(face,1);
end

loc=find(Etype==3310);
if ~isempty(loc)
face=elements(loc,[1,3,2,7,6,5]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[6+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[2,4,1,9,8,5]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[6+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[4,2,3,9,6,10]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[6+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[3,1,4,7,8,10]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[6+loc*0,loc];
flag=flag+size(face,1);
end


loc=find(Etype==338);
if ~isempty(loc)
fsize=4;
face=elements(loc,[4,3,2,1]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[5,6,7,8]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[6,5,1,2]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[2,3,7,6]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[8,7,3,4]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
face=elements(loc,[1,5,8,4]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
end

loc=find(Etype==323);
if ~isempty(loc)
face=elements(loc,[1,2,3]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[3+loc*0,loc];
flag=flag+size(face,1);
end

loc=find(Etype==324);
if ~isempty(loc)
face=elements(loc,[1,2,3,4]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[4+loc*0,loc];
flag=flag+size(face,1);
end

loc=find(Etype==312);
if ~isempty(loc)
face=elements(loc,[1,2]);
faces(flag+1:flag+size(face,1),1:size(face,2))=face;
finfo(flag+1:flag+size(face,1),1:2)=[2+loc*0,loc];
flag=flag+size(face,1);
end
%--------------------------------------------------------------------------
faces(flag+1:end,:)=[];finfo(flag+1:end,:)=[];
%----------------------------------------------------------------------------
[~,b,c]=unique(sort(faces,2),'rows');
fnum=accumarray(c,ones(numel(c),1),[numel(b),1]);
faces=faces(b,:);
finfo=finfo(b,:);
surfaces=faces(fnum==1,:);
sinfov=finfo(fnum==1,:);
surfaces(surfaces==0)=nan;
%--------------------------------------------------------------------------
Mesh.surfaces=surfaces;
Mesh.stype=sinfov(:,1);
Mesh.sinpart=Einpart(sinfov(:,2),1);
Mesh.scenter=Get_geometry_center(Mesh.surfaces,Mesh.nodes);
%--------------------------------------------------------------------------
Mesh.sinfacet=resume_infacet(surfaces,facetnid);
end

function facetnid=save_infacet(surfaces,sinfacet)
fid=unique(sinfacet);
facetnid=cell(numel(fid),1);
for id=1:1:numel(fid)
nid=unique(surfaces(sinfacet==fid(id),:));
facetnid{id}=nid;
end
end

function sinfacet=resume_infacet(surfaces,facetnid)
sinfacet=zeros(size(surfaces,1),1);
for id=1:1:numel(facetnid)
loc0=find(sinfacet==0);
loc=all(ismember(surfaces(loc0,:),facetnid{id}),2);
sinfacet(loc0(loc))=id;
end
end