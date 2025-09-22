%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Mesh=Import_mesh(Mesh,filename)
tic;data=importdata(filename,'',1e15);
%--------------------------------------------------------------------------
maxL=0;
for ii=1:1:size(data,1)
    maxL=max(maxL,size(data{ii},2));
end
maxL=maxL+5;
%------------------------------location---------------------------------------
star=[];partloc=[];instanceloc=[];
surfaceloc=[];stype=[];surface_name={};
for ii=1:1:size(data,1)
    if contains(data{ii},'*Part, name=')
        partloc(end+1,1)=ii;
    end
    if strcmp(data{ii},'*End Part')
        partloc(end,2)=ii;
    end
    if contains(data{ii},'*Instance, name=')
        instanceloc(end+1,1)=ii;
    end
    if strcmp(data{ii},'*End Instance')
        instanceloc(end,2)=ii;
    end
    if contains(data{ii},'*')
        star(end+1,1)=ii;
    end
    if strfind(data{ii},'*Elset, elset=_Cafte')
        surfaceloc(end+1)=ii;
        T=data{ii};T=split(T,",");T=split(T{2},"_");
        surface_name(end+1)=T(2);
        if strcmp(T{end},'SPOS')
            stype(end+1)=1;
        else
            stype(end+1)=str2num(T{end}(2:end));
        end
    end
    if size(data{ii},2)<maxL
        temp=data{ii};
        if ~strcmp(temp(end),',');temp(end+1)=',';end
        temp(size(temp,2)+1:maxL)='0';
        data{ii}=temp;
    end
end
%------------------------import part---------------------------------------
part.name=[];part.node=[];part.etype=[];part.element=[];
for ii=1:1:size(partloc,1)
    loc=star(star>=partloc(ii,1)&star<=partloc(ii,2));
    T=split(data{loc(1)},",");T=split(T{2},"name=");part(ii,1).name=T{end};
    if numel(loc)==2;continue;end
    node=str2num(cell2mat(data(loc(2)+1:loc(3)-1)));
    part(ii,1).node=node(:,2:end-1);
    for jj=3:1:numel(loc)-1
        T=split(data{loc(jj)},",");
        if strcmp(T{1},'*Element')
            T=split(T{2},"type=");
            element=str2num(cell2mat(data(loc(jj)+1:loc(jj+1)-1)));
            part(ii,1).element(end+1:end+size(element,1),1:size(element,2)-2)=element(:,2:end-1);
            part(ii,1).etype(end+1:end+size(element,1),1)=Get_element_type(T{end});
        end
    end
end
%------------------------import instance-----------------------------------
instance.name=[];instance.node=[];instance.etype=[];instance.element=[];
for ii=1:1:size(instanceloc,1)
    loc=star(star>=instanceloc(ii,1)&star<=instanceloc(ii,2));
    T=split(data{loc(1)},",");T1=split(T{2},"name=");instance(ii,1).name=T1{end};
    T2=split(T{3},"part=");ploc=find(strcmp(T2{end},{part.name}));
    instance(ii,1).node=part(ploc).node;instance(ii,1).etype=part(ploc).etype;instance(ii,1).element=part(ploc).element;
    if numel(loc)==2&&loc(2)-loc(1)==1;continue;end
    if numel(loc)>2
        node=str2num_modify(cell2mat(data(loc(2)+1:loc(3)-1)));
        instance(ii,1).node(end+1:end+size(node,1),1:size(node,2)-2)=node(:,2:end-1);
        for jj=3:1:numel(loc)-1
            T=split(data{loc(jj)},",");
            if strcmp(T{1},'*Element')
                T=split(T{2},"type=");
                element=str2num_modify(cell2mat(data(loc(jj)+1:loc(jj+1)-1)));
                instance(ii,1).element(end+1:end+size(element,1),1:size(element,2)-2)=element(:,2:end-1);
                instance(ii,1).etype(end+1:end+size(element,1),1)=Get_element_type(T{end});
            end
        end
    end
    for jj=loc(1)+1:1:loc(2)-1
        tran=str2num(data{jj});tran(end)=[];
        if numel(tran)<4
            instance(ii).node=instance(ii).node+tran(1:size(instance(ii).node,2));
        else
            w=tran([4,5,6])-tran([1,2,3]);theta=tran(7)/180*pi;
            w=w/norm(w);
            Espin=[w(1)^2*(1-cos(theta))+cos(theta),w(1)*w(2)*(1-cos(theta))+w(3)*sin(theta),w(1)*w(3)*(1-cos(theta))-w(2)*sin(theta);
                   w(1)*w(2)*(1-cos(theta))-w(3)*sin(theta),w(2)*w(2)*(1-cos(theta))+cos(theta),w(2)*w(3)*(1-cos(theta))+w(1)*sin(theta)
                   w(1)*w(3)*(1-cos(theta))+w(2)*sin(theta),w(2)*w(3)*(1-cos(theta))-w(1)*sin(theta),w(3)*w(3)*(1-cos(theta))+cos(theta)];
%             instance(ii).node=instance(ii).node-tran([1,2,3]);
            instance(ii).node=(instance(ii).node-tran([1,2,3]))*Espin+tran([1,2,3]);
        end
    end
end
instance(1).nodenumplus=0;instance(1).elementnumplus=0;
for ii=2:1:size(instance,1)
    instance(ii).nodenumplus=instance(ii-1).nodenumplus+size(instance(ii-1).node,1);
    instance(ii).elementnumplus=instance(ii-1).elementnumplus+size(instance(ii-1).element,1);
end
%-----------------------------get node-------------------------------------
nodes=[];ninpart=[];
for ii=1:size(instance,1)
    nodes=[nodes;instance(ii).node];
    ninpart(end+1:end+size(instance(ii).node,1),1)=ii;
end
%---------------------------get  element-----------------------------------
elements=[];etype=[];einpart=[];
for ii=1:size(instance,1)
    enode=instance(ii).element;enode(enode==0)=nan;
    elements(end+1:end+size(enode,1),1:size(enode,2))=enode+instance(ii).nodenumplus;
    etype(end+1:end+size(enode,1),1)=instance(ii).etype;
    einpart(end+1:end+size(enode,1),1)=ii;
end
elements(elements==0)=nan;
%-------------------------get surface--------------------------------------
[~,loc]=unique(surface_name);surface_name=surface_name(sort(loc));
Stype=[];surfaces=[];Sinpart=[];Sinfacet=[];
for ii=1:1:numel(surfaceloc)
    A=surfaceloc(ii);B=star(find(star==A)+1);
    if strfind(data{A},'generate')
        SE=str2num(cell2mat(data(A+1:B-1)));
        SE=SE(1):SE(3):SE(2);
    else
        temp=cell2mat(data(A+1:B-1));
        if size(temp,1)>1
            SE0=str2num(temp(1:end-1,:));
            SE0(:,end)=[];
        else
            SE0=[];
        end
        SE=[SE0(:).',str2num(temp(end,:))];
        SE(end)=[];
    end
    
    T=split(data{A},",");T1=split(T{2},"_");T2=split(T{4},"=");
    fID=find(strcmp(surface_name,T1{2}));
    pID=find(strcmp(T2{end},{instance.name}));
    SE=SE+instance(pID).elementnumplus;
    
    switch etype(SE(1))
        case 334
            facenode=[2,1,3;1,2,4;2,3,4;3,1,4];
        case 338
            facenode=[4,3,2,1;5,6,7,8;6,5,1,2;2,3,7,6;8,7,3,4;1,5,8,4]; 
        case 324
            facenode=[1,2,3,4];
        case {323,223}
            facenode=[1,2,3];
        case 3310
            facenode=[1,3,2,7,6,5;2,4,1,9,8,5;4,2,3,9,6,10;3,1,4,7,8,10];
        otherwise
            error('Unknown element type');
    end
    
    temp_surface=elements(SE,facenode(stype(ii),:));
    surfaces(end+1:end+numel(SE),1:size(facenode,2))=temp_surface;
    Stype(end+1:end+numel(SE),1)=sum(~isnan(temp_surface),2);
    Sinpart(end+1:end+numel(SE),1)=einpart(SE);
    Sinfacet(end+1:end+numel(SE),1)=fID;
end
surfaces(surfaces==0)=nan;
%=================================output===================================
Mesh.nodes=nodes;
Mesh.positions=nodes;
Mesh.ninpart=ninpart;
Mesh.nactivation=true(size(nodes,1),1);
%--------------------------------------------------------------------------
Mesh.elements=elements;
Mesh.etype=etype;
Mesh.einpart=einpart;
Mesh.eactivation=true(size(elements,1),1);
%--------------------------------------------------------------------------
Mesh.surfaces=surfaces;
Mesh.stype=Stype;
Mesh.sinpart=Sinpart;
Mesh.sinfacet=Sinfacet;
%==================================disp====================================
time=toc;fprintf("Model mesh is imported: %fs\n",time);
end
function code=Get_element_type(type)
switch type
    case {'C3D8','C3D8P'}
        code=338;
    case 'C3D8R'
        code=338;
    case {'C3D4','C3D4P','C3D4R'}
        code=334;
    case 'B31'
        code=312;
    case 'S3'
        code=323;
    case 'S4R'
        code=324;
    case 'S4R5'
        code=324;
    case 'S4'
        code=324;
    case 'M3D4'
        code=324;
    case 'CPS3'
        code=223;
    case 'CPE3'
        code=223;
    case 'CPS4'
        code=224;
    case 'C3D10'
        code=3310;
    otherwise
        error(['Unknown element type: ' num2str(type)]);
end
end

function numat=str2num_modify(text)
if isempty(text);numat=[];return;end
num=size(text,1);
n=ceil(num/500000);
num2=size(str2num(text(1,:)),2);
numat=zeros(num,num2);
for i=1:1:n-1
numat((i-1)*500000+1:i*500000,:)=str2num(text((i-1)*500000+1:i*500000,:));
end
numat((n-1)*500000+1:end,:)=str2num(text((n-1)*500000+1:end,:));
end