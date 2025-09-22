%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [mortar_matrix,R_matrix,area_vector,slave_dof]=node_to_surface(slave_surface,master_surface,nodes,plane_orientation)
maxdof=size(nodes,1)*3;
if isempty(slave_surface)||isempty(master_surface)
    slave_dof=zeros(maxdof-maxdof,1);
    area_vector=zeros(maxdof-maxdof,1);
    R_matrix=sparse(maxdof-maxdof,maxdof-maxdof);
    mortar_matrix=sparse(maxdof-maxdof,maxdof);
    return;
end
%---------------------------------------------------------------------------------
master_centers=(nodes(master_surface(:,1),:)+nodes(master_surface(:,2),:)+nodes(master_surface(:,3),:))/3;
slave_node=unique(slave_surface);
slave_area=zeros(size(nodes,1),1);
slave_area_vectors=zeros(size(nodes,1),3);
for ssi=1:1:size(slave_surface,1)
    nid=slave_surface(ssi,1:3);
    slave_fcoor=nodes(nid,1:3);
    normal_vector=cross(slave_fcoor(2,:)-slave_fcoor(1,:),slave_fcoor(3,:)-slave_fcoor(1,:));
    A=norm(normal_vector);
    slave_area(nid)=slave_area(nid)+A/6;
    slave_area_vectors(nid,:)=slave_area_vectors(nid,:)+repmat(normal_vector/6,numel(nid),1);
end
%-----------------------------mortar matrix initialization-----------------
RV=zeros(size(slave_node,1),3);
Area=zeros(size(slave_node,1),1);
DI=zeros(size(slave_node,1),1);
MV=zeros(size(slave_node,1),3);
MJ=zeros(size(slave_node,1),3);
for ssi=1:1:numel(slave_node)
    nid=slave_node(ssi);
    slave_fcoor=nodes(nid,1:3);
    A=slave_area(nid);
    %----------------------------------------------------------------------
    distance=master_centers-repmat(slave_fcoor,size(master_centers,1),1);
    distance=sum(abs(distance),2);[~,sort_id]=sort(distance);
    %----------------------------------------------------------------------
    %figure(100);hold on;fill3(slave_fcoor(:,1),slave_fcoor(:,2),slave_fcoor(:,3),'r')
    for ii=1:1:size(distance,1)
        id=sort_id(ii);
        if distance(id)>sqrt(A)
            break
        end
        nj=master_surface(id,1:3);
        master_fcoor=nodes(nj,1:3);
        master_normal=cross(master_fcoor(2,:)-master_fcoor(1,:),master_fcoor(3,:)-master_fcoor(1,:));
        master_normal=master_normal./norm(master_normal);
        if (slave_fcoor-mean(master_fcoor,1))*master_normal.'>0.01*sqrt(A)
        continue
        end
        R_system=vector_to_system(master_normal).';
        master_fcoorp=master_fcoor*R_system(:,1:2);
        slave_fcoorp=slave_fcoor*R_system(:,1:2);
        isin=inpoly(slave_fcoorp,master_fcoorp);
        %isin=inpolygon(slave_fcoorp(:,1),slave_fcoorp(:,2),master_fcoorp(:,1),master_fcoorp(:,2));
        if isin
            %figure(100);hold on;fill3(master_fcoor(:,1),master_fcoor(:,2),master_fcoor(:,3),'b')
            ip_local=Triangular_inverse_interpolation(master_fcoorp,slave_fcoorp);
            MN=Triangular_interpolation(ip_local);
            DI(ssi)=nid;
            Area(ssi)=A;
            RV(ssi,1:3)=-master_normal;
            MV(ssi,1:3)=MN;
            MJ(ssi,1:3)=nj;
            break
        end
    end
end
%--------------------------------------------------------------------------
loc=DI~=0;
DI=DI(loc);
% RV=RV(loc,:)*plane_orientation;
RV=slave_area_vectors(DI,:)./slave_area(DI);
Area=Area(loc,:);
MV=MV(loc,:);
MJ=MJ(loc,:);
%--------------------------------------------------------------------------
RN=RV./repmat(sqrt(RV(:,1).^2+RV(:,2).^2+RV(:,3).^2),1,3);
RVS=zeros(size(RN,1)*3,3);
for ii=1:1:size(RN,1)
    RVS(ii*3-2:ii*3,:)=vector_to_system(RN(ii,:));
end
%--------------------------------------------------------------------------
DI=[DI*3-2,DI*3-1,DI*3];
%--------------------------------------------------------------------------
area_vector=sparse(DI,1,[Area,Area,Area],maxdof,1);
M_matrix=sparse(DI(:,[1,1,1,2,2,2,3,3,3]),[MJ*3-2,MJ*3-1,MJ*3],[MV,MV,MV],maxdof,maxdof);
R_matrix=sparse(DI(:,[1,1,1,2,2,2,3,3,3]),DI(:,[1,2,3,1,2,3,1,2,3]),[RVS(1:3:end,:),RVS(2:3:end,:),RVS(3:3:end,:)],maxdof,maxdof);
%--------------------------------------------------------------------------
slave_dof=unique(DI(:));slave_dof(slave_dof==0)=[];
R_matrix=R_matrix(slave_dof,slave_dof)*plane_orientation;
area_vector=full(area_vector(slave_dof,:));
D_matrix=sparse(1:numel(slave_dof),slave_dof,1,numel(slave_dof),maxdof);
M_matrix=M_matrix(slave_dof,:);
mortar_matrix=(D_matrix-M_matrix);
end

function flag=inpoly(coor1,coor2)
s=zeros(size(coor1,1),size(coor2,1));
coor2=[coor2;coor2(1,:)];
for i=1:1:size(coor1,1)
    polygon=coor2-repmat(coor1(i,:),size(coor2,1),1);
    for j=1:1:size(polygon,1)-1
        s(i,j)=polygon(j,1)*polygon(j+1,2)-polygon(j+1,1)*polygon(j,2);
    end
end
A=abs(sum(s(1,:)));
flag=Truncated_errors(max(s,[],2).*min(s,[],2)/A,8)>=0;
end