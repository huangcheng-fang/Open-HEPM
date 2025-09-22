%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [mortar_matrix,R_matrix,area_vector,slave_dof]=Surface_to_volume(slave_surface,master_element,Stype,Etype,nodes,isdual)
maxdof=size(nodes,1)*3;
if isempty(slave_surface)||isempty(master_element)
    slave_dof=zeros(maxdof-maxdof,1);
    area_vector=zeros(maxdof-maxdof,1);
    R_matrix=sparse(maxdof-maxdof,maxdof-maxdof);
    mortar_matrix=sparse(maxdof-maxdof,maxdof);
    return;
end
%==============Hammer integral and Surface center==========================
[ip,weight]=Triangular_integral(2);
slave_centers=Get_geometry_center(slave_surface,nodes);
master_centers=Get_geometry_center(master_element,nodes);
%=======================mortar matrix initialization=======================
Value=zeros(size(slave_surface,1)*2000,1);
IndexI=zeros(size(slave_surface,1)*2000,1);
IndexJ=zeros(size(slave_surface,1)*2000,1);
area_vector=zeros(maxdof/3,1);
Rn_vector=zeros(maxdof/3,3);
%----------------------------------------------------------------------
flaglocd2=0;%flagwarning=[1,1];
for ssi=1:1:size(slave_surface,1)
    nid=slave_surface(ssi,:);nid(isnan(nid))=[];
    slave_type=Stype(ssi);
    slave_fcoor=nodes(nid,1:3);
    %----------------------------------------------------------------------
    normal_vector=Get_normal_vector(slave_fcoor);
    R_system=vector_to_system(normal_vector).';
    slave_fcoorp=(slave_fcoor-slave_centers(ssi,1:3))*R_system;
    A=polyarea(slave_fcoorp(:,1),slave_fcoorp(:,2));
    %----------------------------------------------------------------------
    dcoor=master_centers-repmat(slave_centers(ssi,1:3),size(master_centers,1),1);
    distance=sum(abs(dcoor),2);
    active_range=find(distance<sqrt(A)*5);
    distance=distance(active_range,:);
    %----------------------------------------------------------------------
    dsize=numel(nid);
    flagA=A;DK=zeros(dsize);MKS=zeros(dsize,0);mnode=zeros(1,0);
    %figure(1);plot(polyshape(slave_fcoorp(:,1:2)));hold on;
    % figure(100); fill3(slave_fcoor(:,1),slave_fcoor(:,2),slave_fcoor(:,3),'r');hold on
    for ii=1:1:size(distance,1)
        if isempty(distance);break;end
        [~,loc]=min(distance);id=active_range(loc);distance(loc)=[];active_range(loc)=[];
        nj=master_element(id,:);nj(isnan(nj))=[];
        master_type=Etype(id,1);
        master_fcoor=nodes(nj,1:3);
        master_fcoorp=(master_fcoor-repmat(slave_centers(ssi,1:3),size(master_fcoor,1),1))*R_system;
        Vertices=Cross_surface_to_volume(slave_fcoorp,master_fcoorp,master_type,false);
        % figure(100);plot3(Vertices(:,1),Vertices(:,2),Vertices(:,3),'*');patch('Faces',[2,1,3;4,1,2;3,1,4;2,3,4],'Vertices',master_fcoor,'FaceColor','w');hold on
        if size(Vertices,1)>2
            [ConnectivityList,Points]= Overlapping_triangulation(Vertices);
            MK=zeros(dsize,numel(nj));
            mnode=[mnode,nj];
            for ti=1:1:size(ConnectivityList,1)
                tnode_coor=Points(ConnectivityList(ti,:),1:3);
                detJ=(tnode_coor(2,1)-tnode_coor(1,1))*(tnode_coor(3,2)-tnode_coor(1,2))-(tnode_coor(2,2)-tnode_coor(1,2))*(tnode_coor(3,1)-tnode_coor(1,1));
                if abs(detJ/A)<1e-6;continue;end
                if detJ<0;error('Incorrect triangular');end
                flagA=flagA-detJ/2;
                %----------------------------------------------------------------------------------------------
                PN=[1-ip(:,1)-ip(:,2),ip(:,1),ip(:,2)];
                pglobal_coor=PN*tnode_coor;
                slave_ip_local=Inverse_interpolation(slave_type,slave_fcoorp,pglobal_coor);
                master_ip_local=Inverse_interpolation(master_type,master_fcoorp,pglobal_coor);
                %----------------------------------------------------------------------------------------------
                SN=Get_shape_function(slave_type,slave_ip_local);
                MN=Get_shape_function(master_type,master_ip_local);
                %----------------------------------------------------------------------------------------------
                smweights=detJ*weight;
                DK=DK+SN'*diag(smweights)*SN;
                MK=MK+SN'*diag(smweights)*MN;
            end
            MKS=[MKS,MK];
        end
        %---------------------------------------------------------------------------------------------
        if abs(flagA/A)<1e-6;break;end
        %if ii>100&&flagwarning(1);disp('Waring: Exceeded the maximum number of contact pair search');flagwarning(1)=0;end
        %--------------------------------------------------------------------------------------------------
    end
    %------------------------------------------------------------------------------------------------------------
    %if abs(flagA/A)>1e-6&&flagwarning(2);disp('Warning: Slave surface is not completely overlapped');flagwarning(2)=0;end
    if abs(flagA-A)/A<1e-6;continue;end
    %-----------------------------------------------------------------------------------------------------------
    [mnode,~,loc2]=unique(mnode);
    CMK=zeros(size(MKS,1),numel(mnode));
    for ii=1:1:size(MKS,2)
        CMK(:,loc2(ii))=CMK(:,loc2(ii))+MKS(:,ii);
    end
    %----------------------------------------------------------------------
    if isdual
        DF=sum(DK,2);
        dual_Lagrange_coefficient=DK\diag(DF);
        DK=diag(DF);
        CMK=dual_Lagrange_coefficient.'*CMK;
    end
    %----------------------------------------------------------------------
    DM=[DK,-CMK];DMid=[nid,mnode];
    idxI=repmat(nid(:),1,numel(DMid));idxJ=repmat(DMid,numel(nid),1);
    flaglocd1=flaglocd2+1;flaglocd2=flaglocd2+numel(DM);
    Value(flaglocd1:flaglocd2)=DM(:);
    IndexI(flaglocd1:flaglocd2)=idxI(:);
    IndexJ(flaglocd1:flaglocd2)=idxJ(:);
    %----------------------------------------------------------------------
    area=sum(DK,2);
    area_vector(nid)=area_vector(nid)+area;
    %----------------------------------------------------------------------
    Rn_vector(nid,1:3)=Rn_vector(nid,1:3)+area*normal_vector;
end
%==========================Check for memory overflow=======================
if numel(Value)>size(slave_surface,1)*2000
    disp('Waring: Resulting mortar matrix is larger than predefined mortar matrix');
end
%===========================Get mortar matrix===============================
Value(flaglocd2+1:end)=[];IndexI(flaglocd2+1:end)=[];IndexJ(flaglocd2+1:end)=[];
%-----------------------------Extend to 3D---------------------------------
mortar_matrix=sparse(IndexI,IndexJ,Value);
[IndexI,IndexJ,Value]=find(mortar_matrix);
IndexI=IndexI*3;IndexJ=IndexJ*3;
IndexI=[IndexI-2,IndexI-1,IndexI];
IndexJ=[IndexJ-2,IndexJ-1,IndexJ];
Value=repmat(Value,1,3);
%----------------------------Rearrangement---------------------------------
slave_dof=false(maxdof,1);slave_dof(IndexI)=true;slave_dof=find(slave_dof);
index=zeros(maxdof/3,1);index(slave_dof)=1:1:numel(slave_dof);IndexI=index(IndexI);
%----------------------Reassembly and normalization------------------------
area_vector=repmat(area_vector,1,3).';area_vector=area_vector(slave_dof);
mortar_matrix=sparse(IndexI,IndexJ,Value,numel(slave_dof),maxdof);
mortar_matrix=area_vector.\(mortar_matrix);
%============================Get R matrix==================================
Rn_vector=Rn_vector(slave_dof(3:3:end)/3,:);
Rn_vector=sqrt(sum(Rn_vector.^2,2)).\Rn_vector;%normalization
%-----------------------------Extend to 3D---------------------------------
RV=zeros(3,3,size(Rn_vector,1));
for ii=1:1:size(Rn_vector,1)
    RV(:,:,ii)=vector_to_system(Rn_vector(ii,:));
end
%------------------------------Assembly------------------------------------
RI=repmat(reshape(1:1:numel(slave_dof),3,1,[]),1,3,1);
RJ=permute(RI,[2,1,3]);
R_matrix=sparse(RI(:),RJ(:),RV(:));
end