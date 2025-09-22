%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_result(Result,Mesh,n,type,scale,rangePart,LineStyle,Permutation)
if nargin<7
    LineStyle='-';
end
if nargin<8
    Permutation=[1,2,3];
end
dim=size(Mesh.nodes,2);
%--------------------------------------------------------------------------
Nfield=Result.Nfield;
if isfield(Result,'Pfield')
    NEfield=Result.Pfield;
    NEfield.Stress=mean(NEfield.Stress,3);
    NEfield.Strain=mean(NEfield.Strain,3);
    NEfield.Pstrain=mean(NEfield.Pstrain,3);
    NEfield.Parameter=mean(NEfield.Parameter,3);
end
if isfield(Result,'Efield')
    NEfield=Mapping_field_E2N(Result.Efield,Mesh);
end
%--------------------------------------------------------------------------
switch type
    case 'Sx'
        value=NEfield.Stress(1,:);
    case 'Sy'
        value=NEfield.Stress(2,:);
    case 'Sz'
        value=NEfield.Stress(3,:);
    case 'Sxy'
        value=NEfield.Stress(4,:);
    case 'Syz'
        value=NEfield.Stress(5,:);
    case 'Sxz'
        value=NEfield.Stress(6,:);
    case 'Sm'
        value=mean(NEfield.Stress(1:3,:),1);
    case 'S_mises'
        stress=NEfield.Stress;
        J2=((stress(1,:)-stress(2,:)).^2+(stress(2,:)-stress(3,:)).^2+(stress(3,:)-stress(1,:)).^2)/6+(stress(4,:).^2+stress(5,:).^2+stress(6,:).^2);
        value=sqrt(3*J2);
    case 'Min_prin_stress'
        e=NEfield.Stress;
        e(4:6,:,:,:)=e(4:6,:,:,:)/2;
        em=(e(1,:,:,:)+e(2,:,:,:)+e(3,:,:,:))/3;
        e(1:3,:,:,:)=e(1:3,:,:,:)-em;
        J2=0.5*(e(1,:,:,:).^2+e(2,:,:,:).^2+e(3,:,:,:).^2)+e(4,:,:,:).^2+e(5,:,:,:).^2+e(6,:,:,:).^2;
        J3=e(1,:,:,:).*e(2,:,:,:).*e(3,:,:,:)+2.*e(4,:,:,:).*e(5,:,:,:).*e(6,:,:,:)-e(1,:,:,:).*e(5,:,:,:).^2-e(2,:,:,:).*e(6,:,:,:).^2-e(3,:,:,:).*e(4,:,:,:).^2;
        Lode=asin(-1.5*sqrt(3).*J3./(J2).^1.5)/3;
        %--------------------------------------------------------------------------
        value=2*sqrt(J2/3).*sin(Lode-2*pi/3)+em;
    case 'Ex'
        value=NEfield.Strain(1,:);
    case 'Ey'
        value=NEfield.Strain(2,:);
    case 'Ez'
        value=NEfield.Strain(3,:);
    case 'Exy'
        value=NEfield.Strain(4,:);
    case 'Eyz'
        value=NEfield.Strain(5,:);
    case 'Exz'
        value=NEfield.Strain(6,:);
    case 'PEx'
        value=NEfield.Pstrain(1,:);
    case 'PEy'
        value=NEfield.Pstrain(2,:);
    case 'PEz'
        value=NEfield.Pstrain(3,:);
    case 'PExy'
        value=NEfield.Pstrain(4,:);
    case 'PEyz'
        value=NEfield.Pstrain(5,:);
    case 'PExz'
        value=NEfield.Pstrain(6,:);
    case 'PEEQ'
        pstrain=NEfield.Pstrain;
        pe=(pstrain(1,:)-pstrain(2,:)).^2+(pstrain(2,:)-pstrain(3,:)).^2+(pstrain(3,:)-pstrain(1,:)).^2+(pstrain(4,:).^2+pstrain(5,:).^2+pstrain(6,:).^2)*1.5;
        value=sqrt(2/9*pe);
    case 'Max_prin_strain'
        e=NEfield.Strain;
        e(4:6,:,:,:)=e(4:6,:,:,:)/2;
        em=(e(1,:,:,:)+e(2,:,:,:)+e(3,:,:,:))/3;
        e(1:3,:,:,:)=e(1:3,:,:,:)-em;
        J2=0.5*(e(1,:,:,:).^2+e(2,:,:,:).^2+e(3,:,:,:).^2)+e(4,:,:,:).^2+e(5,:,:,:).^2+e(6,:,:,:).^2;
        J3=e(1,:,:,:).*e(2,:,:,:).*e(3,:,:,:)+2.*e(4,:,:,:).*e(5,:,:,:).*e(6,:,:,:)-e(1,:,:,:).*e(5,:,:,:).^2-e(2,:,:,:).*e(6,:,:,:).^2-e(3,:,:,:).*e(4,:,:,:).^2;
        Lode=asin(-1.5*sqrt(3).*J3./(J2).^1.5)/3;
        %--------------------------------------------------------------------------
        value=2*sqrt(J2/3).*sin(Lode+2*pi/3)+em;
    case 'Ux'
        value=Nfield.U(1:dim:end);
    case 'Uy'
        value=Nfield.U(2:dim:end);
    case 'Uz'
        value=Nfield.U(3:dim:end);
    case 'U'
        value=sum(abs(reshape(Nfield.U,3,[])),1);
    case 'Pore_pressure'
        value=P;
    otherwise
        error(['Wrong Keyword: ',type]);
end
%==============================main function===============================
nodes=Mesh.positions+(scale-1)*reshape(Nfield.U,[],size(Mesh.nodes,1)).';
nodes=nodes(:,Permutation);
elements=Mesh.elements;
ninpart=Mesh.ninpart;
surfaces=Mesh.surfaces;
sinpart=Mesh.sinpart;
stype=Mesh.stype;
surfaces(stype==6,4:end)=nan;
einpart=Mesh.einpart;
etype=Mesh.etype;
flag=isempty(rangePart);
figure(n)
clf(n)
if nargin<6||flag
    rangePart=unique(sinpart);
end
for ii=1:1:numel(rangePart)
    patch('Faces',surfaces(sinpart==rangePart(ii),:),'Vertices',nodes,'FaceVertexCData',value(:),'FaceColor','interp','LineWidth',0.001,'LineStyle',LineStyle);
end
if nargin<6||flag
    rangePart=unique(einpart);
end
for ii=1:1:numel(rangePart)
    patch('Faces',elements((etype==223)&einpart==rangePart(ii),1:end-1),'Vertices',nodes,'FaceVertexCData',value(:),'FaceColor','interp','LineWidth',0.001,'LineStyle',LineStyle);
    patch('Faces',elements((etype==312)&einpart==rangePart(ii),1:end-1),'Vertices',nodes,'FaceVertexCData',value(:),'EdgeColor','interp','FaceColor','interp','LineWidth',3);
end
%--------------------------------------------------------------------------
if size(nodes,2)==3
    view(-18,15)
end
axis equal
colorbar
jet_reduce=jet(13);
colormap(jet_reduce(1:12,:))
loc=ismember(ninpart,rangePart);
bar_min_max=[min(value(loc)),max(value(loc))+1e-16];
clim(bar_min_max)
colorbar('ytick',linspace(bar_min_max(1),bar_min_max(2),7));
drawnow
end