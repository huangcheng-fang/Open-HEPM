function CBar=Plot_result_on_particle(Result,Mesh,figureID,vari,range_part,permutation)
%==========================Check input=====================================
bar_min_max=[];cnum=11;
%==============================main function===============================
if isfield(Mesh,'Mesh')
    Mesh=Mesh.Mesh;
end
%--------------------------------------------------------------------------
particles=Mesh.particles;
particle_partID=Mesh.pinpart;
particles=particles(:,permutation);
Particle_Field=Result.Pfield;
vari=Get_var(vari,Particle_Field);
if ~isempty(bar_min_max)
vari(vari<bar_min_max(1))=bar_min_max(1);
vari(vari>bar_min_max(2))=bar_min_max(2);
end
if isempty(range_part)
    range_part=unique(particle_partID);
end
%--------------------------------------------------------------------------
figure(figureID)
clf(figureID)
hold on;
cmap=colormap(jet(cnum));
value=round(mapminmax(vari,1,cnum));
for ii=1:1:numel(range_part)
    pid=particle_partID==range_part(ii);
    scatter3(particles(pid,1),particles(pid,2),particles(pid,3),10,cmap(value(pid),:),'filled')
end
%--------------------------------------------------------------------------
if size(particles,2)==3
    view(-18,15)
end
if isempty(bar_min_max)
    loc=ismember(particle_partID,range_part);
    bar_min_max=[min(vari(loc)),max(vari(loc))+1e-16];
end
axis off
axis equal
jet_reduce=jet(cnum+1);
colormap(jet_reduce(1:cnum,:))
clim(bar_min_max)
CBar=colorbar('ytick',linspace(bar_min_max(1),bar_min_max(2),ceil(cnum/2)+1));
drawnow
end


function value=Get_var(var,Particle_Field)
switch var
    case 'Sx'
        value=Particle_Field.Stress(1,:);
    case 'Sy'
        value=Particle_Field.Stress(2,:);
    case 'Sz'
        value=Particle_Field.Stress(3,:);
    case 'Sxy'
        value=Particle_Field.Stress(4,:);
    case 'Syz'
        value=Particle_Field.Stress(5,:);
    case 'Sxz'
        value=Particle_Field.Stress(6,:);
    case 'Sm'
        value=mean(Particle_Field.Stress(1:3,:),1);
    case 'S_mises'
        stress=Particle_Field.Stress;
        J2=((stress(1,:)-stress(2,:)).^2+(stress(2,:)-stress(3,:)).^2+(stress(3,:)-stress(1,:)).^2)/6+(stress(4,:).^2+stress(5,:).^2+stress(6,:).^2);
        value=sqrt(3*J2);
    case 'Min_prin_stress'
        e=Particle_Field.Stress;
        e(4:6,:,:,:)=e(4:6,:,:,:)/2;
        em=(e(1,:,:,:)+e(2,:,:,:)+e(3,:,:,:))/3;
        e(1:3,:,:,:)=e(1:3,:,:,:)-em;
        J2=0.5*(e(1,:,:,:).^2+e(2,:,:,:).^2+e(3,:,:,:).^2)+e(4,:,:,:).^2+e(5,:,:,:).^2+e(6,:,:,:).^2;
        J3=e(1,:,:,:).*e(2,:,:,:).*e(3,:,:,:)+2.*e(4,:,:,:).*e(5,:,:,:).*e(6,:,:,:)-e(1,:,:,:).*e(5,:,:,:).^2-e(2,:,:,:).*e(6,:,:,:).^2-e(3,:,:,:).*e(4,:,:,:).^2;
        Lode=asin(-1.5*sqrt(3).*J3./(J2).^1.5)/3;
        %--------------------------------------------------------------------------
        value=2*sqrt(J2/3).*sin(Lode-2*pi/3)+em;
    case 'Ex'
        value=Particle_Field.Strain(1,:);
    case 'Ey'
        value=Particle_Field.Strain(2,:);
    case 'Ez'
        value=Particle_Field.Strain(3,:);
    case 'Exy'
        value=Particle_Field.Strain(4,:);
    case 'Eyz'
        value=Particle_Field.Strain(5,:);
    case 'Exz'
        value=Particle_Field.Strain(6,:);
    case 'PEx'
        value=Particle_Field.Plastic_Strain(1,:);
    case 'PEy'
        value=Particle_Field.Plastic_Strain(2,:);
    case 'PEz'
        value=Particle_Field.Plastic_Strain(3,:);
    case 'PExy'
        value=Particle_Field.Plastic_Strain(4,:);
    case 'PEyz'
        value=Particle_Field.Plastic_Strain(5,:);
    case 'PExz'
        value=Particle_Field.Plastic_Strain(6,:);
    case 'PEEQ'
        pstrain=Particle_Field.Pstrain;
        pe=(pstrain(1,:)-pstrain(2,:)).^2+(pstrain(2,:)-pstrain(3,:)).^2+(pstrain(3,:)-pstrain(1,:)).^2+(pstrain(4,:).^2+pstrain(5,:).^2+pstrain(6,:).^2)*1.5;
        value=sqrt(2/9*pe);
    case 'Max_prin_strain'
        e=Particle_Field.Strain;
        e(4:6,:,:,:)=e(4:6,:,:,:)/2;
        em=(e(1,:,:,:)+e(2,:,:,:)+e(3,:,:,:))/3;
        e(1:3,:,:,:)=e(1:3,:,:,:)-em;
        J2=0.5*(e(1,:,:,:).^2+e(2,:,:,:).^2+e(3,:,:,:).^2)+e(4,:,:,:).^2+e(5,:,:,:).^2+e(6,:,:,:).^2;
        J3=e(1,:,:,:).*e(2,:,:,:).*e(3,:,:,:)+2.*e(4,:,:,:).*e(5,:,:,:).*e(6,:,:,:)-e(1,:,:,:).*e(5,:,:,:).^2-e(2,:,:,:).*e(6,:,:,:).^2-e(3,:,:,:).*e(4,:,:,:).^2;
        Lode=asin(-1.5*sqrt(3).*J3./(J2).^1.5)/3;
        %--------------------------------------------------------------------------
        value=2*sqrt(J2/3).*sin(Lode+2*pi/3)+em;
%     case 'Ux'
%         value=Particle_Field.U(1,:);
%     case 'Uy'
%         value=Particle_Field.U(2,:);
%     case 'Uz'
%         value=Particle_Field.U(3,:);
    case 'Vx'
        value=Particle_Field.Velocity(1,:);
    case 'Vy'
        value=Particle_Field.Velocity(2,:);
    case 'Vz'
        value=Particle_Field.Velocity(3,:);
    case 'V'
        value=vecnorm(Particle_Field.Velocity,2,1);
        value=value(:)';
    case 'Fluid_pressure'
        value=Particle_Field.Pressure(1,:);
    otherwise
        error(['Wrong Keyword: ',var]);
end
end