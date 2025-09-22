%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function error=Convergence_analysis_FEM_using_RefResult(Mesh,Result,RefMesh,RefResult,type)
nodes=Mesh.nodes;
elements=Mesh.elements;
Etype=Mesh.etype;
Evolume=Mesh.evolume;
Estress=Result.Estress;
NU=reshape(Result.U,3,[]).';
NP=Result.P;
%==========================================================================
Tet4loc=Etype==334;
Tet4elements=elements(Tet4loc,:);
Tet4detJac=Evolume(Tet4loc);
nid=Tet4elements(:,1:4);
mid=reshape(mean(reshape(nodes(nid.',:).',3,4,[]),2),3,[]).';
x=mid(:,1);y=mid(:,2);z=mid(:,3);
switch type
    case 'stress'
        numerial_value=reshape(Estress(:,:,Tet4loc),6,[]).';
        %----------------------------------------------------------------------
        RefNstress=RefResult.Nstress;
        refx=RefMesh.nodes(:,1);
        refy=RefMesh.nodes(:,2);
        refz=RefMesh.nodes(:,3);
        FS{1}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(1,:,:),[],1));
        FS{2}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(2,:,:),[],1));
        FS{3}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(3,:,:),[],1));
        FS{4}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(4,:,:),[],1));
        FS{5}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(5,:,:),[],1));
        FS{6}=scatteredInterpolant(refx,refy,refz,reshape(RefNstress(6,:,:),[],1));
        AV=zeros(size(numerial_value,1),size(numerial_value,2));
        for ii=1:1:numel(FS)
            AV(:,ii)=FS{ii}(x,y,z);
        end
        %----------------------------------------------------------------------
        error=sqrt(sum((AV-numerial_value).^2.*Tet4detJac,'all')/sum((AV).^2.*Tet4detJac,'all'));
    case 'displacement'
        numerial_value=reshape(mean(reshape(NU(nid.',:).',3,4,[]),2),3,[]).';
        RefNU=reshape(RefResult.U,3,[]).';
        refx=RefMesh.nodes(:,1);
        refy=RefMesh.nodes(:,2);
        refz=RefMesh.nodes(:,3);
        FU{1}=scatteredInterpolant(refx,refy,refz,RefNU(:,1));
        FU{2}=scatteredInterpolant(refx,refy,refz,RefNU(:,2));
        FU{3}=scatteredInterpolant(refx,refy,refz,RefNU(:,3));
        AV=zeros(size(numerial_value,1),size(numerial_value,2));
        for ii=1:1:numel(FU)
            AV(:,ii)=FU{ii}(x,y,z);
        end
        error=sqrt(sum((AV-numerial_value).^2.*Tet4detJac,'all')/sum((AV).^2.*Tet4detJac,'all'));
    case 'pore_pressure'
        numerial_value=mean(NP(nid),2);
        RefNP=RefResult.P;
        refx=RefMesh.nodes(:,1);
        refy=RefMesh.nodes(:,2);
        refz=RefMesh.nodes(:,3);
        FU{1}=scatteredInterpolant(refx,refy,refz,RefNP(:,1));
        AV=zeros(size(numerial_value,1),size(numerial_value,2));
        for ii=1:1:numel(FU)
            AV(:,ii)=FU{ii}(x,y,z);
        end
        error=sqrt(sum((AV-numerial_value).^2.*Tet4detJac,'all')/sum((AV).^2.*Tet4detJac,'all'));
end
end