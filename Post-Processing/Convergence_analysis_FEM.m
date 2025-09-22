%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function error=Convergence_analysis_FEM(Mesh,Result,analytic,type)
nodes=Mesh.nodes;
elements=Mesh.elements;
Etype=Mesh.etype;
% Estress=Result.Estress;
NU=reshape(Result.Nfield.U,3,[]).';
%==========================================================================
element_num=size(elements,1);
Evolume=zeros(element_num,1);
for ei=1:1:element_num
    nid=elements(ei,1:4);
    Evolume(ei)=det(nodes(nid(2:end),:)-nodes(nid(1),:));
end
Evolume=Evolume/6;
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
        AV=zeros(size(numerial_value,1),size(numerial_value,2));
        for ii=1:1:numel(analytic)
            AV(:,ii)=eval(vectorize(analytic{ii}));
        end
        error=sqrt(sum((AV-numerial_value).^2.*Tet4detJac,'all')/sum((AV).^2.*Tet4detJac,'all'));
    case 'displacement'
        numerial_value=reshape(mean(reshape(NU(nid.',:).',3,4,[]),2),3,[]).';
        AV=zeros(size(numerial_value,1),size(numerial_value,2));
        for ii=1:1:numel(analytic)
            AV(:,ii)=eval(vectorize(analytic{ii}));
        end
        error=sqrt(sum((AV-numerial_value).^2.*Tet4detJac,'all')/sum((AV).^2.*Tet4detJac,'all'));
end
end