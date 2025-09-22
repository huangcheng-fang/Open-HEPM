%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [K,D]=Form_USFEM_stiffness_matrix(FE,maxdof)
tic;
[Dof,BT_matrix,DetJac,N_matrix,Dep]=Initial_FE_cell(FE);
NodeDofList2=reshape(1:maxdof/3*6,6,[]);
% load loc
%=========================Form_stiffness_matrix============================
container=Create_parallel_container(numel(FE));
for fi=1:numel(FE)
KValue=pagemtimes(BT_matrix{fi}(:,:,:,:),DetJac{fi});
N=permute(N_matrix{fi},[1,4,3,5,2]);
KValue=pagemtimes(KValue,N);
KValue=permute(sum(KValue,3),[1,2,5,4,3]);
%--------------------------------------------------------------------------
KIndexI=permute(Dof{fi},[1,3,2]);
KIndexJ=reshape(NodeDofList2(:,Dof{fi}(3:3:end,:)/3),[],1,size(KIndexI,3));
KIndexJ=permute(KIndexJ,[2,1,3]);
KIndexI=repmat(KIndexI,1,size(KIndexJ,2),1);
KIndexJ=repmat(KIndexJ,size(KIndexI,1),1,1);
KIndexI=reshape(KIndexI,[],1,size(KValue,3),size(KValue,4));
KIndexJ=reshape(KIndexJ,[],1,size(KValue,3),size(KValue,4));
container(fi,1).c1=reshape(KIndexI(:,:,:,:),[],1);
container(fi,1).c2=reshape(KIndexJ(:,:,:,:),[],1);
container(fi,1).c3=KValue(:);
end
[I,J,V]=Merge_container(container,1);
clear KIndexI KIndexJ KIndexV KValue container
K=sparse(I,J,V,maxdof,maxdof*2);
clear I J V
%==========================Form_volume_vector==============================
volume=zeros(maxdof/3,1);
for fi=1:numel(FE)
H=pagemtimes(N_matrix{fi},'transpose',N_matrix{fi},'none');
H=sum(sum(pagemtimes(H,DetJac{fi}),3),2);
nid=Dof{fi}(3:3:end,:)/3;nid=nid(:);
volume=volume+accumarray(nid,H(:),[maxdof/3,1]);
end
volume=permute(volume,[4,2,3,1]);
%=============================Form_D_matrix================================
container=Create_parallel_container(numel(FE));
for fi=1:numel(FE)
DValue=pagemtimes(Dep{fi},1./volume);
DI=permute(NodeDofList2,[1,3,2]);
DI=repmat(DI,1,size(DI,1));
DJ=pagetranspose(DI);
container(fi,1).c1=DI(:);
container(fi,1).c2=DJ(:);
container(fi,1).c3=DValue(:);
end
[DI,DJ,DV]=Merge_container(container,1);
%=============================Output=======================================
D=sparse(DI,DJ,DV,maxdof*2,maxdof*2);
K=K*D*K';
%==========================================================================
time=toc;fprintf("System stiffness matrix is calculated: %fs\n",time);
end

function [Dof,BT_matrix,DetJac,N_matrix,Dep]=Initial_FE_cell(FE)
numfe=numel(FE);
blank=coder.nullcopy(cell(numfe,1));
Dof=blank;
BT_matrix=blank;
DetJac=blank;
N_matrix=blank;
Dep=blank;
for fi=1:numfe
    Dof{fi}=FE(fi).dof;
    BT_matrix{fi}=reshape(pagetranspose(FE(fi).B_matrix),[],1,size(FE(fi).B_matrix,3),size(FE(fi).B_matrix,4));
    DetJac{fi}=FE(fi).detJac;
    N_matrix{fi}=FE(fi).N_matrix;
    Dep{fi}=FE(fi).Dep;
end
end