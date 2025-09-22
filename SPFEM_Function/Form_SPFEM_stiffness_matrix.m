%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function K=Form_SPFEM_stiffness_matrix(FEcell,FEUdof,Dmatrix,Fai,node_volume,maxdof)
tic
CB_matrix=FEcell.B_matrix;
CEid=FEcell.Element_id;
FEcell_num=numel(CEid);
%=========================Calculate_B_matrix===============================
container=Create_parallel_container(FEcell_num);
for fi=1:FEcell_num
KValue=CB_matrix{fi};
KIndexI=reshape((CEid{fi}*6-[5,4,3,2,1,0])',6,1,1,[]);
KIndexI=repmat(KIndexI,1,size(CB_matrix{fi},2),1,1);
KIndexJ=reshape(FEUdof{fi},size(FEUdof{fi},1),1,1,[]);
KIndexJ=pagetranspose(repmat(KIndexJ,1,size(KIndexI,1),1,1));
container(fi,1).c1=KIndexI(:);
container(fi,1).c2=KIndexJ(:);
container(fi,1).c3=KValue(:);
end
[I,J,V]=Merge_container(container,1);
%--------------------------------------------------------------------------
B=sparse(I,J,V,size(Fai{1},2)*6,maxdof);
%=============================Form_D_matrix================================
node_volume=node_volume/6;
NodeDofList2=reshape(1:maxdof/3*6,6,[]);
DV=pagemtimes(Dmatrix,permute(full(node_volume),[4,2,3,1]));
DI=permute(NodeDofList2,[1,3,2]);
DI=repmat(DI,1,size(DI,1));
DJ=pagetranspose(DI);
%====================Calculate_smoothing_factor_matrix=====================
K=sparse(maxdof,maxdof);
for ipi=1:1:numel(Fai)
    [I,J,V]=find(Fai{ipi});
    I=I*6;J=J*6;
    fai=sparse(I,J,V,size(Fai{ipi},1)*6,size(Fai{ipi},2)*6);
    for ii=1:1:5
        fai=fai+sparse(I-ii,J-ii,V,size(fai,1),size(fai,2));
    end
    SB=fai*B;
    D=sparse(DI(:),DJ(:),reshape(DV(:,:,ipi,:),[],1),maxdof*2,maxdof*2);
    K=K+SB'*D*SB;
end
end
