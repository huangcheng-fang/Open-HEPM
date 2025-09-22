%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function NResult=Mapping_variable_E2N(Result)
Nnum=size(Result.U,1)/3;
FEinfo=Result.FEinfo;
Estress=Result.Estress;
Estrain=Result.Estrain;
Epstrain=Result.Epstrain;
%==========================================================================
container=Create_parallel_container(numel(FEinfo));
for fi=1:1:numel(FEinfo)
    eid=FEinfo(fi).eid;
    nid=FEinfo(fi).dof(3:3:end,:)/3;
    detJac=FEinfo(fi).detJac;
    N_matrix=FEinfo(fi).N_matrix;
    IPnum=size(FEinfo(fi).detJac,3);
    %----------------------------------------------------------------------
    N_matrix=permute(pagemtimes(N_matrix,detJac),[3,2,4,1]);
    stress=permute(Estress(:,:,1:IPnum,eid),[1,3,4,2]);
    strain=permute(Estrain(:,:,1:IPnum,eid),[1,3,4,2]);
    pstrain=permute(Epstrain(:,:,1:IPnum,eid),[1,3,4,2]);
    %----------------------------------------------------------------------
    stress=pagemtimes(stress,N_matrix);
    strain=pagemtimes(strain,N_matrix);
    pstrain=pagemtimes(pstrain,N_matrix);
    volume=sum(N_matrix,1);
    IndexJ=repmat(permute(nid,[3,1,2]),6,1,1);
    IndexI=repmat((1:6)',1,size(IndexJ,2),size(IndexJ,3));
    %----------------------------------------------------------------------
    container(fi,1).c1=IndexI(:);
    container(fi,1).c2=IndexJ(:);
    container(fi,1).c3=stress(:);
    container(fi,1).c4=strain(:);
    container(fi,1).c5=pstrain(:);
    container(fi,1).c6=volume(:);
end
[I,J,stress,strain,pstrain,volume]=Merge_container(container,1);
%============================Normalization=================================
NResult.U=Result.U;NResult.P=Result.P;NResult.NCstress=Result.NCstress;
if isempty(I)
    NResult.Nstress=zeros([6,Nnum]);
    NResult.Nstrain=zeros([6,Nnum]);
    NResult.Npstrain=zeros([6,Nnum]);
else
    Nvolume=accumarray([I(6:6:end)/6,J(6:6:end)],volume,[1,Nnum]);
    NResult.Nstress=accumarray([I,J],stress,[6,Nnum])./Nvolume;
    NResult.Nstrain=accumarray([I,J],strain,[6,Nnum])./Nvolume;
    NResult.Npstrain=accumarray([I,J],pstrain,[6,Nnum])./Nvolume;
end
end