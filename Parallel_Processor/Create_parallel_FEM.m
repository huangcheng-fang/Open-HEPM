%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function parFE=Create_parallel_FEM(FE)
num=100000;num0=num-1;
parnum=zeros(numel(FE),1);
for fi=1:1:numel(FE)
    enum=numel(FE(fi).eid);
    parnum(fi)=ceil(enum/num);
end

parFE(sum(parnum,'all'),1).type=zeros(enum*0,enum*0);
pari=0;
for fi=1:1:numel(FE)
    for jj=1:1:parnum(fi)-1
        pari=pari+1;
        parFE(pari).type=FE(fi).type;
        parFE(pari).R_matrix=FE(fi).R_matrix(:,:,:,jj*num-num0:jj*num);
        parFE(pari).B_matrix=FE(fi).B_matrix(:,:,:,jj*num-num0:jj*num);
        parFE(pari).detJac=FE(fi).detJac(:,:,:,jj*num-num0:jj*num);
        parFE(pari).dof=FE(fi).dof(:,jj*num-num0:jj*num);
        parFE(pari).eid=FE(fi).eid(jj*num-num0:jj*num,1);
        parFE(pari).N_matrix=FE(fi).N_matrix;
        parFE(pari).material=FE(fi).material;
    end
    pari=pari+1;
    parFE(pari).type=FE(fi).type;
    parFE(pari).B_matrix=FE(fi).B_matrix(:,:,:,parnum(fi)*num-num0:end);
    parFE(pari).R_matrix=FE(fi).R_matrix(:,:,:,parnum(fi)*num-num0:end);
    parFE(pari).detJac=FE(fi).detJac(:,:,:,parnum(fi)*num-num0:end);
    parFE(pari).dof=FE(fi).dof(:,parnum(fi)*num-num0:end);
    parFE(pari).eid=FE(fi).eid(parnum(fi)*num-num0:end,1);
    parFE(pari).N_matrix=FE(fi).N_matrix;
    parFE(pari).material=FE(fi).material;
end
end