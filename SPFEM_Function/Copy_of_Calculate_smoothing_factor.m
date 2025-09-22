%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [Fai,node_volume]=Calculate_smoothing_factor(FEcell,Nactivation,element_num)
CDetJac=FEcell.DetJac;
CElements=FEcell.Elements;
CElement_id=FEcell.Element_id;
CNodes=FEcell.Nodes;
CN_matrix=FEcell.N_matrix;
node_num=size(Nactivation,1);
%=======================Calculate_smoothing_factor=========================
container1=Create_parallel_container(numel(CElement_id));
container2=Create_parallel_container(numel(CElement_id));
for fi=1:numel(CElement_id)
B=FEcell.B_matrix{fi};
B=pagetranspose([B(1,1:3:end,:,:);B(2,2:3:end,:,:);B(3,3:3:end,:,:)]);
Gauss_point=pagemtimes(CN_matrix{fi},CNodes{fi});
dxyz=Gauss_point-CNodes{fi};
dxyz=[B.*dxyz,B];
N=repmat(pagetranspose(CN_matrix{fi}),1,1,1,size(dxyz,4));

KIndexI=pagetranspose(CElements{fi});
KIndexJ=repmat(permute(CElement_id{fi},[4,2,3,1]),size(KIndexI,1),1,1,1);
KValue1=pagemtimes(CDetJac{fi},N);
KValue2=pagemtimes(CDetJac{fi},dxyz);

container1(fi,1).c1=KIndexI(:);
container1(fi,1).c2=KIndexJ(:);
container1(fi,1).c3=reshape(KValue1(:,1,:,:),[],1);

container2(fi,1).c1=reshape(KValue2(:,1,:,:),[],1);
container2(fi,1).c2=reshape(KValue2(:,2,:,:),[],1);
container2(fi,1).c3=reshape(KValue2(:,3,:,:),[],1);
container2(fi,1).c4=reshape(KValue2(:,4,:,:),[],1);
container2(fi,1).c5=reshape(KValue2(:,5,:,:),[],1);
container2(fi,1).c6=reshape(KValue2(:,6,:,:),[],1);
end
[I,J,fai]=Merge_container(container1,1);
[Vxx,Vyy,Vzz,faiX,faiY,faiZ]=Merge_container(container2,1);
%--------------------------------------------------------------------------
fai=sparse(I,J,fai,node_num,element_num);
node_volume=sum(fai,2);
node_volume(~Nactivation)=1e-12;
fai=diag(node_volume)\fai;
r=sqrt(1.2)*(0.75/pi).^(1/3)*(node_volume).^(1/3)*0;
%--------------------------------------------------------------------------
Vxx=accumarray(I,Vxx,[node_num,1]);
Vyy=accumarray(I,Vyy,[node_num,1]);
Vzz=accumarray(I,Vzz,[node_num,1]);
Vxx(~Nactivation)=1e-12;
Vyy(~Nactivation)=1e-12;
Vzz(~Nactivation)=1e-12;

faiX=sparse(I,J,faiX,node_num,element_num);
faiY=sparse(I,J,faiY,node_num,element_num);
faiZ=sparse(I,J,faiZ,node_num,element_num);
dfaidx=diag(sparse(Vxx))\(faiX-diag(sum(faiX,2))*fai);
dfaidy=diag(sparse(Vyy))\(faiY-diag(sum(faiY,2))*fai);
dfaidz=diag(sparse(Vzz))\(faiZ-diag(sum(faiZ,2))*fai);

dfaidx=diag(r)*dfaidx;
dfaidy=diag(r)*dfaidy;
dfaidz=diag(r)*dfaidz;
%================================Output====================================
Fai{1}=fai+dfaidx;
Fai{2}=fai-dfaidx;
Fai{3}=fai+dfaidy;
Fai{4}=fai-dfaidy;
Fai{5}=fai+dfaidz;
Fai{6}=fai-dfaidz;
end