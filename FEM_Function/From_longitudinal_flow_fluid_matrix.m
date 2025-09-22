%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [KF,KCP]=From_longitudinal_flow_fluid_matrix(Tri3D,ED_permeability,maxdof)
R=Tri3D.R_matrix;
SDof=Tri3D.dof;
DetJac=Tri3D.detJac;
SB_matrix=Tri3D.B_matrix;
FDof=SDof(3:3:end,:)/3;
FB_matrix=[SB_matrix(1,1:2:end,:,:);SB_matrix(2,2:2:end,:,:);];
EIPnum=size(FB_matrix,3);
Sdofdim=size(SB_matrix,2);
Sdofdim3D=size(SDof,1);
Fdofdim=size(FB_matrix,2);
DetJac=reshape(DetJac,1,1,EIPnum,[]);
%--------------------------------------------------------------------------
[ip,~]=Triangular_integral(1);
N_matrix=Triangular_interpolation(ip);
%===========================Form_stiffness_matrix==========================
EDepJac=pagemtimes(ED_permeability,DetJac);
KValue=pagemtimes(FB_matrix,'transpose',EDepJac,'none');
KValue=pagemtimes(KValue,FB_matrix);
KValue=reshape(sum(KValue,3),Fdofdim,Fdofdim,[]);
%--------------------------------------------------------------------------
KIndexI=reshape(FDof,Fdofdim,1,[]);
KIndexI=repmat(KIndexI,1,Fdofdim,1);
KIndexJ=pagetranspose(KIndexI);
%===========================Form_coupling_matrix===========================
CValue=pagemtimes(DetJac,N_matrix);
CValue=pagemtimes(sum(SB_matrix(1:2,:,:,:),1),'transpose',CValue,'none');
CValue=reshape(sum(CValue,3),Sdofdim,Fdofdim,[]);
CValue=pagemtimes(R,'transpose',CValue,'none');
%--------------------------------------------------------------------------
CIndexI=reshape(SDof,Sdofdim3D,1,[]);
CIndexI=repmat(CIndexI,1,Fdofdim,1);
CIndexJ=reshape(FDof,1,Fdofdim,[]);
CIndexJ=repmat(CIndexJ,Sdofdim3D,1,1);
%===============================Output=====================================
KF=sparse(KIndexI(:),KIndexJ(:),KValue(:),maxdof/3,maxdof/3);
KCP=sparse(CIndexI(:),CIndexJ(:),CValue(:),maxdof,maxdof/3);
end