%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function K=From_longitudinal_flow_solid_matrix(Tri3D,EDe,maxdof)
R=Tri3D.R_matrix;
SDof=Tri3D.dof;
DetJac=Tri3D.detJac;
SB_matrix=Tri3D.B_matrix;
EIPnum=size(SB_matrix,3);
Sdofdim=size(SB_matrix,2);
Sdofdim3D=size(SDof,1);
DetJac=reshape(DetJac,1,1,EIPnum,[]);
%===========================Form_stiffness_matrix==========================
EDepJac=pagemtimes(EDe,DetJac);
KValue=pagemtimes(SB_matrix,'transpose',EDepJac,'none');
KValue=pagemtimes(KValue,SB_matrix);
KValue=reshape(sum(KValue,3),Sdofdim,Sdofdim,[]);
KValue=pagemtimes(KValue,R);
KValue=pagemtimes(R,'transpose',KValue,'none');
%--------------------------------------------------------------------------
KIndexI=reshape(SDof,Sdofdim3D,1,[]);
KIndexI=repmat(KIndexI,1,Sdofdim3D,1);
KIndexJ=pagetranspose(KIndexI);
%===============================Output=====================================
K=sparse(KIndexI(:),KIndexJ(:),KValue(:),maxdof,maxdof);
end