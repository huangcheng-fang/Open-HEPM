%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @PolyU
%Email: Huangcheng.Fang@polyu.edu.hk;valy_f@bjtu.edu.cn
%Ref:https://doi.org/10.1016/0045-7825(90)90152-C
function [D_eq,stress,pstrain,voidratio,pc,is_plastic]=Modified_Cam_Clay3D(D_elastic,kappa,lamda,M,stress,pstrain,dstrain,voidratio,pc)
stress=-stress;dstrain=-dstrain;
%--------------------------------------------------------------------------
% voidratio=max(exp(sum(-dstrain(1:3))).*(1+voidratio)-1,0);
%--------------------------------------------------------------------------
p=(stress(1)+stress(2)+stress(3))/3;
if p<0; error('The average stress is not allowed to be tensile stress in Modified_Cam_Clay');end
%--------------------------------------------------------------------------
sigma_tr=stress+D_elastic*dstrain;
[p_tr,q_tr,s_tr]=Get_stress_info(sigma_tr);
F=q_tr^2/M^2+p_tr*(p_tr-pc);
if F/pc^2<1e-3
    D_eq=D_elastic;
    stress=-sigma_tr;
    is_plastic=0;
    return
end
%--------------------------------------------------------------------------
%%here (in theta), i use an updated void ratio. In provided Ref, the void ratio of the previous iteration is used.
theta=(1+voidratio)/(lamda-kappa);
G=D_elastic(6,6);K=D_elastic(1,2)+2*G/3;%%%
fai=0;pc0=pc;p=p_tr;q=q_tr;
iternum=0;itermax=1e5;
while abs(F/pc^2)>1e-10
    dFdfai=Get_dFdfai(fai,K,G,M,theta,p,q,pc);
    fai=fai-F/dFdfai;
    Gcoeff=theta*fai/(1+2*fai*K);
    H=pc0*exp(Gcoeff*(2*p_tr-pc))-pc;
    while abs(H/pc)>1e-10
        dGdpc=-Gcoeff*(H+pc)-1;
        pc=pc-H/dGdpc;
        H=pc0*exp(Gcoeff*(2*p_tr-pc))-pc;
        iternum=iternum+1;
        if iternum>itermax;break;end
    end
    p=(p_tr+fai*K*pc)/(1+2*fai*K);
    q=q_tr/(1+6*G*fai/M^2);
    F=q^2/M^2+p*(p-pc);
    iternum=iternum+1;
    if iternum>itermax;error('Material nonlinear iteration cannot converge in Modified Cam-Clay model');end
end
%--------------------------------------------------------------------------
dFdp=2*p-pc;dFdq=2*q/M^2;s_tr(4:6)=s_tr(4:6)*2;
dF=dFdp*[1/3;1/3;1/3;0;0;0]+dFdq*(1.5)/q_tr*s_tr;
dFdep=-p*pc*theta*ones(6,1);
dpstrain=fai*dF;
stress=stress+D_elastic*(dstrain-dpstrain);
D_eq=D_elastic-0.5*(D_elastic*dF)*(dF'*D_elastic)/(-(dFdep')*dF+dF'*D_elastic*dF);
stress=-stress;dpstrain=-dpstrain;is_plastic=1;
pstrain=pstrain+dpstrain;
end

function [p,q,s]=Get_stress_info(s)
p=(s(1)+s(2)+s(3))/3;
s(1)=s(1)-p;s(2)=s(2)-p;s(3)=s(3)-p;
J2=0.5*(s(1).^2+s(2).^2+s(3).^2)+s(4).^2+s(5).^2+s(6).^2;
q=sqrt(3*J2);
end

function dFdfai=Get_dFdfai(fai,K,G,M,theta,p,q,pc)
dFdp=2*p-pc;dFdq=2*q/M^2;dFdpc=-p;
temp=(1+(2*K+theta*pc)*fai);
dpdfai=-K*dFdp/temp;
dqdfai=-q/(fai+M^2/(6*G));
dpcdfai=theta*pc*dFdp/temp;
dFdfai=dFdp*dpdfai+dFdq*dqdfai+dFdpc*dpcdfai;
end