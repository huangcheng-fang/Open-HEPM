%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @BJTU
%Implicit return algorithm; Backward Euler Method
%Last update @2024/1/24
%Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
%Ref:KRISTIAN KRABBENHØFT, 2002, BASIC COMPUTATIONAL PLASTICITY
%==========================Input parameter=================================
%elasPara:6x6 matrix,[D_elastic]
%plasPara:2x1 vector, [friction angle; dilation angle]; Unit: rad
%hardPara:nx3 vector,[Equivalent plastic strain,Cohesion,tangent];%%perfectly-plastic:[0,c,0]
%histPara: [];
%stress:6x1 vector,[Sx;Sy;Sz;Sxy;Syz;Szx]
%pstrain(plastic strain):6x1 vector,[PEx;PEy;PEz;PExy;PEyz;PEzx]
%dstrain(strain increment):6x1 vector,[dEx;dEy;dEz;dExy;dEyz;dEzx]
%==========================Output parameter================================
%Dct(Consistent tangent matrix):6x6 matrix
%stress:6x1 vector,[Sx;Sy;Sz;Sxy;Syz;Szx]
%pstrain(plastic strain):6x1 vector,[PEx;PEy;PEz;PExy;PEyz;PEzx]
%histPara: [];
function [Dct,stress,pstrain,histPara]=Von_Mises_new(elasPara,plasPara,hardPara,histPara,stress,pstrain,dstrain)
D_elastic=elasPara;
[sigmarlim,~]=Hardening(hardPara,pstrain);
%===========================Elastic prediction=============================
dsigma=D_elastic*dstrain;
[s,~,sigmabar]=Get_stress_info(stress+dsigma);
%=============================Elastic output===============================
F=sigmabar-sigmarlim;refF=sigmarlim;
if F<refF*1e-5||norm(dstrain)==0
    Dct=D_elastic;
    stress=stress+dsigma;
    return;
end
%======================Plastic constitutive integral=======================
iternum=0;tol=1e-10;eye6=eye(6);dLamda=0;
while true
    %Calculate derivative
    [dGdSig,ddGddSig]=Get_potential_derivative(s,sigmabar);
    dpstrain=dLamda*dGdSig;

    % Calculate residuals
    dsigma_corrected=D_elastic*(dstrain-dpstrain);
    sigma_res=dsigma_corrected-dsigma;
    [sigmarlim_corrected,H]=Hardening(hardPara,pstrain+dpstrain);
    sigmarlim_res=sigmarlim_corrected-sigmarlim;  
    % Calculate tangent matrix
    dpeqdLamda=Get_equivalent_strain(dGdSig);
    k=[eye6+dLamda*D_elastic*ddGddSig,  D_elastic*dGdSig,         [0;0;0;0;0;0];
                             dGdSig.',                 0,                    -1;
                          0,0,0,0,0,0,     -H*dpeqdLamda,                    1];
    
    % Convergence criterion
    if abs(F)<=refF*tol&&norm(sigma_res)<=tol*norm(dsigma)&&abs(sigmarlim_res)<=tol*abs(sigmarlim_corrected); break;end

    %Update variables
    x=k\[sigma_res;-F;sigmarlim_res];
    dsigma=dsigma+x(1:6);
    dLamda=dLamda+x(7);
    sigmarlim=sigmarlim+x(8);

    %Calculate yield value
    [s,~,sigmabar]=Get_stress_info(stress+dsigma);
    F=sigmabar-sigmarlim;refF=sigmarlim;

    %count
    iternum=iternum+1;
    if iternum>10000
        error('Material nonlinear iteration cannot converge in UMat')
    end
end
%========================Plastic output====================================
stress=stress+dsigma;
pstrain=pstrain+dpstrain;
Dct=k\[D_elastic;zeros(2,6)];%Consistent tangent matrix
Dct=Dct(1:6,1:6);
end

function [s,sm,sbar]=Get_stress_info(s)
sm=(s(1)+s(2)+s(3))/3;
s(1)=s(1)-sm;s(2)=s(2)-sm;s(3)=s(3)-sm;
J2=0.5*(s(1).^2+s(2).^2+s(3).^2)+s(4).^2+s(5).^2+s(6).^2;
sbar=sqrt(3*J2);
end

function [sigmarlim,H]=Hardening(hardpara,pstrain)
pe=Get_equivalent_strain(pstrain);
[~,loc]=min(abs(hardpara(:,1)-pe));
if hardpara(loc,1)-pe>0
    loc=loc-1;
end
H=hardpara(loc,3);sigmarlim=hardpara(loc,2)+H*(pe-hardpara(loc,1));
end

function [dGdSig,ddGddSig]=Get_potential_derivative(s,sigmabar)
dJ2dsigma=[s(1);s(2);s(3);2*s(4);2*s(5);2*s(6)];
dsdsigma=[ 2/3, -1/3, -1/3, 0, 0, 0;
          -1/3,  2/3, -1/3, 0, 0, 0;
          -1/3, -1/3,  2/3, 0, 0, 0;
             0,    0,    0, 2, 0, 0;
             0,    0,    0, 0, 2, 0;
             0,    0,    0, 0, 0, 2];
dGdSig=[s(1);s(2);s(3);2*s(4);2*s(5);2*s(6)]/(2/3*sigmabar);
ddGddSig=(-2.25*(sigmabar)^-3*dJ2dsigma)*dJ2dsigma.'+1.5/sigmabar*dsdsigma;
end

function pe=Get_equivalent_strain(strain)
pe=(strain(1)-strain(2)).^2+(strain(2)-strain(3)).^2+(strain(3)-strain(1)).^2+(strain(4).^2+strain(5).^2+strain(6).^2)*1.5;
pe=sqrt(2/9*pe);
end

