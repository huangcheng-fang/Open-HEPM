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
%elasPara(D_elastic):6x6 matrix
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
function [Dct,stress,pstrain,histPara]=Mohr_Coulomb3D(elasPara,plasPara,hardPara,histPara,stress,pstrain,dstrain)
D_elastic=elasPara;fai=plasPara(1);pfai=plasPara(2);m=0;
%--------------------------------------------------------------------------
theta_T=25/180*pi;
A1=[-2.93057555085368	 8.48875837836269	-4.67585018301484;
    -2.93057555085368	-8.48875837836269	-4.67585018301484];
A2=[-3.93747122467738	 8.32143144099294	-4.65632790876395;
     3.93747122467738	 8.32143144099294	 4.65632790876395];
A=A1+A2*sin(fai);
%===========================Elastic prediction=============================
dsigma=D_elastic*dstrain;
[~,sigmam,sigmabar,~,K,~,~,~]=Get_stress_info(stress+dsigma,fai,theta_T,A);
[c,~]=Hardening(hardPara,Get_equivalent_strain(pstrain));
F=sigmam*sin(fai)+sqrt((sigmabar*K)^2+(m*c*cos(fai))^2)-c*cos(fai);
refF=abs(F+c*cos(fai))+c*cos(fai);
%=============================Elastic output===============================
if F<refF*1e-5
    Dct=D_elastic;
    stress=stress+dsigma;
    return;
end
%===========================Plastic correction=============================
tol=1e-5;eye6=eye(6);dLamda=0;t=0;dt=1;flag=true;dsigma=zeros(6,1);dpstrain=zeros(6,1);
dpstrain_stored=zeros(6,1);dLamda_stored=0;c_stored=c;k=1;
while true
if flag
    if t==1;break;end
    t=t+dt;
    dpstrain_stored=dpstrain;
    dLamda_stored=dLamda;
    c_stored=c;
else
    dt=dt/2;
    t=t-dt;
    dpstrain=dpstrain_stored;
    dLamda=dLamda_stored;
    c=c_stored;
    if dt<0.1
        error('Material nonlinear iteration cannot converge in Mohr_Coulomb')
    end
end
dsigma=D_elastic*(t*dstrain-dpstrain);
flag=true;iternum=0;
while true
    %Calculate yield value
    [s,sigmam,sigmabar,Lode,K,dKdtheta,ddKddtheta,J3]=Get_stress_info(stress+dsigma,fai,theta_T,A);
    F=sigmam*sin(fai)+sqrt((sigmabar*K)^2+(m*c*cos(fai))^2)-c*cos(fai);
    refF=abs(F+c*cos(fai))+c*cos(fai);

    %Calculate derivative
    dFdSig=Get_yield_derivative(fai,c,s,sigmabar,Lode,K,m,dKdtheta);
    [dGdSig,ddGddSig]=Get_potential_derivative(pfai,c,s,sigmabar,Lode,K,m,dKdtheta,ddKddtheta,J3);
    dpstrain=dLamda*dGdSig;
    
    % Hardening
    [c_corrected,H]=Hardening(hardPara,Get_equivalent_strain(pstrain+dpstrain));
    
    % Calculate residuals
    sigma_res=D_elastic*(t*dstrain-dpstrain)-dsigma;
    c_res=c_corrected-c;

    % Calculate tangent matrix
    dpeqdLamda=Get_equivalent_strain(dGdSig);
    dFdc=c*(m*cos(fai))^2/sqrt((sigmabar*K)^2+(m*c*cos(fai))^2)-cos(fai);
    k=[eye6+dLamda*D_elastic*ddGddSig,  D_elastic*dGdSig,         [0;0;0;0;0;0];
                             dFdSig.',                 0,                  dFdc;
                          0,0,0,0,0,0,     -H*dpeqdLamda,                    1];

    % Convergence criterion
    if abs(F)/refF<=tol&&norm(sigma_res)/norm(dsigma)<=tol&&abs(c_res)/abs(c_corrected)<=tol;break;end
    
    %Update variables
    x=k\[sigma_res;-F;c_res];
    dsigma=dsigma+x(1:6);
    dLamda=dLamda+x(7);
    c=c+x(8);
    
    %count
    iternum=iternum+1;
    if iternum>200
       flag=false;
       break
    end
end
end
%========================Plastic output====================================
stress=stress+dsigma;
pstrain=pstrain+dpstrain;
Dct=k\[D_elastic;zeros(2,6)];%Consistent tangent matrix
Dct=Dct(1:6,1:6);
% Dc=(eye6+dLamda*D_elastic*ddGddSig)\D_elastic;
% Dct=Dc-(Dc*dGdSig)*(dFdSig'*Dc)/(H*cos(fai)*dpeqdLamda+dFdSig'*Dc*dGdSig);
end

function [s,sm,sbar,Lode,K,dKdtheta,ddKddtheta,J3]=Get_stress_info(s,fai,theta_T,A)
sm=(s(1)+s(2)+s(3))/3;
s(1)=s(1)-sm;s(2)=s(2)-sm;s(3)=s(3)-sm;
J2=0.5*(s(1).^2+s(2).^2+s(3).^2)+s(4).^2+s(5).^2+s(6).^2;
J3=s(1)*s(2)*s(3)+2*s(4)*s(5)*s(6)-s(1)*s(5).^2-s(2)*s(6).^2-s(3)*s(4).^2;
sbar=sqrt(J2);
if sbar~=0
    temp=-1.5*sqrt(3.0)*J3/sbar^3;
    temp=sign(temp)*min(abs(temp),1.0);
    Lode=real(asin(temp)/3.0);
else
    Lode=0;sbar=1e-32;
end

if Lode>theta_T
    K=A(1,1)+A(1,2)*sin(3.0*Lode)+A(1,3)*sin(3.0*Lode).^2;
    dKdtheta=3.0*A(1,2)*cos(3.0*Lode)+A(1,3)*3*sin(6.0*Lode);
    ddKddtheta=-9.0*A(1,2)*sin(3.0*Lode)+A(1,3)*18*cos(6.0*Lode);
elseif Lode<-theta_T
    K=A(2,1)+A(2,2)*sin(3.0*Lode)+A(2,3)*sin(3.0*Lode).^2;
    dKdtheta=3.0*A(2,2)*cos(3.0*Lode)+A(2,3)*3*sin(6.0*Lode);
    ddKddtheta=-9.0*A(2,2)*sin(3.0*Lode)+A(2,3)*18*cos(6.0*Lode);
else
    K=cos(Lode)-1.0/sqrt(3.0)*sin(fai)*sin(Lode);
    dKdtheta=-sin(Lode)-1.0/sqrt(3.0)*sin(fai)*cos(Lode);
    ddKddtheta=-cos(Lode)+1.0/sqrt(3.0)*sin(fai)*sin(Lode);
end
end

function dF=Get_yield_derivative(fai,c,s,sigmabar,Lode,K,m,dKdtheta)
J2=sigmabar^2;
dsigmam=[1/3;1/3;1/3;0.0;0.0;0.0];
dJ2=[s(1);s(2);s(3);2*s(4);2*s(5);2*s(6);];
dJ3=[s(2)*s(3)-s(5)*s(5)+J2/3;
    s(1)*s(3)-s(6)*s(6)+J2/3;
    s(1)*s(2)-s(4)*s(4)+J2/3;
    (s(5)*s(6)-s(3)*s(4))*2;
    (s(6)*s(4)-s(1)*s(5))*2;
    (s(4)*s(5)-s(2)*s(6))*2;];
dsigmabar=dJ2/(2*sigmabar);
%==========================================================================
C1=sin(fai);
C2=K-tan(3.0*Lode)*dKdtheta;
C3=-sqrt(3.0)/2.0/cos(3.0*Lode)/J2*dKdtheta;
tempA=sqrt((sigmabar*K)^2+(m*c*cos(fai))^2);
A=sigmabar*K/tempA;
dF=C1*dsigmam+A*(C2*dsigmabar+C3*dJ3);
end

function [dG,ddG]=Get_potential_derivative(pfai,c,s,sigmabar,Lode,K,m,dKdtheta,ddKddtheta,J3)
J2=sigmabar^2;
dsigmam=[1/3;1/3;1/3;0.0;0.0;0.0];
dJ2=[s(1);s(2);s(3);2*s(4);2*s(5);2*s(6);];
dJ3=[s(2)*s(3)-s(5)*s(5)+J2/3;
    s(1)*s(3)-s(6)*s(6)+J2/3;
    s(1)*s(2)-s(4)*s(4)+J2/3;
    (s(5)*s(6)-s(3)*s(4))*2;
    (s(6)*s(4)-s(1)*s(5))*2;
    (s(4)*s(5)-s(2)*s(6))*2;];
dsigmabar=dJ2/(2*sigmabar);
ddJ2=[ 2/3, -1/3, -1/3, 0, 0, 0;
      -1/3,  2/3, -1/3, 0, 0, 0;
      -1/3, -1/3,  2/3, 0, 0, 0;
         0,    0,    0, 2, 0, 0;
         0,    0,    0, 0, 2, 0;
         0,    0,    0, 0, 0, 2];
ddJ3=[s(1)-s(2)-s(3),         2*s(3),          2*s(2),  2*s(4), -4*s(5),  2*s(6);
              2*s(3), s(2)-s(1)-s(3),          2*s(1),  2*s(4),  2*s(5), -4*s(6);
              2*s(2),         2*s(1),  s(3)-s(1)-s(2), -4*s(4),  2*s(5),  2*s(6);
              2*s(4),         2*s(4),         -4*s(4), -6*s(3),  6*s(6),  6*s(5);
             -4*s(5),         2*s(5),          2*s(5),  6*s(6), -6*s(1),  6*s(4);
              2*s(6),        -4*s(6),          2*s(6),  6*s(5),  6*s(4), -6*s(2);];
ddJ3=ddJ3/3;
ddsigmabar=1/(2*sigmabar)*ddJ2-dJ2*dJ2.'*(0.25*sigmabar^(-3));
dtheta=-sqrt(3)/2/sigmabar^3/cos(3*Lode)*(dJ3-3*J3/sigmabar*dsigmabar);
%==========================================================================
C1=sin(pfai);
C2=K-tan(3.0*Lode)*dKdtheta;
C3=-sqrt(3.0)/2.0/cos(3.0*Lode)/J2*dKdtheta;
tempA=sqrt((sigmabar*K)^2+(m*c*cos(pfai))^2);
A=sigmabar*K/tempA;
tempdG=A*(C2*dsigmabar+C3*dJ3);
dG=C1*dsigmam+tempdG;
%==========================================================================
dA=(1-A^2)/tempA*(dsigmabar*K+sigmabar*dKdtheta*dtheta);
dC2=dtheta*(dKdtheta-ddKddtheta*tan(3*Lode)-3*dKdtheta*sec(3*Lode)^2);
dC3=sqrt(3)/2/cos(3*Lode)/sigmabar^2*(-dtheta*(ddKddtheta+3*dKdtheta*tan(3*Lode))+2/sigmabar*dKdtheta*dsigmabar);
ddG=dC2*dsigmabar'+C2*ddsigmabar+dC3*dJ3'+C3*ddJ3;
ddG=tempdG*dA'+A*ddG;
end

function pe=Get_equivalent_strain(strain)
pe=(strain(1)-strain(2)).^2+(strain(2)-strain(3)).^2+(strain(3)-strain(1)).^2+(strain(4).^2+strain(5).^2+strain(6).^2)*1.5;
pe=sqrt(2/9*pe);
end

function [sigmarlim,H]=Hardening(hardpara,pe)
[~,loc]=min(abs(hardpara(:,1)-pe));
if hardpara(loc,1)-pe>0
    loc=loc-1;
end
H=hardpara(loc,3);sigmarlim=hardpara(loc,2)+H*(pe-hardpara(loc,1));
end


% theta_T=-25/180*pi;fai=30/180*pi;
% k=[1,sin(3.0*theta_T),sin(3.0*theta_T).^2;
%    0,3.0*cos(3.0*theta_T),3*sin(6.0*theta_T);
%    0,-9.0*sin(3.0*theta_T),18*cos(6.0*theta_T)];
% F=[cos(theta_T),-1.0/sqrt(3.0)*sin(theta_T);
%    -sin(theta_T),-1.0/sqrt(3.0)*cos(theta_T);
%    -cos(theta_T),+1.0/sqrt(3.0)*sin(theta_T);];
% x=(k\F)';A=x(1);B=x(2);C=x(3);