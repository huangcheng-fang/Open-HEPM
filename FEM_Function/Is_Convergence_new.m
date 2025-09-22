%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function is_convergence=Is_Convergence_new(dU,stepU,Fext,Fint,constraint,CE,iternum,Tolerance)
maxdof=size(dU,1);
%==========================================================================
Fint=Fint+CE.Dual_Lagrange.Fcon+CE.Penalty.vector;
Fun=Fext-Fint;Fun(constraint(:,1))=0;
Fun=Fun(1:maxdof);
%==========================================================================
U_error=norm(dU)/(norm(stepU)+1e-32);
F_error=norm(Fun)/(norm(Fext)+norm(Fint)+1e-32);
if (U_error<Tolerance.Utol)&&(F_error<Tolerance.Ftol)%&&(P_error<Tolerance.Utol)&&(Flux_error<Tolerance.Ftol)
    is_convergence=true;
else
    is_convergence=false;
end
%===========================Output===========================================
U_error = sprintf('%0.3e',U_error);
F_error = sprintf('%0.3e',F_error);
disp(['Displacement and Stress norm errors in iternum=',num2str(iternum),' are ',U_error,' and ',F_error])
if iternum>Tolerance.MaxIterNum&&~is_convergence
    is_convergence=true;
    warning('Exceed maximum number of iterations ')
end
end