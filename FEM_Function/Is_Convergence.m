%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function is_convergence=Is_Convergence(Result,Fext,Fint,Flux_ext,Flux_int,constraint,CE,iternum,Tolerance)
dU=Result.dU;
stepU=Result.stepU;
dP=Result.dP;
P=Result.P;
maxdof=size(dU,1);
% constraint=[constraint;CE.Dual_Lagrange.constraint];
%==========================================================================
FHM=[Fext-Fint;Flux_int-Flux_ext];FHM(constraint(:,1))=0;
FHM=FHM-CE.Dual_Lagrange.Fcon-CE.Penalty.vector-CE.interface_storage_vector;
Fun=FHM(1:maxdof);Flux_un=FHM(maxdof+1:end);
%==========================================================================
U_error=norm(dU)/(norm(stepU)+1e-32);
P_error=norm(dP)/(norm(P)+1e-32);
F_error=norm(Fun)/(norm(Fext)+norm(Fint)+1e-32);
Flux_error=norm(Flux_un)/(norm(Flux_ext)+norm(Flux_int)+1e-32);
if norm(Flux_un)<1e-12;Flux_error=0;end
if (U_error<Tolerance.Utol)&&(F_error<Tolerance.Ftol)%&&(P_error<Tolerance.Utol)&&(Flux_error<Tolerance.Ftol)
    is_convergence=true;
else
    is_convergence=false;
end
%===========================Output===========================================
U_error = sprintf('%0.3e',U_error);
P_error = sprintf('%0.3e',P_error);
F_error = sprintf('%0.3e',F_error);
Flux_error = sprintf('%0.3e',Flux_error);
disp(['Displacement (Pore pressure) and Stress (Flux) norm errors in iternum=',num2str(iternum),' are ',U_error,' (',P_error,')',' and ',F_error,' (',Flux_error,')'])
if iternum>Tolerance.MaxIterNum&&~is_convergence
    error('Exceed maximum number of iterations ')
end
end