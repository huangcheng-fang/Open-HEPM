%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Flux_couple=Update_Flux_couple(Result,FE,maxdof)
Evstrain=Result.Edvstrain;
Flux_couple=zeros(maxdof/3,1);
for fi=1:1:numel(FE)
    if isempty(FE(fi).D_conductivity);continue;end
    N_matrix=FE(fi).N_matrix;
    DetJac=FE(fi).detJac;
    eid=FE(fi).eid;
    DOF=FE(fi).dof(3:3:end,:)/3;
    IPnum=size(DetJac,3);
    %--------------------------------------------------------------
    Evjac=pagemtimes(DetJac,Evstrain(:,:,1:IPnum,eid));
    Flux_element=sum(pagemtimes(N_matrix,'transpose',Evjac,'none'),3);
    Flux_couple=Flux_couple+accumarray(DOF(:),Flux_element(:),[maxdof/3,1]);
end
end