%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Key=Keywords()
%------------------------------framework----------------------------------
Key.FEM=0;
Key.PFEM=0;
Key.NS_PFEM=1;
Key.ENS_PFEM=2;
Key.SNS_PFEM=3;
%--------------------------------category-------------------------------------
Key.solid='solid';
Key.structure='structure';
%-----------------------------elastic type---------------------------------
Key.linear=1;
Key.porous_elastic=2;
Key.smeared_crack=3;
%-----------------------------section type---------------------------------
Key.none=0;
Key.link=1;
Key.plate_strain=2;
Key.plate_stress=3;
%-----------------------------plastic type---------------------------------
Key.none=0;
Key.mohr_coulomb=1;
Key.modified_cam_clay=2;
Key.mises=3;
Key.drucker_prager=4;
Key.tresca=5;
%-------------------------------fluid type---------------------------------
Key.none=0;
Key.darcy=1;
%---------------------------interface flow type----------------------------
Key.interface_cubic=0;
Key.interface_darcy=1;
%-----------------------------fluid pressure-------------------------------
Key.fluid_pressure=1;
%------------------------------fluid flux----------------------------------
Key.nodal_flux=1;
Key.surface_flux=2;
Key.body_flux=3;
%------------------------------Dirichlet-----------------------------------
Key.displacement=1;
%--------------------------------load type---------------------------------
Key.concentrated_force=1;
Key.surface_traction=2;
Key.pressure=3;
Key.body_force=4;
end