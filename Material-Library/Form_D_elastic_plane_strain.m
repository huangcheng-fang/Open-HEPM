%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function D_elastic=Form_D_elastic_plane_strain(elastic_modulus,poisson_ratio)
LAM=elastic_modulus/(1+poisson_ratio);
a=(1-poisson_ratio)/(1-2*poisson_ratio)*LAM;
b=poisson_ratio/(1-2*poisson_ratio)*LAM;
c=0.5*LAM;
D_elastic=[a, b, b, 0;
           b, a, b, 0;
           b, b, a, 0;
           0, 0, 0, c];
end