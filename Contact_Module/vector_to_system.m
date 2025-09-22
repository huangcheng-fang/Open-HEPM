%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function R=vector_to_system(Normal_vector)
if abs(Normal_vector(1))<0.999847695156391
    Tangential_vector2=cross(Normal_vector,[1,0,0]);
else
    Tangential_vector2=cross(Normal_vector,[0,0,1]);
end
Tangential_vector2=Tangential_vector2/norm(Tangential_vector2);
Tangential_vector1=cross(Tangential_vector2,Normal_vector);
R=[Tangential_vector1;Tangential_vector2;Normal_vector];
end

% function R=vector_to_system(Normal_vector)
% Tangential_vector1=[-Normal_vector(2),Normal_vector([3,1])];
% Tangential_vector2=cross(Normal_vector,Tangential_vector1);
% Tangential_vector2=Tangential_vector2/norm(Tangential_vector2);
% Tangential_vector1=cross(Tangential_vector2,Normal_vector);
% R=[Tangential_vector1;Tangential_vector2;Normal_vector];
% end