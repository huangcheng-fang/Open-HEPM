%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Normal_vector=Get_normal_vector(plane0)
if size(plane0,1)<3
    error('Number of vertices less than 3')
end

if size(plane0,1)>3
    plane=[plane0;plane0(1:2,:)];
else
    plane=plane0;
end
Normal_vectors=zeros(size(plane,1)-2,3);
for i=1:1:size(plane,1)-2
    line1=plane(i+1,:)-plane(i,:);
    line2=plane(i+2,:)-plane(i+1,:);
    Normal_vectors(i,:)=cross(line1,line2);
    Normal_vectors(i,:)=Normal_vectors(i,:)/norm(Normal_vectors(i,:));
end
Normal_vector=sum(Normal_vectors,1)/size(Normal_vectors,1);
Normal_vector=Normal_vector/norm(Normal_vector);
%--------------------------------------------------------------------------
end