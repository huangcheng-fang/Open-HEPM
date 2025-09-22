%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
%By Fang Huangcheng @BJTU 
%Email: valy_f@bjtu.edu.cn
%Last update @2023/9/18
%rotation: if surface and volume have been rotated and translated to the midpoint, then rotation=false
function Intersection=Cross_surface_to_volume(surface,volume,type,rotation)
surface=Truncated_errors(surface,12);
volume=Truncated_errors(volume,12);
if rotation
    Normal_vector=Get_normal_vector(surface);
    R=vector_to_system(Normal_vector).';
    mid=mean(surface,1);
    surface=(surface-mid)*R;
    volume=(volume-repmat(mid,size(volume,1),1))*R;
end
switch type
    case 338
        line_Index=[1,2;1,4;1,5;2,3;2,6;3,4;3,7;4,8;5,6;5,8;6,7;7,8];
    case 334
        line_Index=[1,2;1,3;1,4;2,3;2,4;3,4];
    otherwise
        error('Unknown element type')
end
%--------------------------------------------------------------------------
Vertices=zeros(size(line_Index,1),3);
flag=false(size(line_Index,1),1);
for i=1:1:size(line_Index,1)
    Lv=volume(line_Index(i,:),1:3);
    t=Lv(1,3)/(Lv(1,3)-Lv(2,3));
    if t>=-1e-5&&t<=1+1e-5
        Vertices(i,:)=Lv(2,:)*t+(1-t)*Lv(1,:);
        flag(i)=true;
    end
end
Vertices=Vertices(flag,:);
if size(Vertices,1)<3
    Intersection=zeros(0,3);
    return
end
%--------------------------------------------------------------------------
Vertices=Truncated_errors(Vertices,12);
Vertices=unique(Vertices,'rows');
%--------------------------------------------------------------------------
vmid=mean(Vertices,1);
Ver=Vertices-repmat(vmid,size(Vertices,1),1);
[~,loc]=sort(cart2pol(Ver(:,1),Ver(:,2)));
Vertices=Vertices(loc,:);
% fill3(Vertices(:,1),Vertices(:,2),Vertices(:,3),'b')
%--------------------------------------------------------------------------
Intersection=Cross_surface_to_surface(surface(:,1:2),Vertices(:,1:2));
if rotation
    Intersection=[Intersection,zeros(size(Intersection,1),1)]*R.'+repmat(mid,size(Intersection,1),1);
else
    Intersection=[Intersection,zeros(size(Intersection,1),1)];
end
% plot3(Intersection(:,1),Intersection(:,2),Intersection(:,3),'o')
end