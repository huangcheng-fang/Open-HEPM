%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Intersection=Cross_line_to_volume(line,volume,type,rotation)
line=Truncated_errors(line,12);
volume=Truncated_errors(volume,12);
if rotation
    Normal_vector=line(2,:)-line(1,:);
    R=vector_to_system(Normal_vector/norm(Normal_vector)).';
    mid=mean(line,1);
    line=(line-mid)*R;
    volume=(volume-repmat(mid,size(volume,1),1))*R;
end
%--------------------------------------------------------------------------
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
    t=Lv(1,1)/(Lv(1,1)-Lv(2,1));
    if t>=-1e-8&&t<=1+1e-8
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
Vertices0=Vertices-repmat(mean(Vertices,1),size(Vertices,1),1);
[~,loc]=sort(cart2pol(Vertices0(:,2),Vertices0(:,3)));
Vertices=Vertices(loc,:);
%--------------------------------------------------------------------------
isin=inpoly(line(:,2:3),Vertices(:,2:3));
Intersection=line(isin,:);
Vertices=[Vertices;Vertices(1,:)];
for j=1:1:size(Vertices,1)-1
    t=Vertices(j,2)/(Vertices(j,2)-Vertices(j+1,2));
    if t>=-1e-8&&t<=1+1e-8
        point=Vertices(j+1,:)*t+(1-t)*Vertices(j,:);
        t=(point(3)-line(1,3))/(line(2,3)-line(1,3));
        if t>=-1e-8&&t<=1+1e-8
            Intersection=[Intersection;point];
        end
    end
end
Intersection=Truncated_errors(Intersection,12);
Intersection=unique(Intersection,'rows');
%--------------------------------------------------------------------------
if rotation
    Intersection=Intersection*R.'+repmat(mid,size(Intersection,1),1);
end
end

function flag=inpoly(coor1,coor2)
s=zeros(size(coor1,1),size(coor2,1));
coor2=[coor2;coor2(1,:)];
for i=1:1:size(coor1,1)
    polygon=coor2-repmat(coor1(i,:),size(coor2,1),1);
    for j=1:1:size(polygon,1)-1
        s(i,j)=polygon(j,1)*polygon(j+1,2)-polygon(j+1,1)*polygon(j,2);
    end
end
A=sqrt(abs(sum(s(1,:))));s=sign(s).*sqrt(abs(s));
flag=Truncated_errors(max(s,[],2).*min(s,[],2)/A,10)>=0;
end