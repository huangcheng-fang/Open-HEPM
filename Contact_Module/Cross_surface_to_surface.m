%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Vertices=Cross_surface_to_surface(coor1,coor2)
%--------------------------------------------------------------------------
mid1=mean(coor1,1);mid2=mean(coor2,1);
rcoor1=coor1-repmat(mid1,size(coor1,1),1);
rcoor2=coor2-repmat(mid2,size(coor2,1),1);
r1=max(sqrt(rcoor1(:,1).^2+rcoor1(:,2).^2));
r2=max(sqrt(rcoor2(:,1).^2+rcoor2(:,2).^2));
if norm(mid1-mid2)-r1-r2>0
    Vertices=zeros(0,2);
    return
end
%--------------------------------------------------------------------------
isin1=inpoly(coor1,coor2);
isin2=inpoly(coor2,coor1);
% isin1=inpolygon(coor1(:,1),coor1(:,2),coor2(:,1),coor2(:,2));
% isin2=inpolygon(coor2(:,1),coor2(:,2),coor1(:,1),coor1(:,2));
Vertices=[coor1(isin1,:);coor2(isin2,:);];
coor2=[coor2;coor2(1,:)];
coor1=[coor1;coor1(1,:)];
for i=1:1:size(coor1,1)-1
    for j=1:1:size(coor2,1)-1
        k=[coor1(i+1,1)-coor1(i,1),coor2(j,1)-coor2(j+1,1);coor1(i+1,2)-coor1(i,2),coor2(j,2)-coor2(j+1,2)];
        if abs(det(k))<1e-10;continue;end
        t=k\[coor2(j,1)-coor1(i,1);coor2(j,2)-coor1(i,2)];
         if (t(1)>=0&&t(1)<=1)&&(t(2)>=0&&t(2)<=1)
             Vertices=[Vertices;(coor1(i+1,1)-coor1(i,1))*t(1)+coor1(i,1),(coor1(i+1,2)-coor1(i,2))*t(1)+coor1(i,2)];
         end
    end
end
Vertices=Truncated_errors(Vertices,12);
Vertices=unique(Vertices,'rows');
vmid=mean(Vertices,1);
Ver=(Vertices-repmat(vmid,size(Vertices,1),1));
[~,loc]=sort(cart2pol(Ver(:,1),Ver(:,2)));
Vertices=Vertices(loc,:);
% hold on
% plot(coor1(:,1),coor1(:,2))
% plot(coor2(:,1),coor2(:,2))
% plot(Vertices(:,1),Vertices(:,2),'*')
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
A=abs(sum(s(1,:)));
flag=Truncated_errors(max(s,[],2).*min(s,[],2)/A,12)>=0;
end

function A=mean(B,C)
A=sum(B,C)/size(B,C);
end