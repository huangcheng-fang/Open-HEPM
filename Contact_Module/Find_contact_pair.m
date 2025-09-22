%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function [pair,orient]=Find_contact_pair(slave_part,master_part,nodes,surfaces,scenter,sinpart,sinfacet,tolgap)
N=unique(sinfacet);
Patch=struct('inpart',0,'no',0,'surface',0,'center',0,'mid',0);
coder.varsize('Patch');
Patch(1)=[];
for i=1:1:size(N,1)
    loc=find(sinfacet==N(i));
    P.inpart=sinpart(loc(1),1);
    P.no=N(i);
    P.surface=surfaces(loc,:);
    P.center=scenter(loc,:);
    P.mid=[0,0,0];
    P.mid=mean(P.center,1);
    Patch=[Patch;P];
end
%-------------------------------------------------------------------------------
pair=zeros(0,2);coder.varsize('pair');
sortpair=zeros(0,2);coder.varsize('sortpair');
for ii=1:1:numel(Patch)
    if ~ismember(Patch(ii).inpart,slave_part)
        continue
    end
    mid=Patch(ii).mid;
    centeri=Patch(ii).center;
    
    for jj=1:1:numel(Patch)
        if ~ismember(Patch(jj).inpart,master_part)||ii==jj||ismember(sort([ii,jj]),sortpair,'rows')
            continue
        end
        centerj=Patch(jj).center;
        dist=centerj-repmat(mid,size(centerj,1),1);
        distance=dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
        [~,locj]=min(distance);
        dist=centeri-repmat(centerj(locj,:),size(centeri,1),1);
        distance=dist(:,1).^2+dist(:,2).^2+dist(:,3).^2;
        [~,locis]=mink(distance,10);
        
        for kk=1:1:numel(locis)
            loci=locis(kk);
            facei=Patch(ii).surface(loci,:);facei(isnan(facei))=[];
            Nvectori=Get_normal_vector(nodes(facei,1:3));
            facej=Patch(jj).surface(locj,:);facej(isnan(facej))=[];
            Nvectorj=Get_normal_vector(nodes(facej,1:3));

            CS=vector_to_system(Nvectori).';
            f1=nodes(facei,1:3)*CS(:,1:2);
            f2=nodes(facej,1:3)*CS(:,1:2);
            Vertices=Cross_surface_to_surface(f1,f2);
            A=polyarea(f1(:,1),f1(:,2));
            if size(Vertices,1)<3||polyarea(Vertices(:,1),Vertices(:,2))/A<1e-4
                continue
            else
                break
            end
        end
        if size(Vertices,1)<3||polyarea(Vertices(:,1),Vertices(:,2))/A<1e-4
            continue
        end
        d=abs((centeri(loci,:)-centerj(locj,:))*Nvectori')/sqrt(A);
        if d<tolgap%&&Nvectori*Nvectorj'<0
            if Nvectori*Nvectorj'>0
               oriflag=-1;
            else
                oriflag=1;
            end
            pair=[pair;ii*oriflag,jj];
            sortpair=[sortpair;sort([ii,jj])];
        end
    end
end
%-------------------------------------------------------------------------------
if ~isempty(intersect(slave_part(:),master_part(:)))
    for ii=1:1:size(pair,1)
        mid1=Patch(pair(ii,1)).mid;mid2=Patch(pair(ii,2)).mid;
        dist1=sum((Patch(pair(ii,1)).center-repmat(mid1,size(Patch(pair(ii,1)).center,1),1)).^2,2);
        dist2=sum((Patch(pair(ii,2)).center-repmat(mid2,size(Patch(pair(ii,2)).center,1),1)).^2,2);
        [~,loc]=mink(dist1,9);r=dist1(loc(end));
        if sum(dist1<=r)<sum(dist2<=r)
            pair(ii,1:2)=pair(ii,[2,1]);
        end
        pair(ii,1:2)=[Patch(pair(ii,1)).no,Patch(pair(ii,2)).no];
    end
end
%-------------------------------------------------------------------------------
orient=sign(pair(:,1));
pair=abs(pair);
end

function Normal_vector=Get_normal_vector(plane0)
if size(plane0,1)<3
    error('Number of vertices less than 3')
end

if size(plane0,1)>3
    plane=[plane0;plane0(1:2,:)];
else
    Normal_vector=cross(plane0(2,:)-plane0(1,:),plane0(3,:)-plane0(1,:));
    Normal_vector=Normal_vector/norm(Normal_vector);
    return
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
end