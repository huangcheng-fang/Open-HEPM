%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function isin=In_polyhedron(polyhedron_nodes,polyhedron,points)
Nconnection=Get_node_connection(polyhedron);
%==========================surface normal==================================
line1=polyhedron_nodes(polyhedron(:,2),:)-polyhedron_nodes(polyhedron(:,1),:);
line2=polyhedron_nodes(polyhedron(:,3),:)-polyhedron_nodes(polyhedron(:,1),:);
surface_normal=cross(line1,line2,2);
surface_area_sqrt=sqrt(sum(surface_normal.^2,2));
surface_normal=surface_area_sqrt.\surface_normal;
surface_area_sqrt=sqrt(surface_area_sqrt);
%=======================find nearest surface node==========================
loc=true(size(polyhedron_nodes,1),1);
loc(polyhedron)=false;
polyhedron_nodes(loc,:)=nan;
kdTree=KDTreeSearcher(polyhedron_nodes);
s_nids=knnsearch(kdTree,points,'k',5);
%====================remove outer element==================================
isin=false(size(points,1),1);
for iii=1:1:size(s_nids,2)
    for ei=1:1:size(points,1)
        flag=0;
        s_nid=s_nids(ei,iii);
        sid=unique(Nconnection(s_nid,1:end-1));
        sid(sid==0)=[];
        snid=polyhedron(sid(1),:);
        ref_point=surface_normal(sid(1),:)*(1e-8*surface_area_sqrt(sid(1)))+mean(polyhedron_nodes(snid,:),1)-polyhedron_nodes(s_nid,:);
        emid_point=points(ei,:)-polyhedron_nodes(s_nid,:);
        for i=1:1:numel(sid)
            d1=(surface_normal(sid(i),:)*ref_point.');
            d2=(surface_normal(sid(i),:)*emid_point.');

            if abs(d2)/surface_area_sqrt(sid(1))<1e-12
                flag=0;break
            end
            if d1*d2>=0
                continue
            end
            d1=abs(d1);d2=abs(d2);

            intersection=(ref_point*d2+emid_point*d1)/(d1+d2);
            snid=polyhedron(sid(i),:);
            loc=snid~=s_nid;
            surface_line=polyhedron_nodes(snid(loc),:)-polyhedron_nodes(s_nid,:);

            surface_line=sqrt(sum(surface_line.^2,2)).\surface_line;
            intersection=norm(intersection).\intersection;
            if all(surface_line(1,:)*surface_line(2,:).'<surface_line*intersection')
                flag=flag+1;
            end
        end
        if mod(flag,2)==0
            isin(ei)=true;
        end
    end
end
end