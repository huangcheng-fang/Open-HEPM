%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Lpoints=Line_grid_intersect(lines,nodes,node_cubeID)
if isempty(node_cubeID)
    node_cubeID=ceil(nodes);
end
%==========================================================================
Lpoints=zeros(1e6,3);flag=1;
for li=1:1:size(lines,1)
    line_coor_li=nodes(lines(li,:),:);
    cut_cellID=node_cubeID(lines(li,:),:);
    ID_range=[min(cut_cellID,[],1);max(cut_cellID,[],1)];
    for api=1:1:3
        for ii=ID_range(1,api):1:ID_range(2,api)-1
            local_coor=(ii-line_coor_li(1,api))/(line_coor_li(2,api)-line_coor_li(1,api));
            if local_coor<=1.000001&&local_coor>=-0.000001
                Lpoint=(line_coor_li(2,:)-line_coor_li(1,:))*local_coor+line_coor_li(1,:);
                Lpoint(api)=ii;
                Lpoints(flag,:)=Lpoint;
                flag=flag+1;
            end
        end
    end
end
if flag>1e6
    error('Exceeding the preset array size');
end
Lpoints(flag:end,:)=[];
end