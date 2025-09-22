%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Adjust_contact_surface(Model,varargin)
%==========================Check input=====================================
tolerance=1e-8;
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'Slave_set'
            slave_set=varargin{vi+1};
        case 'Master_set'
            master_set=varargin{vi+1};
        case 'Tolerance'
            tolerance=varargin{vi+1};
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
nodes=Model.Mesh.nodes;
Sset=Model.Set.surface_set;
surfaces=Model.Mesh.surfaces;
slave=Merge_cell(Sset,slave_set);
master=Merge_cell(Sset,master_set);
slave_surface=surfaces(slave,:);
master_surface=surfaces(master,:);
err=1;flag=0;
while err>tolerance
[contact_matrix,R_matrix,~,slave_dof]=surface_to_surface(slave_surface,master_surface,nodes,1,true);
coor=nodes.';coor=coor(:);
gap=R_matrix*(contact_matrix*coor);
gap(1:3:end)=0;gap(2:3:end)=0;
coor(slave_dof)=coor(slave_dof)-R_matrix.'*gap;
nodes=reshape(coor,3,[]).';
err=norm(gap)/norm(coor(slave_dof));
flag=flag+1;
if flag>20
    waring('Unable to adjust to the given error')
end
end
%==========================================================================
Model.Mesh.nodes=nodes;
end