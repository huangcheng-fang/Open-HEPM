%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Plot_contact(ContactPair,nodes,P,type)
%==============================main function===============================
for ci=1:1:numel(ContactPair)
    slave_dof=ContactPair(ci).slave_dof;
    nid=slave_dof(3:3:end)/3;
    x=nodes(nid,1);y=nodes(nid,2);z=nodes(nid,3);
    switch type
        case 'contact_status'
            value=ContactPair(ci).contact_status;
            value=[num2str(value(1:3:end)),repmat(',',numel(x),1),num2str(value(2:3:end)),repmat(',',numel(x),1),num2str(value(3:3:end))];
        case 'contact_stress'
            value=ContactPair(ci).contact_stress;
            value=[num2str(value(3:3:end),3)];
        case 'contact_gap'
            value=ContactPair(ci).contact_gap(3:3:end);
            value=num2str(value);
        case 'contact_area'
            value=ContactPair(ci).area_vector;
            value=[num2str(value(1:3:end)),repmat(',',numel(x),1),num2str(value(2:2:end))];
        case 'node_number'
            value=ContactPair(ci).slave_dof(2:2:end)/2;
            value=num2str(value);
        case 'normal_vector'
            [~,~,normal_vector]=find(ContactPair(ci).R_matrix(3:3:end,:));
            normal_vector=reshape(normal_vector,2,[]).'.*abs(ContactPair(ci).area_vector(3:3:end));
            plot3(x+normal_vector(:,1),y+normal_vector(:,2),z+normal_vector(:,3),'.r')
            x=zeros(0,1);y=zeros(0,1);value=char(0);
        otherwise
            error('Unknown type')
    end
    text(x,y,z,value)
end
end