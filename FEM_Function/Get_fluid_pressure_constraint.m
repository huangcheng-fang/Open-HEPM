%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function FPC=Get_fluid_pressure_constraint(Boundary,Node_set,nodes,P)
if numel(Boundary)==0;FPC=zeros(0,2);return;end
%--------------------------------------------------------------------------
container=Create_parallel_container(numel(Boundary));
for bi=1:1:numel(Boundary)
    setid=Boundary(bi).set;
    nid=Merge_cell(Node_set,setid);
    value=value_or_expr(Boundary(bi).value,Boundary(bi).expression,nodes(nid,:));
    dof=nid;
    container(bi).c1=[dof(:),value(:)];
end
FPC=Merge_container(container,1);
[~,loc]=unique(FPC(:,1),'last');
FPC=FPC(loc,:);
FPC(:,2)=FPC(:,2)-P(FPC(:,1));
FPC(:,1)=FPC(:,1)+numel(nodes);
end

function V=value_or_expr(value,expression,nodes)
if numel(expression)~=numel(value)
    error('The size of Expression should be same with Value');
end
x=nodes(:,1);y=nodes(:,2);
V=zeros(numel(x),numel(expression));
if size(nodes,2)==3;z=nodes(:,3);else;z=zeros(numel(x),1);end
for epi=1:1:numel(expression)
    if isempty(expression{epi})
        V(:,epi)=repmat(value(epi),size(nodes,1),1);
    else
        V(:,epi)=real(Simple_calculator(expression{epi},x,y,z));
    end
end
end