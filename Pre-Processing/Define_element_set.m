%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Model=Define_element_set(Model,setID,varargin)
%==========================Check input=====================================
rangeX=[-inf,inf];rangeY=[-inf,inf];rangeZ=[-inf,inf];rangePart=[];R=eye(3);
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'RangeX'
            rangeX=varargin{vi+1};
        case 'RangeY'
            rangeY=varargin{vi+1};
        case 'RangeZ'
            rangeZ=varargin{vi+1};
        case 'RangePart'
            rangePart=varargin{vi+1};
        case 'Rotation'
            theta=varargin{vi+1}/180*pi;
            thetax=theta(1);thetay=theta(2);thetaz=theta(3);
            a=cos(thetax);b=sin(thetax);
            Rx=[1,0,0;0,a,b;0,-b,a];
            a=cos(thetay);b=sin(thetay);
            Ry=[a,0,-b;0,1,0;b,0,a];
            a=cos(thetaz);b=sin(thetaz);
            Rz=[a,b,0;-b,a,0;0,0,1];
            R=Rx*Ry*Rz;
        otherwise
            disp(['Unknow input type is ignored:',varargin{vi}])
    end
end
%========================main function=====================================
Mesh=Model.Mesh;
%--------------------------------------------------------------------------
if ~isempty(rangePart)
    loc=find(ismember(Mesh.einpart,rangePart));
else
    loc=(1:1:size(Mesh.einpart,1)).';
end
%--------------------------------------------------------------------------
coor=Model.Mesh.ecenter(loc,:);
coor=coor*R(1:size(coor,2),1:size(coor,2));
lim1=[rangeX(1),rangeY(1),rangeZ(1)];
lim2=[rangeX(2),rangeY(2),rangeZ(2)];
coor1=coor-lim1(:,1:size(coor,2));
coor2=coor-lim2(:,1:size(coor,2));
flag=coor1.*coor2;
Model.Set.element_set{setID,1}=loc(flag(:,1)<=0&flag(:,2)<=0&flag(:,end)<=0);
%--------------------------------------------------------------------------
end