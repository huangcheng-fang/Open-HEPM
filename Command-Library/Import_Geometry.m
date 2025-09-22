%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Geometry=Import_Geometry(varargin)
%==========================Check input=====================================
fileName=[];scalePara=[];translatePara=[];rotatePara={};plotFlag=false;
for vi=1:2:numel(varargin)
    if isempty(varargin{vi});continue;end
    switch varargin{vi}
        case 'File'
            fileName=varargin{vi+1};
        case 'Scale'
            scalePara=varargin{vi+1};
        case 'Translate'
            translatePara=varargin{vi+1};
        case 'Rotate'
            rotatePara=varargin{vi+1};
        case 'Plot'
            plotFlag=varargin{vi+1};
        otherwise
            error(['Unknow input:',varargin{vi}])
    end
end
%========================main function=====================================
Geometry=importGeometry(fileName);
if ~isempty(scalePara)
    Geometry=scale(Geometry,scalePara);
end
if ~isempty(translatePara)
    Geometry=translate(Geometry,translatePara);
end
for i=1:1:size(rotatePara,1)
    Geometry=rotate(Geometry,rotatePara(i,4),[0,0,0],rotatePara(i,1:3));
end
if plotFlag
    pdegplot(Geometry,"VertexLabels","on","EdgeLabels","on","FaceLabels","on");
end
end