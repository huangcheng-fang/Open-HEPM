%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function value=Simple_calculator(expr,x,y,z)
value=eval(vectorize(expr));
return
expr=convert(expr);
value=complex(zeros(numel(x),1));
for ii=1:1:numel(x)
    temp=expr;
    while contains(temp,'(')
    loc1=find(temp=='(',1,'last');loc1=loc1(1);
    loc2=loc1+find(temp(loc1+1:end)==')',1,'first');loc2=loc2(1);
    pure_expr=temp(loc1+1:loc2-1);
    v=pure_math(pure_expr,x(ii),y(ii),z(ii));
    temp=[temp(1:loc1-1),sprintf('%0.16f',v),temp(loc2+1:end)];
    end
    if strcmp(temp(1:2),'--');temp(1:2)=[];end
    value(ii)=pure_math(temp,x(ii),y(ii),z(ii));
end
value=reshape(value,size(x,1),size(x,2));
end

function value=pure_math(expr,x,y,z)
num=numel(expr);
data=complex(zeros(1,num*0+50));
math=zeros(1,num*0+50);
%-------------------Separating data from operators-------------------------
flag=1;ni=1;
for ii=1:1:num
    if contains('+-*/^SCTP',expr(ii))
        switch expr(flag:ii-1)
            case 'x'
                data(1,ni)=x;
            case 'y'
                data(1,ni)=y;
            case 'z'
                data(1,ni)=z;
            otherwise
                temp=str2double(expr(flag:ii-1));
                data(1,ni)=temp;
                if isnan(temp)
                    if expr(ii)=='-'
                        continue
                    else
                        error(['Incorrect expression:',expr(flag:ii-1)])
                    end
                end
        end
        switch expr(ii)
            case '+'
                math(1,ni)=1;
            case '-'
                math(1,ni)=-1;
            case '*'
                math(1,ni)=2;
            case '/'
                math(1,ni)=-2;
            case '^'
                math(1,ni)=3;
            case 'S'
                math(1,ni)=4;
            case 'C'
                math(1,ni)=5;
            case 'T'
                math(1,ni)=6;
            case 'P'
                math(1,ni)=7;
            otherwise
                error(['Unsupported operator:',expr(ii)])
        end
        flag=ii+1;ni=ni+1;
    end
end
switch expr(flag:end)
    case 'x'
        data(1,ni)=x;
    case 'y'
        data(1,ni)=y;
    case 'z'
        data(1,ni)=z;
    otherwise
        temp=str2double(expr(flag:end));
        data(1,ni)=temp;
        if isnan(temp)
            error(['Incorrect expression:',expr(flag:end)])
        end
end
%---------------------------calculate expr---------------------------------
math(ni:end)=[];
if isempty(math)
    value=data(1,1);
    return
end
value=zeros(numel(x),1);
while ~isempty(math)
    [~,loc]=max(abs(math));
    switch math(loc)
        case 1
            value=data(1,loc)+data(1,loc+1);
        case -1
            value=data(1,loc)-data(1,loc+1);
        case 2
            value=data(1,loc)*data(1,loc+1);
        case -2
            value=data(1,loc)/data(1,loc+1);
        case 3
            value=data(1,loc)^data(1,loc+1);
        case 4
            value=sin(data(1,loc+1));
        case 5
            value=cos(data(1,loc+1));
        case 6
            value=tan(data(1,loc+1));
        case 7
            value=exp(data(1,loc+1));
        otherwise
            error('wrong')
    end
    math(loc)=[];
    data(loc)=[];
    data(loc)=value;
end
end


function expr=convert(expr)
expr=strrep(expr,'sin','0S');
expr=strrep(expr,'cos','0C');
expr=strrep(expr,'tan','0T');
expr=strrep(expr,'exp','0P');
expr=strrep(expr,' ','');
expr=strrep(expr,'pi',sprintf('%0.16f',pi));
end