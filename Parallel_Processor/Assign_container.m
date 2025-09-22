%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function container=Assign_container(container,C1,C2,C3,C4,C5,C6)
Tnum=size(C1,1);
Cnum=size(container,1);
Index=Split_parallel_data(Tnum,Cnum);
if ~isempty(C6)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
        container(ci).c2=C2(idx,:);
        container(ci).c3=C3(idx,:);
        container(ci).c4=C4(idx,:);
        container(ci).c5=C5(idx,:);
        container(ci).c6=C6(idx,:);
    end
    return
end
if ~isempty(C5)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
        container(ci).c2=C2(idx,:);
        container(ci).c3=C3(idx,:);
        container(ci).c4=C4(idx,:);
        container(ci).c5=C5(idx,:);
    end
    return
end
if ~isempty(C4)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
        container(ci).c2=C2(idx,:);
        container(ci).c3=C3(idx,:);
        container(ci).c4=C4(idx,:);
    end
    return
end
if ~isempty(C3)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
        container(ci).c2=C2(idx,:);
        container(ci).c3=C3(idx,:);
    end
    return
end
if ~isempty(C2)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
        container(ci).c2=C2(idx,:);
    end
    return
end
if ~isempty(C1)
    for ci=1:1:Cnum
        idx=Index(ci,1):Index(ci,2);
        container(ci).c1=C1(idx,:);
    end
    return
end
end