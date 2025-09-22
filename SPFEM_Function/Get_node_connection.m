%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Nconnection=Get_node_connection(elements)
elements(elements==0)=nan;
Num=max(elements,[],'all');
%==========================================================================
inelement=zeros(Num,200);
node_connect_element_num=zeros(Num,1);
for ei=1:1:size(elements,1)
    element=elements(ei,:);
    element(isnan(element))=[];
    inelement(element.'+Num*node_connect_element_num(element))=ei;
    node_connect_element_num(element)=node_connect_element_num(element)+1;
end
%================================Output====================================
Nconnection=[inelement(:,1:max(node_connect_element_num)),node_connect_element_num];
end