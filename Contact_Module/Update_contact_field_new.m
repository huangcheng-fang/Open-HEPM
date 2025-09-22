%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Update_contact_field_new(Contact_pair,Contact_condition,position,stepU,locking)
coor=position.';coor=coor(:);
for ci=1:1:numel(Contact_pair)
    contact_matrix=Contact_pair(ci).contact_matrix;
    R_matrix=Contact_pair(ci).R_matrix;
    slave_dof=Contact_pair(ci).slave_dof;
    clength=sqrt(Contact_pair(ci).area_vector(3:3:end));
    contact_status0=Contact_pair(ci).contact_status;
    contact_stress0=reshape(Contact_pair(ci).contact_stress,3,[]).';
    mu=Contact_condition(ci).parameter.friction_coeff;
    kn=Contact_condition(ci).parameter.Emodulus;
    kt=Contact_condition(ci).parameter.Gmodulus;
    %======================================================================
    contact_gap=R_matrix*(contact_matrix*stepU);
    contact_gap(3:3:end)=R_matrix(3:3:end,:)*(contact_matrix*coor);
    gapN=contact_gap(3:3:end);gapT1=contact_gap(1:3:end);gapT2=contact_gap(2:3:end);
%     gapN(abs(gapN)<1E-4)=0;
    %======================================================================
    if strcmp(Contact_condition(ci).constraint,'tie')
        contact_stress=[kt;kt;kn].*reshape(R_matrix*(contact_matrix*stepU),3,[]);
        Contact_pair(ci).contact_gap=contact_gap;
        Contact_pair(ci).contact_stress=contact_stress(:);
        Contact_pair(ci).contact_status=true(numel(slave_dof),1);
        continue
    end
    %======================================================================
    if strcmp(Contact_condition(ci).constraint,'frictionless')
        stressN=max(kn*gapN,0);
        statusN=gapN>=0;
        contact_stress=[zeros(size(stressN,1),2), stressN].';contact_stress=contact_stress(:);
        contact_status=[false(size(stressN,1),2), statusN].';contact_status=contact_status(:);
        Contact_pair(ci).contact_gap=contact_gap;
        Contact_pair(ci).contact_stress=contact_stress;
        Contact_pair(ci).contact_status=contact_status;
        continue
    end
    %======================================================================
    if norm([gapT1,gapT2])==0
        statusN=gapN>-1e-5*clength;
        gapN(statusN)=0;
        contact_stress0=zeros(numel(statusN),3);
        direction=zeros(numel(statusN),2);
    else
        statusN=contact_status0(3:3:end);
        statusT=contact_status0(1:3:end);
        direction=[gapT1,gapT2];
        direction(~statusT,:)=contact_stress0(~statusT,1:2);
        temp=vecnorm(direction,2,2);direction(temp==0,:)=[gapT1(temp==0,:),gapT2(temp==0,:)];
        direction=vecnorm(direction,2,2).\direction;
    end
    directionFlag=sum([gapT1,gapT2].*contact_stress0(:,1:2),2)>=0;
    %======================================================================
    stressN=kn*gapN;
    stressN(~statusN)=0;
    stressT1=kt*gapT1;
    stressT2=kt*gapT2;
    stressT=sqrt(stressT1.^2+stressT2.^2);
    %======================================================================
    if ~locking
        statusN=gapN>=0;
        statusT=statusN;
        statusT(statusN)=max(mu*stressN(statusN),0)-stressT(statusN).*directionFlag(statusN)>=0;
    else
        statusT=contact_status0(1:3:end);
        statusN=contact_status0(3:3:end);
    end
    if mu==0
        statusT(:)=false;
    end
    %======================================================================   
    meanL=mean(clength);
    remove_loc=clength/meanL<1e-2;
    statusT(remove_loc)=false;statusN(remove_loc)=false;
    %======================================================================
    maxfriction=mu*max(stressN(~statusT),0);
    stressT1(~statusT)=maxfriction.*direction(~statusT,1);
    stressT2(~statusT)=maxfriction.*direction(~statusT,2);
    stressN=max(kn*gapN,0);
    %--------------------------------------------------------------
    % figure(10);hold off;plot(contact_status0(1:3:end))
    % figure(10);hold on;plot(statusT)
    % figure(11);hold off;plot(contact_stress0(:,1))
    % figure(11);hold on;plot(stressT1)
    % figure(12);hold off;plot(slave_dof(3:3:end)/3,contact_stress0(:,3))
    % figure(12);hold on;plot(slave_dof(3:3:end)/3,stressN)
    %--------------------------------------------------------------
    contact_stress=[stressT1, stressT2, stressN].';contact_stress=contact_stress(:);
    contact_status=[statusT , statusT , statusN].';contact_status=contact_status(:);
    %======================================================================
    Contact_pair(ci).status_change=~isequal(Contact_pair(ci).contact_status,contact_status);
    Contact_pair(ci).contact_gap=contact_gap;
    Contact_pair(ci).contact_stress=contact_stress;
    Contact_pair(ci).contact_status=contact_status;
    % Contact_pair(ci).NCstress=zeros(size(position,1),1);
    % Contact_pair(ci).NCstatus=false(size(position,1),1);
    % Contact_pair(ci).NCstress(slave_dof)=R_matrix.'*contact_stress.*Contact_pair(ci).area_vector;
    % Contact_pair(ci).NCstatus(slave_dof)=Contact_pair(ci).contact_status;
end
end