%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Contact_pair=Update_contact_field(Contact_pair,Contact_condition,position,stepU,locking)
coor=position.';coor=coor(:);
for ci=1:1:numel(Contact_pair)
    contact_matrix=Contact_pair(ci).contact_matrix;
    R_matrix=Contact_pair(ci).R_matrix;
    slave_dof=Contact_pair(ci).slave_dof;
    clength=sqrt(Contact_pair(ci).area_vector(3:3:end));
    contact_status0=Contact_pair(ci).contact_status;
    mu=Contact_condition(ci).parameter.friction_coeff;
    kn=Contact_condition(ci).parameter.Emodulus;
    kt=Contact_condition(ci).parameter.Gmodulus;
    %======================================================================
    contact_gap=R_matrix*(contact_matrix*stepU);
    contact_gap(3:3:end)=R_matrix(3:3:end,:)*(contact_matrix*coor);
    gapT1=contact_gap(1:3:end);gapT2=contact_gap(2:3:end);
    directionFlag=gapT1>=0;
    %======================================================================
    direction=[gapT1,gapT2];
    direction=vecnorm(direction,2,2).\direction;
    direction(abs(direction)<1e-3)=0;
    loc=isnan(direction(:,1));direction(loc,:)=0;
    R1=permute(direction,[3,2,1]);R1(1,1,loc)=1;
    RR_matrix=zeros(3,3,size(R1,3));
    RR_matrix(1,1:2,:)=R1;
    RR_matrix(2,[2,1],:)=R1;
    RR_matrix(2,1,:)=-RR_matrix(2,1,:);
    RR_matrix(3,3,:)=1;
    IndexI=reshape(1:1:numel(slave_dof),3,1,[]);
    IndexI=repmat(IndexI,1,3,1);
    IndexJ=pagetranspose(IndexI);
    RR_matrix=sparse(IndexI(:),IndexJ(:),RR_matrix(:));
    contact_gap=RR_matrix*contact_gap;
    gapN=contact_gap(3:3:end);gapT=contact_gap(1:3:end);
    %======================================================================
    if strcmp(Contact_condition(ci).constraint,'tie')
        Contact_pair(ci).contact_gap=contact_gap;
        Contact_pair(ci).contact_stress=contact_stress;
        Contact_pair(ci).contact_status=true(numel(slave_dof),1);
        continue
    end
    %======================================================================
    if isempty(contact_status0)
        statusN=gapN>-1e-5*clength;
        gapN(statusN)=0;
    else
        statusN=contact_status0(3:3:end);
    end
    %======================================================================
    stressN=kn*gapN;
    stressN(~statusN)=0;
    stressT=kt*gapT;
    %======================================================================
    if ~locking
        statusN=gapN>=0;
        statusT=statusN;
        statusT(statusN)=max(mu*stressN(statusN),0)-stressT(statusN).*directionFlag(statusN)>=0;
    else
        statusT=contact_status0(1:3:end);
        statusN=contact_status0(3:3:end);
    end
    %======================================================================   
    meanL=mean(clength);
    remove_loc=clength/meanL<1e-2;
    statusT(remove_loc)=false;statusN(remove_loc)=false;
    %======================================================================
    stressT(~statusT)=mu*max(stressN(~statusT),0);
    stressN=max(kn*gapN,0);
    %--------------------------------------------------------------
    % figure(10);hold off;plot(contact_status0(1:3:end))
    % figure(10);hold on;plot(statusT)
    % figure(11);hold off;plot(Contact_pair(ci).contact_stress(1:3:end))
    % figure(11);hold on;plot(stressT)
    % figure(12);hold off;plot(contact_stress0(:,3))
    % figure(12);hold on;plot(stressN)
    %--------------------------------------------------------------
    contact_stress=[stressT, zeros(numel(stressT),1), stressN].';contact_stress=contact_stress(:);
    contact_status=[statusT, true(numel(stressT),1),  statusN].';contact_status=contact_status(:);
    %======================================================================
    Contact_pair(ci).status_change=~isequal(Contact_pair(ci).contact_status,contact_status);
    Contact_pair(ci).contact_gap=contact_gap;
    Contact_pair(ci).contact_stress=contact_stress;
    Contact_pair(ci).contact_status=contact_status;
    Contact_pair(ci).R_matrix=RR_matrix*R_matrix;
end
end