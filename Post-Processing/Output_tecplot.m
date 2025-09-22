%% Function Header Comment
% Developer: FANG Huangcheng @ PolyU
% Last updated: 2025-01-01
% Email: valy_f@bjtu.edu.cn;huangcheng.fang@polyu.edu.hk
% Website: https://www.researchgate.net/profile/Huangcheng-Fang
% Please do not remove this Header Comment under any circumstances, such as using or modifying this code, or convert this code to another programming language
function Output_tecplot(Model,Result,filename1,scale)
Nfield=Result.Nfield;
if isfield(Result,'Pfield')
    NEfield=Result.Pfield;
    NEfield.Stress=mean(NEfield.Stress,3);
    NEfield.Strain=mean(NEfield.Strain,3);
    NEfield.Pstrain=mean(NEfield.Pstrain,3);
    NEfield.Parameter=mean(NEfield.Parameter,3);
end
if isfield(Result,'Efield')
    NEfield=Mapping_field_E2N(Result.Efield,Model.Mesh);
end
%--------------------------------------------------------------------------
Mesh=Model.Mesh;
Etype=Mesh.etype;
Einpart=Mesh.einpart;
elements=Mesh.elements;
U=reshape(Result.Nfield.U,size(Mesh.nodes,2),[]).';
nodes=Mesh.nodes+(scale-1)*U;
Nstress=reshape(NEfield.Stress,6,[]).';
Nstrain=reshape(NEfield.Strain,6,[]).';
Npstrain=reshape(NEfield.Pstrain,6,[]).';
Npore_pressure=Result.Nfield.U*0;
%--------------------------------------------------------------------------
outputfile1=['Output-Files\',filename1];
fid=fopen(outputfile1,'w');
fprintf(fid,'%s\n',['Title=',filename1(1:end-4)]);
PN=unique(Einpart);
for ii=1:1:numel(PN)
    partid=PN(ii);
    EID=Einpart==partid;
    %----------------------------------------------------------------------
    Tet=elements(Etype==334&EID,1:4);
    nodeID=unique(Tet);
    nodeID(isnan(nodeID))=[];
    Index=zeros(max(nodeID),1);
    Index(nodeID,1)=1:size(nodeID,1);
    if ~isempty(Tet)
        fprintf(fid,'%s\n','Variables=x,y,z,Ux,Uy,Uz,Ex,Ey,Ez,Exy,Eyz,Exz,Sx,Sy,Sz,Sxy,Syz,Sxz,PEx,PEy,PEZ,PExy,PEyz,PExz,Pore_pressure,PartID');
        fprintf(fid,'%s\n',['Zone n=',num2str(size(nodeID,1)),',e=',num2str(size(Tet,1)),',f=fepoint,et=tetrahedron']);
        fpnode=[nodes(nodeID,1:3),U(nodeID,1:3),Nstrain(nodeID,1:6),Nstress(nodeID,1:6),Npstrain(nodeID,1:6),Npstrain(nodeID,1:1),repmat(partid,numel(nodeID),1)];
        fpnode(isnan(fpnode))=0;
        fprintf(fid,'%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g\n',fpnode.');
        fprintf(fid,'%8d%8d%8d%8d\n',Index(Tet(:,1:4).'));
    end
    %----------------------------------------------------------------------
    Hex=elements(Etype==338&EID,1:4);
    nodeID=unique(Hex);
    nodeID(isnan(nodeID))=[];
    Index=zeros(max(nodeID),1);
    Index(nodeID,1)=1:size(nodeID,1);
    if ~isempty(Hex)
        fprintf(fid,'%s\n','Variables=x,y,z,Ux,Uy,Uz,Ex,Ey,Ez,Exy,Eyz,Exz,Sx,Sy,Sz,Sxy,Syz,Sxz,PEx,PEy,PEZ,PExy,PEyz,PExz,Pore_pressure,PartID');
        fprintf(fid,'%s\n',['Zone n=',num2str(size(nodeID,1)),',e=',num2str(size(Hex,1)),',f=fepoint,et=brick']);
        fpnode=[nodes(nodeID,1:3),U(nodeID,1:3),Nstrain(nodeID,1:6),Nstress(nodeID,1:6),Npstrain(nodeID,1:6),Npore_pressure(nodeID),repmat(partid,numel(nodeID),1)];
        fprintf(fid,'%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g%16.8g\n',fpnode.');
        fprintf(fid,'%8d%8d%8d%8d%8d%8d%8d%8d\n',Index(Hex(:,1:8).'));
    end
end
fclose(fid);
%--------------------------------------------------------------------------
system(['tec360.exe ',outputfile1,'&']);
end