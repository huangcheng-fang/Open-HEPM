clc;clear all;close all;
Key=Keywords();
Model=Creat_Model();
Geometry=Import_Geometry('File',"Cavity expansion2.stp",'Scale',1000);
Model=Generate_Mesh(Model,Geometry,1,"Hmin",3,"Hface",{[0.2,0.2,0.2,0.5;]},"Hgrad",1.2,'Order',1);
Plot_model(Model,2);

u=0;
RF=0;
Ur=0;dUr=0.05;

tempU=0;
Model.Set.node_set{4}=[];
addNodes=[];
Geometry=Generate_geometry(Model.Mesh,1,addNodes);
Model=Generate_Mesh(Model,Geometry,1,"Hmin",3,"Hmax",3,"Hface",{[0.2,0.2,0.2,2*(0+2)*pi/4/20;]},"Hgrad",1.1,'Order',1);


Model=Define_surface_set(Model,1,'RangeR',[0,0,0,1.8,2.2+Ur],'RangeX',[0.001,inf],'RangeY',[0.001,inf],'RangeZ',[0.001,inf]);
% Model=Define_surface_set(Model,2,'RangeR',[0,0,0,4.9,5.1],'RangeX',[0.001,inf],'RangeY',[0.001,inf],'RangeZ',[0.001,inf]);

Model=Define_node_set(Model,1,'RangeX',[-0.01,0.01]);
Model=Define_node_set(Model,2,'RangeY',[-0.01,0.01]);
Model=Define_node_set(Model,3,'RangeZ',[-0.01,0.01]);
Model.Set.node_set{4}=unique(Model.Mesh.surfaces(Model.Set.surface_set{1},:));

bnodes=Model.Mesh.nodes(Model.Set.node_set{4},:);
br=vecnorm(bnodes,2,2);
normal_vector=br.\bnodes;
bnodes=bnodes+(Ur+2-br).*normal_vector;
Model.Mesh.nodes(Model.Set.node_set{4},:)=bnodes;
Model.Mesh.positions(Model.Set.node_set{4},:)=bnodes;


Model=Define_element_set(Model,1,'RangePart',1);
Model=Define_particle_set(Model,1,'RangePart',1);

elasticity=struct('type',Key.linear,'parameter',[500e3,0.3]);
Model=Define_material_property(Model,1,'Eset',[1],'Pset',[1],'Elasticity',elasticity);


Model=Define_boundary_condition(Model,1,'Type',Keywords().displacement,'Nset',[1],'Expression',{0},'Direction',[1]);
Model=Define_boundary_condition(Model,2,'Type',Keywords().displacement,'Nset',[2],'Expression',{0},'Direction',[2]);
Model=Define_boundary_condition(Model,3,'Type',Keywords().displacement,'Nset',[3],'Expression',{0},'Direction',[3]);


Tolerance=struct('Utol',1e-3,'Ftol',1e-3,'PCGtol',1e-8,'MaxIterNum',100);

Result=Initialize_FEM_Solver(Model,'GeometricNonlinear',1,'Tolerance',Tolerance);
for i=1:1:80
incrRatio=1;dincrRatio=1;
% Model=Define_load_condition(Model,4,'Type',Key.pressure,'Sset',[1],'Expression',{incrRatio*Ur},'Direction',[1]);
Model=Define_boundary_condition(Model,4,'Type',Key.displacement,'Nset',[4],'Expression',{[num2str(Ur+dUr),'*x/sqrt(x^2+y^2+z^2)'],[num2str(Ur+dUr),'*y/sqrt(x^2+y^2+z^2)'],[num2str(Ur+dUr),'*z/sqrt(x^2+y^2+z^2)']},'Direction',[1,2,3]);
[Result,Model]=Submit_To_FEM_Solver_Implicit(Result,Model,'TimeStep',100000);
Plot_result(Result,Model.Mesh,23,'Sx',1,1);drawnow
Ur=Ur+dUr;
end