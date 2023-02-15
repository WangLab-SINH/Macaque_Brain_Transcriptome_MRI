%function PLS_bootstrap(response_var_file, predictor_var_file, output_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the PLS bootstrap function with the following arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% response_var_file ------ full path to the PLS_MRI_response_vars.csv file
%%%                           that is created by the NSPN_CorticalMyelination
%%%                           wrapper script
%%% predictor_var_file ----- full path to the PLS_gene_predictor_vars.csv file
%%%                           that is provided as raw data
%%% output_dir ------------- where to save the PLS_geneWeights and PLS_ROIscores
%%%                           files (for PLS1 and PLS2 separately)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Created by Petra Vertes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Running PLS')

%import response variables
%importdata(response_var_file);
importdata('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/vsd4_ncx_rmbatch_97regions_matchedROI173_ct.csv');

%unwrap and tidy MRI response variable names
ROIname=ans.textdata(:,1);
ResponseVarNames=ans.textdata(1,:);
ResponseVarNames=ans.textdata(1,2);
ResponseVarNames=ans.textdata(1,:);
ResponseVarNames(1)=[];
ROIname(1)=[];
%and store the response variables in matrix Y
MRIdata=ans.data;
clear ans

%import predictor variables
%indata=importdata(predictor_var_file);
indata=importdata('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/vsd4_ncx_rmbatch_97regions_meanExp.csv');
GENEdata=indata.data;
%GENEdata(1,:)=[];
genes=indata.textdata;
genes=genes(2:length(genes));
geneindex=1:length(genes);
clear indata

cd('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/')
%number of bootstrap iterations
bootnum=10000;

%DO PLS in 2 dimensions (with 2 components)
X=GENEdata';
Y=zscore(MRIdata);
Y=MRIdata;

dim=30;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
% Plot the percent of variance explained in the response variable as a function of the number of components
figure('visible','off');
plot(1:30,cumsum(100*PCTVAR(2,:)),'-bo');
xlabel('Number of PLS components');
ylabel('Percent Variance Explained in y');
saveas(gcf, 'PCTVAR.pdf');
% write PCTVAR into file
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PCTVAR.csv'),PCTVAR);


dim=5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);

%store regions IDs and weights in descending order of weight for both
%components
[R1,p1]=corr([XS(:,1),XS(:,2),XS(:,3),XS(:,4),XS(:,5)],MRIdata);

%align PLS components with desired direction%
% if R1(1,2)<0
%     stats.W(:,1)=-1*stats.W(:,1);
%     XS(:,1)=-1*XS(:,1);
% end
% if R1(2,4)<0
%     stats.W(:,2)=-1*stats.W(:,2);
%     XS(:,2)=-1*XS(:,2);
% end

[PLS1w,x1] = sort(stats.W(:,1),'descend');
PLS1ids=genes(x1);
geneindex1=geneindex(x1);
[PLS2w,x2] = sort(stats.W(:,2),'descend');
PLS2ids=genes(x2);
geneindex2=geneindex(x2);
[PLS3w,x3] = sort(stats.W(:,3),'descend');
PLS3ids=genes(x3);
geneindex3=geneindex(x3);
[PLS4w,x4] = sort(stats.W(:,4),'descend');
PLS4ids=genes(x4);
geneindex4=geneindex(x4);
[PLS5w,x5] = sort(stats.W(:,5),'descend');
PLS5ids=genes(x5);
geneindex5=geneindex(x5);

%print out results
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173.PLS1_ROIscores.csv'),XS(:,1));
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173.PLS2_ROIscores.csv'),XS(:,2));
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173.PLS3_ROIscores.csv'),XS(:,3));
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173.PLS4_ROIscores.csv'),XS(:,4));
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173.PLS5_ROIscores.csv'),XS(:,5));





%define variables for storing the (ordered) weights from all bootstrap runs
PLS1weights=[];
PLS2weights=[];
PLS3weights=[];
PLS4weights=[];
PLS5weights=[];

%start bootstrap
disp('  Bootstrapping - could take a while')
pctvar_m1=zeros(bootnum,5);
for i=1:bootnum
    disp(i)
    myresample = randsample(size(X,1),size(X,1),1);
    res(i,:)=myresample; %store resampling out of interest
    Xr=X(myresample,:); % define X for resampled regions
    Yr=Y(myresample,:); % define Y for resampled regions
    dim=5;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(Xr,Yr,dim); %perform PLS for resampled data
    pctvar_m1(i,:)=PCTVAR(2,:);
    temp=stats.W(:,1);%extract PLS1 weights
    newW=temp(x1); %order the newly obtained weights the same way as initial PLS
    if corr(PLS1w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS1weights=[PLS1weights,newW];%store (ordered) weights from this bootstrap run

    temp=stats.W(:,2);%extract PLS2 weights
    newW=temp(x2); %order the newly obtained weights the same way as initial PLS
    if corr(PLS2w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS2weights=[PLS2weights,newW]; %store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,3);%extract PLS3 weights
    newW=temp(x3); %order the newly obtained weights the same way as initial PLS
    if corr(PLS3w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS3weights=[PLS3weights,newW]; %store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,4);%extract PLS4 weights
    newW=temp(x4); %order the newly obtained weights the same way as initial PLS
    if corr(PLS4w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS4weights=[PLS4weights,newW]; %store (ordered) weights from this bootstrap run
    
    temp=stats.W(:,5);%extract PLS5 weights
    newW=temp(x5); %order the newly obtained weights the same way as initial PLS
    if corr(PLS5w,newW)<0 % the sign of PLS components is arbitrary - make sure this aligns between runs
        newW=-1*newW;
    end
    PLS5weights=[PLS5weights,newW]; %store (ordered) weights from this bootstrap run
    
end

csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_pctvar_m1.csv'),pctvar_m1);


%get standard deviation of weights from bootstrap runs
PLS1sw=std(PLS1weights');
PLS2sw=std(PLS2weights');
PLS3sw=std(PLS3weights');
PLS4sw=std(PLS4weights');
PLS5sw=std(PLS5weights');


%get bootstrap weights
temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';
temp3=PLS3w./PLS3sw';
temp4=PLS4w./PLS4sw';
temp5=PLS5w./PLS5sw';



%order bootstrap weights (Z) and names of regions
[Z1 ind1]=sort(temp1,'descend');
PLS1=PLS1ids(ind1);
geneindex1=geneindex1(ind1);

[Z2 ind2]=sort(temp2,'descend');
PLS2=PLS2ids(ind2);
geneindex2=geneindex2(ind2);

[Z3 ind3]=sort(temp3,'descend');
PLS3=PLS3ids(ind3);
geneindex3=geneindex3(ind3);

[Z4 ind4]=sort(temp4,'descend');
PLS4=PLS4ids(ind4);
geneindex4=geneindex4(ind4);

[Z5 ind5]=sort(temp5,'descend');
PLS5=PLS5ids(ind5);
geneindex5=geneindex5(ind5);


% get bootstrap weights distribution
PLS1weightsZ=[];
PLS2weightsZ=[];
PLS3weightsZ=[];
PLS4weightsZ=[];
PLS5weightsZ=[];

for i=1:bootnum
    % PLS1
    temp=PLS1weights(:,i);
    newW=temp;
    newW=newW./PLS1sw';
    PLS1weightsZ=[PLS1weightsZ,newW];
    
    % PLS2
    temp=PLS2weights(:,i);
    newW=temp;
    newW=newW./PLS2sw';
    PLS2weightsZ=[PLS2weightsZ,newW];
    
    % PLS3
    temp=PLS3weights(:,i);
    newW=temp;
    newW=newW./PLS3sw';
    PLS3weightsZ=[PLS3weightsZ,newW];
    
    % PLS4
    temp=PLS4weights(:,i);
    newW=temp;
    newW=newW./PLS4sw';
    PLS4weightsZ=[PLS4weightsZ,newW];
    
    % PLS5
    temp=PLS5weights(:,i);
    newW=temp;
    newW=newW./PLS5sw';
    PLS5weightsZ=[PLS5weightsZ,newW];
    
end    

% get p value
PLS1pvalue=[];
PLS2pvalue=[];
PLS3pvalue=[];
PLS4pvalue=[];
PLS5pvalue=[];


temp1=PLS1w./PLS1sw';
temp2=PLS2w./PLS2sw';
temp3=PLS3w./PLS3sw';
temp4=PLS4w./PLS4sw';
temp5=PLS5w./PLS5sw';


for i=1:length(genes)
    %PLS1
    a1=sum(PLS1weightsZ(i,:)>temp1(i));
    a2=sum(PLS1weightsZ(i,:)<temp1(i));
    a=min(a1,a2);
    a=a/1000;
    PLS1pvalue=[PLS1pvalue,a];
    
    %PLS2
    a1=sum(PLS2weightsZ(i,:)>temp2(i));
    a2=sum(PLS2weightsZ(i,:)<temp2(i));
    a=min(a1,a2);
    a=a/1000;
    PLS2pvalue=[PLS2pvalue,a];
    
    %PLS3
    a1=sum(PLS3weightsZ(i,:)>temp3(i));
    a2=sum(PLS3weightsZ(i,:)<temp3(i));
    a=min(a1,a2);
    a=a/1000;
    PLS3pvalue=[PLS3pvalue,a];
    
    %PLS4
    a1=sum(PLS4weightsZ(i,:)>temp4(i));
    a2=sum(PLS4weightsZ(i,:)<temp4(i));
    a=min(a1,a2);
    a=a/1000;
    PLS4pvalue=[PLS4pvalue,a];
    
    %PLS5
    a1=sum(PLS5weightsZ(i,:)>temp5(i));
    a2=sum(PLS5weightsZ(i,:)<temp5(i));
    a=min(a1,a2);
    a=a/1000;
    PLS5pvalue=[PLS5pvalue,a];
        
end    

PLS1pvalue=PLS1pvalue(ind1);
PLS2pvalue=PLS2pvalue(ind2);
PLS3pvalue=PLS3pvalue(ind3);
PLS4pvalue=PLS4pvalue(ind4);
PLS5pvalue=PLS5pvalue(ind5);



%print out results
fid1 = fopen(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS1_geneWeights.csv'),'w');
for i=1:length(genes)
  fprintf(fid1,'%s, %d, %f, %f\n', PLS1{i}, geneindex1(i), Z1(i), PLS1pvalue(i) );
end
fclose(fid1);

fid2 = fopen(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS2_geneWeights.csv'),'w');
for i=1:length(genes)
  fprintf(fid2,'%s, %d, %f, %f\n', PLS2{i}, geneindex2(i), Z2(i), PLS2pvalue(i) );
end
fclose(fid2);

fid3 = fopen(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS3_geneWeights.csv'),'w');
for i=1:length(genes)
  fprintf(fid3,'%s, %d, %f, %f\n', PLS3{i}, geneindex3(i), Z3(i), PLS3pvalue(i) );
end
fclose(fid3);

fid4 = fopen(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS4_geneWeights.csv'),'w');
for i=1:length(genes)
  fprintf(fid4,'%s, %d, %f, %f\n', PLS4{i}, geneindex4(i), Z4(i), PLS4pvalue(i)  );
end
fclose(fid4);

fid5 = fopen(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS5_geneWeights.csv'),'w');
for i=1:length(genes)
  fprintf(fid5,'%s, %d, %f, %f\n', PLS5{i}, geneindex5(i), Z5(i), PLS5pvalue(i)  );
end
fclose(fid5);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DO PLS in 5 dimensions (with 5 components)
X=GENEdata';
%Y=zscore(MRIdata);
Y=MRIdata;
dim=5;
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats]=plsregress(X,Y,dim);
temp=cumsum(100*PCTVAR(2,1:dim));
Rsquared = temp(dim);

%align PLS components with desired direction%
[R1,p1]=corr([XS(:,1),XS(:,2),XS(:,3),XS(:,4),XS(:,5)],MRIdata);

% if R1(1,2)<0
%     XS(:,1)=-1*XS(:,1);
% end
% if R1(2,4)<0
%     XS(:,2)=-1*XS(:,2);
% end

%calculate correlations of PLS components with MRI variables
% [R1,p1]=corr(XS(:,1),MRIdata);
% [R2,p2]=corr(XS(:,2),MRIdata);
% a=[R1',p1',R2',p2'];


%assess significance of PLS result
pctvar_m2=zeros(10000,5);
for j=1:bootnum
    disp(j)
    order=randperm(size(Y,1));
    Yp=Y(order,:);
    [XLr,YLr,XSr,YSr,BETAr,PCTVARr,MSEr,statsr]=plsregress(X,Yp,dim);
    pctvar_m2(j,:)=PCTVARr(2,:);
    temp=cumsum(100*PCTVARr(2,1:dim));
    Rsq(j) = temp(dim);
end
p=length(find(Rsq>=Rsquared))/j;
csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_pctvar_m2.csv'),pctvar_m2);

% plot histogram
figure('visible','off');
hist(Rsq,30)
hold on
plot(Rsquared,20,'.r','MarkerSize',15)
set(gca,'Fontsize',14)
xlabel('R squared','FontSize',14);
ylabel('Permuted runs','FontSize',14);
mytitle=sprintf('%s%f','p=',p);
title(mytitle)
saveas(gcf, 'p_hist.pdf');

save('vsd4_ncx_rmbatch_97regions_roi173.mat');
% 
% %save stats
% myStats=[PCTVAR; p, j];
% csvwrite(fullfile('F:/Lab_info/wanglab/My_project/crab-eating macaque/results/results_20201202/latest_pls_20210330/roi173_97NCXareas/pls_results/vsd4_ncx_rmbatch_97regions_roi173_PLS_stats.csv'),myStats);
