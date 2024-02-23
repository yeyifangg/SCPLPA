
clear;
clc;

load('dataset5430.mat');    
 %LL:lncRNA functional similarities
 %DD:disease similarities
 
 DD=DD2;
 LL=MM;
  LD=miRNA_disease_Y;
  
  
  gama=0.9;
beta=0.9;
omega=0.6;
 
 zeta=1-omega;
  
    
 nl=size(LL,1);         %nm: No. of lnRNAs
 nd=size(DD,1);         %nd: No. of diseases

ddfs=getSimilarityDisease_3(LD',DD);% 
LL1=miRNASS( LD, DD );
llfs=getSimilarityRNA_1(LD',LL1);

LPALD= LPA(ddfs, llfs, LD,gama,beta);

  ldpl= PROJECTT(llfs,LPALD);
  ldpd= PROJECTT(ddfs,LPALD');
  

weight=ldpl*omega+ldpd'*zeta;
prevalue=weight;
prevalueldpl=ldpl;
prevalueldpd=ldpd;

y_train =LD;

nfolds =5430; nruns=1;
crossval_idx = crossvalind('Kfold',LD(:),nfolds);

for fold = 1:nfolds
    
    y_train = LD;
    test_idx  = find(crossval_idx==fold);
    y_train(test_idx) = 0;
                    
                    ddfsnew=getSimilarityDisease_3(y_train',DD);
                      LL1new=miRNASS(y_train, DD );
                   llfsnew=getSimilarityRNA_1(y_train',LL1new);
                 
                 
                 LPALDnew=LPA(ddfsnew, llfsnew, y_train,gama,beta);
                 
                  ldplnew= PROJECTT(llfsnew,LPALDnew);
                  ldpdnew= PROJECTT(ddfsnew,LPALDnew');
                  weightnew=ldplnew*omega+ldpdnew'*zeta;              
            
                prevalue(test_idx)= weightnew(test_idx);
                  prevalueldpl(test_idx)= ldplnew(test_idx);
                   prevalueldpd(test_idx)= ldpdnew(test_idx);
              
    end
 

% [X_1,Y_1,tpr,aupr_1] = perfcurve(LD(:), prevalue(:),1, 'xCrit', 'reca', 'yCrit', 'prec');
%  AUC=calculatewutuquzhong(LD,prevalue)
 AUCldpl=calculatewutuquzhong(LD,prevalueldpl)
 AUCldpd=calculatewutuquzhong(LD,prevalueldpd')

[X,Y,THRE,AUC,OPTROCPT,SUBY,SUBYNAMES] = perfcurve(LD(:), prevalue(:),1);
AUC=AUC
plot(X,Y)
save result prevalue prevalueldpl prevalueldpd



  
  
  
  
  
  
  
  
  
