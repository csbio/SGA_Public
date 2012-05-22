%%
% MULTI_CLASS_LDA
%
% Authors: Chad Myers (cmyers@cs.umn.edu) 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%
% Last revision: 2010-07-19
%
%%

function [Xnorm,D,V,V2] =  multi_class_lda(X,classes,perc)

saveX=X;
a=sum(abs(X),1);
b=sum(abs(X),2);

X = X(b>0,a>0);
X(X==0)=normrnd(nanmean(X(:)),nanstd(X(:)),sum(sum(X==0)),1);

classes = classes(b>0);

class_labels = unique(classes);
dataMean = nanmean(X,1);
Sw=zeros(size(X,2));
Sb=zeros(size(X,2));


for i=1:length(class_labels)
  
  ind = find(classes == class_labels(i));
  
  if numel(ind) == 0  % empty class: ignore
  %if numel(ind) < 3  % orphan class: ignore
    continue;
  end
  
  classMean = nanmean(X(ind,:));
  
  Sw = Sw + cov(X(ind,:),1);

  Sb = Sb + numel(ind)*(classMean-dataMean)' ...
                      *(classMean-dataMean)  ;
   
end

eig_mat = pinv(Sw)*Sb;


if(5/perc < 20)
    [U,D,V] = svds(eig_mat,5/perc);

    a = diag(D)/max(diag(D));
    stopind = max(find(a >= perc));

    count=2;
    while(length(stopind) == 0)
        [U,D,V] = svds(eig_mat,(5/perc)*count);

        a = diag(D)/max(diag(D));
        stopind = max(find(a >= perc));    
        count = count+1;
    end
else
    [U,D,V] = svd(eig_mat);
    a = diag(D)/max(diag(D));
    stopind = max(find(a >= perc));    
end

N = V(:,1:stopind);

Xnorm = saveX;
a=sum(abs(Xnorm),1);
b=sum(abs(Xnorm),2);

%option 1: subtract projection of each vector
Xnorm(b>0,a>0)= Xnorm(b>0,a>0) - Xnorm(b>0,a>0)*N*N';
