%%
% MULTI_CLASS_LDA
%
% Authors: Chad Myers (cmyers@cs.umn.edu) 
%          Anastasia Baryshnikova (a.baryshnikova@utoronto.ca)
%          Benjamin VanderSluis (bvander@cs.umn.edu)
%
% Last revision: 2012-06-06
%
%%

function [batch_effect] =  multi_class_lda(X,classes,sv_thresh)

% X is plates by array positions (i.e. N x 1536)
batch_effect=X;
batch_effect(isnan(batch_effect))=0;

% ignore any plate or array position with no data
% and replace remaining nans with random data
valid_arrays = ~all(isnan(X),1);
valid_plates = ~all(isnan(X),2);
X = X(valid_plates, valid_arrays); 
X(isnan(X))=normrnd(nanmean(X(:)),nanstd(X(:)),sum(isnan(X(:))),1);

classes = classes(valid_plates);
class_labels = unique(classes);
dataMean = mean(X,1);
Sw=zeros(size(X,2));
Sb=zeros(size(X,2));


for i=1:length(class_labels)
  
  ind = find(classes == class_labels(i));
  
  if numel(ind) < 3  % orphan class: ignore
    continue;
  end
  
  classMean = mean(X(ind,:));
  
  Sw = Sw + cov(X(ind,:),1);
  Sb = Sb + numel(ind)*(classMean-dataMean)' ...
                      *(classMean-dataMean)  ;
end

eig_mat = pinv(Sw)*Sb;

% If the first pass doesn't give us enough singular
% values, go back and get more
stopind = 0;
scaled_E = [];
count = 1;
while(stopind == length(scaled_E))
    [U,E,V] = svds(eig_mat,ceil((5/sv_thresh)*count));

    scaled_E = diag(E)/max(diag(E));
    stopind = max(find(scaled_E >= sv_thresh));    
    count = count+1;
end

% project onto batch subspace 
N = V(:,1:stopind);
batch_effect(valid_plates, valid_arrays)= ...
  batch_effect(valid_plates, valid_arrays)*N*N';
