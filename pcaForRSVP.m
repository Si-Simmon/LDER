function [train_pca,test_pca,dataset_coef_inv] = pcaForRSVP(train,test,PCA_parameter)
%%

[mtrain,ntrain] = size(train);
[mtest,ntest] = size(test);
dataset = train;
[dataset_coef,dataset_score,dataset_latent] = pca(dataset); 
train_pca_1 = dataset_score(1:mtrain,:);
train_pca= train * inv(dataset_coef');
test_pca= test * inv(dataset_coef');
dataset_coef_inv = inv(dataset_coef');
train_pca = train_pca(:,1:PCA_parameter);
test_pca = test_pca(:,1:PCA_parameter);
