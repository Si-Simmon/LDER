function [class_1,POSTERIOR_test_1,class_2,POSTERIOR_test_2] = sthcp(train_t,train_nt,trainData,testData_2,testData_3,ch_num,train_label,trainNum,testNum)

PCA_parameter = 7;


[Q,Q_rest] = CG_CSP_train(train_t,train_nt,ch_num/2);
for i = 1:trainNum
    trainData_csp(:,:,i) = Q * trainData(:,:,i);
end
for i = 1:testNum
    testData_csp(:,:,i) = Q * testData_2(:,:,i);
end
for i = 1:testNum
    testData_csp_2(:,:,i) = Q * testData_3(:,:,i);
end

for i=1:ch_num
    for j=1:trainNum
        eval(['data_singelch_train_',num2str(i),'(j,:)=trainData_csp(i,:,j);']);
        eval(['data_singelch.data_singelch_train_',num2str(i),' = data_singelch_train_',num2str(i),';']);
    end
end

for i=1:ch_num
    for j=1:testNum
        eval(['data_singelch_test_',num2str(i),'(j,:)=testData_csp(i,:,j);']);
    end
end

for i=1:ch_num
    for j=1:testNum
        eval(['data_singelch_test3_',num2str(i),'(j,:)=testData_csp_2(i,:,j);']);
    end
end



for i=1:ch_num
    eval(['[train_pca_',num2str(i),',test_pca_',num2str(i),',dataset_coef_inv_',num2str(i),']=pcaForRSVP(data_singelch_train_',num2str(i),',data_singelch_test_',num2str(i),',PCA_parameter);']);
    eval(['dataset_coef.dataset_coef_inv_',num2str(i),' = dataset_coef_inv_',num2str(i),';']);
end

for i=1:ch_num
    eval(['[train_pca_',num2str(i),',test3_pca_',num2str(i),',dataset_coef_inv_',num2str(i),']=pcaForRSVP(data_singelch_train_',num2str(i),',data_singelch_test3_',num2str(i),',PCA_parameter);']);
    eval(['dataset_coef.dataset_coef_inv_',num2str(i),' = dataset_coef_inv_',num2str(i),';']);
end




for i=0:ch_num-1
    eval(['train_pca(:,1+PCA_parameter*i:PCA_parameter+PCA_parameter*(i))=train_pca_',num2str(i+1),';']);
end

for i=0:ch_num-1
    eval(['test_pca(:,1+PCA_parameter*i:PCA_parameter+PCA_parameter*(i))=test_pca_',num2str(i+1),';']);
end

for i=0:ch_num-1
    eval(['test3_pca(:,1+PCA_parameter*i:PCA_parameter+PCA_parameter*(i))=test3_pca_',num2str(i+1),';']);
end


[class_1,err_1,POSTERIOR_test_1,logp_1,coeff_1] =classify_ccc(test_pca,train_pca,train_label','diagLinear');
[class_2,err_2,POSTERIOR_test_2,logp_2,coeff_2] =classify_ccc(test3_pca,train_pca,train_label','diagLinear');

end

