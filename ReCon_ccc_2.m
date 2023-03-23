
function results = ReCon_ccc_2(data,data_,cfg)
disp('Start ..');
results.cfg = cfg;
data = data(:)';
srate = cfg.srate;
epoch_twd = cfg.epoch_twd;
base_twd = cfg.base_twd;
resync_twd = cfg.resync_twd;
% lowpass_freq = cfg.lowpass_freq;
latencies = cfg.latencies;
latencies_ns = cfg.latencies_ns;
if ~isfield(cfg,'fig_visible') cfg.fig_visible = 'on';end
if ~isfield(cfg,'glb') cfg.glb = 0;end
%handling error
for hand_e = 1:1
if base_twd(1)< epoch_twd(1) || base_twd(2) > epoch_twd(2)
    ME = MException('baseline:exceed', ...
        'baseline window should not exceed epoch time window');
    throw(ME)
end
if resync_twd(2) - resync_twd(1) < 50
    ME = MException('template:window', ...
        'template time widnow should not be shorter than 50 ms');
    throw(ME)
end
if epoch_twd(1) > -100
    ME = MException('epoch:window', ...
        'left end of epoch time window should be <= 100 ms ');
    throw(ME)
end
end

sample_twd = fix(epoch_twd(1)*srate/1000):fix(epoch_twd(2)*srate/1000);

regressor = zeros(length(sample_twd),length(data));
for j = 1:length(sample_twd)
    temp = latencies + sample_twd(j);
    temp(temp<1|temp>length(data)) = [];
    regressor(j,temp) = 1;
end
ERP = regressor*data'/length(latencies);

%baseline
sample_base_twd = (round((base_twd(1)-epoch_twd(1))*srate/1000)+1):(round((base_twd(2)-epoch_twd(1))*srate/1000)+1);

ERP = ERP - mean(ERP(sample_base_twd));


single_trial = [];
%get single trials
for j = 1:length(latencies)
    sample2take = latencies(j) + sample_twd;
    sample2take(sample2take<1) = 1;
    sample2take(sample2take>length(data)) = length(data);
    single_trial(:,j) = data(sample2take);
    single_trial(:,j) = single_trial(:,j) - ...
        mean(single_trial(sample_base_twd,j));
end

sample_resync_twd = (round((resync_twd(1)-epoch_twd(1))*srate/1000)+1):(round((resync_twd(2)-epoch_twd(1))*srate/1000));
template = ERP(sample_resync_twd);


for i = 1:size(data_,1)
regressor_ = zeros(length(sample_twd),size(data_,2));
for j = 1:length(sample_twd)
    temp_ = latencies_ns + sample_twd(j);
    temp_(temp_<1|temp_>length(data_(i,:))) = [];
    regressor_(j,temp_) = 1;
end
template_1(i,:) = regressor_*data_(i,:)'/length(latencies_ns);
ERP_(i,:) = template_1(i,:) - mean(template_1(i,sample_base_twd));
end
template_2 = ERP_(:,sample_resync_twd);

for j = 1:length(latencies)
    sample2take = latencies(j) + sample_resync_twd + round(epoch_twd(1)*srate/1000);
    sample2take(sample2take<1) = 1;
    sample2take(sample2take>size(data_,2)) = size(data_,2);
    data2match = data_(:,sample2take);
    est_latency(j) = Estimate_lantency_ccc(template_2,data_(:,latencies(j)+51:latencies(j)+epoch_twd(end)*srate/1000-20));
end
est_latency = round(est_latency - resync_twd(1)*srate/1000 +50);

regressor2 = zeros(length(sample_twd),length(data));
latencies2 = latencies + est_latency;
for j = 1:length(sample_twd)
    latencies_temp = latencies2 + sample_twd(j);
    latencies_temp(latencies_temp<1|latencies_temp>length(data)) = [];
    regressor2(j,latencies_temp) = 1;
end

ERPs = [ERP,zeros(length(ERP),1)];
if std(est_latency) ~= 0
invR = inv([regressor;regressor2]*[regressor;regressor2]')*[regressor;regressor2];
ERPs = invR*data';
ERPs = reshape(ERPs,length(ERP),2);

%re-baseline 
for j = 1:size(ERPs,2)
    ERPs(:,j) = ERPs(:,j) - mean(ERPs(sample_base_twd,j));
end
end
    
ERP_re = sum(ERPs,2);

ST = ERPs(:,1)'*regressor;
ST2 = ERPs(:,2)'*regressor2;
ST_re = ERP_re'*regressor;
data_re = data - ST - ST2 + ST_re;
single_trial_re = [];
%get single trials
for j = 1:length(latencies)
    sample2take = latencies(j) + sample_twd;
    sample2take(sample2take<1) = 1;
    sample2take(sample2take>length(data)) = length(data);
    single_trial_re(:,j) = data_re(sample2take);
    single_trial_re(:,j) = single_trial_re(:,j) - ...
    mean(single_trial_re(sample_base_twd,j));
end

results.original_data = data;
results.resync_data = data_re;
results.est_latency = est_latency;
results.decomp_ERPs = ERPs;
results.original_ERP = ERP;
results.resync_ERP = ERP_re;

    


