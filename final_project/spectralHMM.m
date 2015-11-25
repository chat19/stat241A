data = Quandl.get('YAHOO/INDEX_GSPC','start_date','2010-11-23','end_date','2015-11-23','transformation','rdiff');
Close = data.Close;
obs = Close.data;

% pulls in SP500 data in using quandl
interval_len=floor(length(data)/10);

agg_returns=zeros(interval_len,1);
for i=1:interval_len
    j=10*i-9;
    agg_returns(i)=prod(1+obs(j:(j+9)));
end

clear i j;

agg_returns = agg_returns - 1;
% the above code generates a list of rates of increase for 10 day
% intervals. do this to avoid having really small observation values

% num_discrete_obs is the number of discrete observations
num_discrete_obs = 6;

timestep = (1:interval_len)';
hmm = [timestep HMM_obs(agg_returns,num_discrete_obs)];
agg_returns = [timestep agg_returns];

%create empirical estimates of P1, P21, P3_x_1
%first randomly sample observation triples
%num_samples is the number of samples
num_samples = 100;
%also we truncate our observation array by the last two elements in order
%to sample the end triple
indices_to_sample = timestep(1:end-2);

sampled_indices=randsample(indices_to_sample,num_samples,true);
triples=zeros(num_samples,3);
for i=1:num_samples
    j=sampled_indices(i);
    triples(i,:)=hmm(j:j+2,2);
end
