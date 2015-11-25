%If data is not already loaded, uncomment this
% Authorize use of Quandl

%[data, headers] = Quandl.get('YAHOO/INDEX_GSPC','start_date','2010-11-23','end_date','2015-11-23','transformation','rdiff');

Close = data.Close;
obs = Close.data;

% the above pulls in SP500 data in using quandl
interval_len=floor(length(data)/10);

agg_returns=zeros(interval_len,1);
for i=1:interval_len
    j=10*i-9;
    agg_returns(i)=prod(1+obs(j:(j+9)));
end

agg_returns = agg_returns - 1;
% the above code generates a list of rates of increase for 10 day
% intervals. do this to avoid having really small observation values

% num_discrete_obs is the number of discrete observations
num_discrete_obs = 6;

timestep = (1:interval_len)';
hmm = [timestep HMM_obs(agg_returns,num_discrete_obs)];
agg_returns = [timestep agg_returns];
%hmm table records our DISCRETE observations for each timestep
%Use hmm for spectral method

%create empirical estimates of P1, P21, P3_x_1
%first randomly sample observation triples
%num_samples is the number of samples
num_samples = 100;
%also we truncate our observation array by the last two elements in order
%to sample the end triple
indices_to_sample = timestep(1:end-2);

rng('default');

sampled_indices=randsample(indices_to_sample,num_samples,true);
triples=zeros(num_samples,3);
for i=1:num_samples
    j=sampled_indices(i);
    triples(i,:)=hmm(j:j+2,2);
end

%now that we have our samples of observation triples, generate P1, P2,
%P3_x,1
P1 = zeros(1,num_discrete_obs);
%initialize P1
for i=1:num_discrete_obs
    P1(i) = sum(triples(:,1)==i)/num_samples;
end

P2 = zeros(num_discrete_obs);
%initialize P2
for i=1:num_discrete_obs
    for j=1:num_discrete_obs
        P2(i,j)=sum(ismember(triples(:,1:2),[i j],'rows'))/num_samples;
    end
end

%Create our P3_x_1 matrices
