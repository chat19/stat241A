% quandl_authorize
% 
% %If data is not already loaded, uncomment this
% %Authorize use of Quandl
% 
% SP500_data_collect;

% the above pulls in SP500 data in using quandl
% given our dataset and the number of states
function [ B, b_1, b_inf ] = HMM_calculate_params(data,num_states,num_discrete_obs,num_samples)

% num_discrete_obs is the number of discrete observations
% num_discrete_obs = 6;

% num_samples is the number of samples
% num_samples = 1000;

% this is our number of states
% num_states = 3;

num_days=10;
agg_returns = aggregate(data,num_days);
% generates a list of rates of increase for 10 day
% intervals. do this to avoid having really small observation values

timestep = (1:floor(length(data)/num_days))';
hmm = [timestep HMM_discretize(agg_returns,num_discrete_obs)];
% agg_returns = [timestep agg_returns];

% hmm table records our DISCRETE observations for each timestep
% Use hmm for spectral method

% create empirical estimates of P1, P21, P3_x_1
% first randomly sample observation triples
% also we truncate our observation array by the last two elements in order
% to sample the end triple
indices_to_sample = timestep(1:end-2);

rng('default');

sampled_indices=randsample(indices_to_sample,num_samples,true);
triples=zeros(num_samples,3);
for i=1:num_samples
    j=sampled_indices(i);
    triples(i,:)=hmm(j:j+2,2);
end

% now that we have our samples of observation triples, generate P1, P2_1,
% P3_x,1
P1 = zeros(1,num_discrete_obs)';
% initialize P1
for i=1:num_discrete_obs
    P1(i) = sum(triples(:,1)==i)/num_samples;
end

P2 = zeros(num_discrete_obs);
% initialize P2
for i=1:num_discrete_obs
    for j=1:num_discrete_obs
        P2(i,j)=sum(ismember(triples(:,1:2),[j i],'rows'))/num_samples;
    end
end

% Create our P3_x_1 matrices
% initialize P3 as a multidimensional array
% P3{i} corresponds to P3_x_1
H = zeros(num_discrete_obs,num_discrete_obs,num_discrete_obs);
P3 = cell(1,num_discrete_obs);
for x=1:num_discrete_obs
    for i=1:num_discrete_obs
        for j=1:num_discrete_obs
            H(i,j,x)=sum(ismember(triples,[j x i],'rows'))/num_samples;
        end
    end
    P3{x}=H(:,:,x);
end

clear H;

% now we compute model parameters
[A,~,~]=svd(P2);

% this is our number of states
% num_states=3;

U=A(:,1:num_states);

b_1 = U'*P1;
b_inf = pinv(P2' * U) * P1;
B = cell(1,num_discrete_obs);
for i=1:num_discrete_obs
    B{i} = U' * P3{i} * pinv(U'*P2);
end
