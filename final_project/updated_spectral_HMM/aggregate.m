function [ agg_returns ] = aggregate( data,num_days )
% Generates a list of rates of increase for num_days day
% intervals. Do this to avoid having really small observation values

interval_len=floor(length(data)/num_days);

agg_returns=zeros(interval_len,1);
for i=1:interval_len
    j=num_days*i - (num_days-1);
    agg_returns(i)=prod(1+data(j:(j+num_days-1)));
end

agg_returns = agg_returns - 1;

end

