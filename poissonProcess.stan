data {
    int n; // number of data points for m (samples)
    int m[n]; // number of mutations per sample in short time period (t)
    int N; // number of data points for M
    int M[N]; // number of mutations since origin (long time period, T)
    
    // prior parameters:
    real<lower=0> a;
    real<lower=0> b;
    real<lower=0> min_permitted_lower_bound;
    real<lower=0> max_permitted_upper_bound;
    real<lower=0> expon_rate;
}
parameters {
    real<lower=0> rate; // rate parameter of the Poisson process (mutations per unit time [year])
    real<lower=0> t; // short time interval in which mutations are generated (years)
    real<lower=min_permitted_lower_bound,
         upper=max_permitted_upper_bound> T; // long time interval (years)
}
model {
    // prior for 't' is Beta value ( [0,1) )
    t ~ beta(a, b);
    
    // prior for 'T' is Uniform ( [min_permitted, max_permitted) )
    T ~ uniform(min_permitted_lower_bound, max_permitted_upper_bound);
    
    // prior on mutation rate is Exponential ( [0, inf) )
    // NB: Expectation(mutation rate) = 1/expon_rate
    rate ~ exponential(expon_rate);
    
    // likelihood of 'm' is Poisson, rate scaled by 't'
    m ~ poisson(rate * t);
    
    // likelihood of 'M' is Poisson, rate scaled by 'T'
    M ~ poisson(rate * T);
}
