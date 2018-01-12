data {
    int n; // number of data points for m (samples)
    int m[n]; // number of mutations per sample in short time period (t)
    int M; // number of mutations since origin of CTVT, until MRCA of all samples (long time period, T)
    
    // prior parameters:
    real<lower=0> min_permitted_lower_bound;
    real<lower=0> max_permitted_upper_bound;
    real<lower=0> expon_rate;
}
transformed data {
    real PUPPY_AGE = 10.0 / 12; // age of puppy at time of sampling - 10 months
    real PARENT_AGE = 2.5; // age of parent dog at time of sampling - 2.5 years (excluding 6 months prior to sexual maturity)
}
parameters {
    real<lower=0,
         upper=(PARENT_AGE - PUPPY_AGE)> t_mrca_raw; // raw vector later used to construct length of time (years) before present to MRCA of all samples
    
    vector<lower=0,
           upper=10.0/12>[n] t; // length of time (years) before present to each sample's MRCA
    
    real<lower=0> rate; // rate parameter of the Poisson process (mutations per unit time [year])
    
    real<lower=min_permitted_lower_bound,
         upper=max_permitted_upper_bound> t_origin; // length of time (years) between origin of CTVT and MRCA of all samples
}
transformed parameters {
    real t_mrca = t_mrca_raw + PUPPY_AGE;
    vector<lower=0, upper=3>[n] i = t_mrca - t; // time intervals between MRCA of all samples, and MRCAs of each sample
}
model {
    t_mrca_raw ~ exponential(18) T[, PARENT_AGE - PUPPY_AGE];  // 95% of prior weight is on 2 months prior to birth (reflecting our belief that the parent was infected in the most recent heat cycle)
    for (j in 1:n) {
        t[j] ~ exponential(3.6) T[, PUPPY_AGE]; // 95% of prior weight is on the interval 0 - PUPPY_AGE
    }
    
    // prior for 'T' is Uniform ( [min_permitted, max_permitted) )
    t_origin ~ uniform(min_permitted_lower_bound, max_permitted_upper_bound);
    
    // prior on mutation rate is Exponential ( [0, inf) )
    // NB: Expectation(mutation rate) = 1/expon_rate
    rate ~ exponential(expon_rate);
    
    // likelihood of 'm' is Poisson, rate scaled by 't'
    for (j in 1:n) {
        m[j] ~ poisson(rate * i[j]);
    }
    
    // likelihood of 'M' is Poisson, rate scaled by 'T'
    M ~ poisson(rate * t_origin);
}
