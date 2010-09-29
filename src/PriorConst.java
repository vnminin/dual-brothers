interface PriorConst {

    /**
     * 0.5*ln(2*pi)
     */
    static final double normal_const = 0.91894; 
    
    /**
     * Hyperprior mean of mean of alpha [alpha_mean ~ Normal(alpha_mean_mean, alpha_mean_variance)]
     */
    static final double alpha_mean_mean = 0;

    
    /** Hyperprior variance of mean alpha [alpha_mean ~ Normal(alpha_mean_mean, alpha_mean_variance)]
     *
     */
    static final double alpha_mean_variance = 1000;

    
    /**
     * Hyperprior shape of precision of alpha, [precision ~ Gamma(shape, scale)]
     */
    static final double alpha_precision_shape = 0.01;

    
    /** Hyperprior scale of precision of alpha, [precision ~ Gamma(shape, scale)]
     *
     */
    static final double alpha_precision_scale = 0.01;

    /**
     * Hyperprior mean of mean of mu [mu_mean ~ Normal(mu_mean_mean, mu_mean_variance)]
     */
    static final double mu_mean_mean = 0;

    
    /** Hyperprior variance of mean mu [mu_mean ~ Normal(mu_mean_mean, mu_mean_variance)]
     *
     */
    static final double mu_mean_variance = 1000;

    
    /**
     * Hyperprior shape of precision of mu, [precision ~ Gamma(shape, scale)]
     */
    static final double mu_precision_shape = 0.01;

    
    /** Hyperprior scale of precision of alpha, [precision ~ Gamma(shape, scale)]
     *
     */
    static final double mu_precision_scale = 0.01;

    
}
