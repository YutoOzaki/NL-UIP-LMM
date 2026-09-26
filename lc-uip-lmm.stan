data {
  int<lower=2> p;
  int<lower=1> N;
  int<lower=1> M;
  int<lower=1> n[M];
  
  vector[p] m;
  int<lower=1> q;
  int<lower=1> r;
  
  matrix[2*N, p] X;
  vector[2*N] y;
  
  matrix[p, p] XX;
  real yy;
  vector[p] Xy;
  matrix[p, p] XX_p;
  real yy_p;
  vector[p] Xy_p;
  matrix[p, p] XX_g[M];
  real yy_g[M];
  matrix[M, p] Xy_g;
}

transformed data{
  real ln2pi = log(2*pi());
  real d_1 = 0.01;
  vector[2] d_v = [2*N, p^2]';
  real d_2 = max(d_v);
}

parameters {
  vector[p-1] v_other;
  real<lower=0> v_q;
  
  real<lower=0> sgm;
  real<lower=0> s_1;
  real<lower=0> s_2;
  real<lower=0, upper=1> u;
}

transformed parameters {
  real g = (1 - u)/u;
  
  // linear constraint
  vector[p] beta;
  {
  int k = 1;
  for(i in 1:p) {
    if(i != q){
      beta[i] = v_other[k];
      k += 1;
    }
  }
  
  beta[q] = beta[r] + v_q;
  }
}

model {
  // priors for variance components
  sgm ~ student_t(2, 0, 1000);
  s_1 ~ student_t(2, 0, 1000);
  s_2 ~ student_t(2, 0, 1000);
  u ~ beta(d_1, d_1*d_2);

  // Analytical and efficient likelihood evaluation
  real a = sgm^2;
  real b = s_1^2;
  real c = s_2^2;
  
  real lndetV = N*log(a) + (N-M)*log(a + 2*c);
    
  real w1 = 1/a;
  real w2 = c/(a*(a + 2*c));
    
  matrix[p, p] XVX = w1*XX - w2*XX_p;
  vector[p] XVy = w1*Xy - w2*Xy_p;
  real yVy = w1*yy - w2*yy_p;
    
  for(i in 1:M) {
    real denom = a + 2*c + 2*n[i]*b;
      
    lndetV += log(denom);
      
    real w3 = b/((a + 2*c)*denom);
    XVX -= w3*XX_g[i];
    XVy -= w3*to_vector(Xy_g[i, ]);
    yVy -= w3*yy_g[i];
  }
  
  // prior for beta
  matrix[p, p] Lmd = XVX/(2*N*g);
  beta ~ multi_normal_prec(m, Lmd);

  // evaluate likelihood
  real Q = quad_form_sym(XVX, beta) - 2*dot_product(beta, XVy) + yVy;
  real lpdf_likelihood_anal = -0.5*(2*N)*ln2pi - 0.5*lndetV - 0.5*Q;
  
  target += lpdf_likelihood_anal;
}

