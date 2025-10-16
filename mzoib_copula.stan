functions {
  real clamp01(real u) {
    return fmin(1 - 1e-12, fmax(1e-12, u));
  }

  // ---- Gaussian copula pieces ----
  // z = Phi^{-1}(u)
  real z_of_u(real u_in) {
    return inv_Phi(clamp01(u_in));
  }

  // c(u,v; rho) for continuous–continuous
  real gauss_c(real u_in, real v_in, real rho) {
    real u = clamp01(u_in);
    real v = clamp01(v_in);
    real zu = inv_Phi(u);
    real zv = inv_Phi(v);
    real omr2 = 1 - rho * rho;

    // c = 1/sqrt(1-rho^2) * exp( - (zu^2 - 2rho zu zv + zv^2)/(2(1-rho^2)) + (zu^2+zv^2)/2 )
    real quad = (zu*zu - 2*rho*zu*zv + zv*zv) / (2*omr2) - 0.5*(zu*zu + zv*zv); // ← 修正这里

    return inv_sqrt(omr2) * exp(-quad);
  }

  // ∂C/∂u(u,v; rho) = Phi( (zv - rho*zu)/sqrt(1-rho^2) )
  real gauss_dC_du(real u_in, real v_in, real rho) {
    real u = clamp01(u_in);
    real v = clamp01(v_in);
    real zu = inv_Phi(u);
    real zv = inv_Phi(v);
    real s = (zv - rho * zu) / sqrt(1 - rho * rho);
    return Phi(s);
  }

  // ∂C/∂v(u,v; rho) = Phi( (zu - rho*zv)/sqrt(1-rho^2) )
  real gauss_dC_dv(real u_in, real v_in, real rho) {
    real u = clamp01(u_in);
    real v = clamp01(v_in);
    real zu = inv_Phi(u);
    real zv = inv_Phi(v);
    real s = (zu - rho * zv) / sqrt(1 - rho * rho);
    return Phi(s);
  }

  // ---- discrete–discrete rectangle mass via 5-point Gauss-Legendre on u ----
  real gauss_rect_mass_dd(real uL_in, real uR_in, real vL_in, real vR_in, real rho) {
    real uL = clamp01(uL_in);
    real uR = clamp01(uR_in);
    real vL = clamp01(vL_in);
    real vR = clamp01(vR_in);
    if (uR <= uL) return 0.0;
    // 5-point Gauss–Legendre nodes/weights on [-1,1]
    array[5] real x;
    array[5] real w;
    x[1] = -0.9061798459386640;
    x[2] = -0.5384693101056831;
    x[3] =  0.0;
    x[4] =  0.5384693101056831;
    x[5] =  0.9061798459386640;

    w[1] = 0.2369268850561891;
    w[2] = 0.4786286704993665;
    w[3] = 0.5688888888888889;
    w[4] = 0.4786286704993665;
    w[5] = 0.2369268850561891;

    real mid = 0.5 * (uL + uR);
    real halflen = 0.5 * (uR - uL);
    real acc = 0.0;
    for (i in 1:5) {
      real s = mid + halflen * x[i];
      real integrand = gauss_dC_du(s, vR, rho) - gauss_dC_du(s, vL, rho);
      integrand = fmax(integrand, 0.0);
      acc += w[i] * integrand;
    }
    real mass = halflen * acc;
    return fmax(mass, 1e-300);
  }
}
data {
  int<lower=1> N;
  int<lower=1> P1;
  int<lower=1> P2;
  int<lower=1> P3;
  int<lower=1> P4;
  matrix[N, P1] X1;
  matrix[N, P2] X2;
  matrix[N, P3] X3;
  matrix[N, P4] X4;
  vector<lower=0,upper=1>[N] y;
}
parameters {
  vector[P1] beta1;   // logit(p1)=X1*beta1
  vector[P2] beta2;   // logit(p2)=X2*beta2
  vector[P3] beta3;   // logit(v)=X3*beta3
  vector[P4] beta4;   // log(phi)=X4*beta4
  real rho_raw;       // unconstrained correlation
}
transformed parameters {
  vector[N] p1 = inv_logit(X1 * beta1);
  vector[N] p2 = inv_logit(X2 * beta2);
  vector[N] v  = inv_logit(X3 * beta3);

  vector[N] eta_phi = X4 * beta4;
  vector<lower=1e-6>[N] phi = log1p_exp(eta_phi) + 1e-4;

  vector[N] mu;
  for (n in 1:N) {
    real p1n = fmin(fmax(p1[n], 1e-6), 1 - 1e-6);  
    real p2n = fmin(fmax(p2[n], 1e-6), 1 - 1e-6);
    real denom = fmax(1e-6, 1.0 - p2n);
    real num   = (v[n] / p1n) - p2n;
    real mu_raw = num / denom;
    mu[n] = fmin(fmax(mu_raw, 1e-4), 1 - 1e-4); 

  }

  real rho = tanh(rho_raw);
}
model {
  // priors
  beta1 ~ normal(0, 2);
  beta2 ~ normal(0, 2);
  beta3 ~ normal(0, 0.5);
  beta4 ~ normal(0, 0.5);
  rho_raw ~ normal(0, 0.1);

  // marginal ZOIB likelihood
  for (n in 1:N) {
    if (y[n] == 0) {
      target += bernoulli_lpmf(0 | p1[n]);
    } else if (y[n] == 1) {
      target += bernoulli_lpmf(1 | p1[n]) + bernoulli_lpmf(1 | p2[n]);
    } else {
      target += bernoulli_lpmf(1 | p1[n]);
      target += bernoulli_lpmf(0 | p2[n]);
      target += beta_lpdf(y[n] | mu[n] * phi[n], (1 - mu[n]) * phi[n]);
    }
  }

  // pair-copula (Gaussian) contributions
  {
    vector[N] u;
    vector[N] u_left;
    vector[N] u_right;

    for (n in 1:N) {
      if (y[n] == 0) {
        u[n]       = clamp01(1 - p1[n]);
        u_left[n]  = 0.0;
        u_right[n] = clamp01(1 - p1[n]);
      } else if (y[n] == 1) {
        u[n]       = 1.0;
        u_left[n]  = clamp01(1 - p1[n] * p2[n]);
        u_right[n] = 1.0;
      } else {
	real a = fmax(mu[n] * phi[n], 1e-4);
	real b = fmax((1 - mu[n]) * phi[n], 1e-4);
        real Fbeta = beta_cdf(y[n] | a, b);
        u[n] = clamp01(1 - p1[n] + p1[n] * (1 - p2[n]) * Fbeta);
        u_left[n]  = 0.0;
        u_right[n] = 0.0;
      }
    }

    real rho_eff = fmin(0.9999, fmax(-0.9999, rho));

    for (n in 2:N) {
      if (y[n] > 0 && y[n] < 1 && y[n-1] > 0 && y[n-1] < 1) {
        // continuous–continuous
        target += log( gauss_c(u[n], u[n-1], rho_eff) );
      } else if ((y[n] == 0 || y[n] == 1) && (y[n-1] == 0 || y[n-1] == 1)) {
        // discrete–discrete via 1D quadrature on u
        real mass = gauss_rect_mass_dd(u_left[n],  u_right[n],
                                       u_left[n-1],u_right[n-1], rho_eff);
        target += log(mass);
      } else if ((y[n] > 0 && y[n] < 1) && (y[n-1] == 0 || y[n-1] == 1)) {
        // continuous–discrete: dC/dv(u_t, v_right) - dC/dv(u_t, v_left)
        real part = gauss_dC_dv(u[n], u_right[n-1], rho_eff)
                  - gauss_dC_dv(u[n], u_left[n-1],  rho_eff);
        target += log(fmax(part, 1e-300));
      } else { // (y[n] == 0 or 1) && (y[n-1] in (0,1)) : discrete–continuous
        // dC/du(u_right, u_{t-1}) - dC/du(u_left, u_{t-1})
        real part = gauss_dC_du(u_right[n], u[n-1], rho_eff)
                  - gauss_dC_du(u_left[n],  u[n-1], rho_eff);
        target += log(fmax(part, 1e-300));
      }
    }
  }
}
generated quantities {
  vector[N] y_rep;
  vector[N] log_lik;

  // ---- 首先重建 p1,p2,mu,phi 已有，无需重复 ----
  // 生成 y_rep（你已有）...
  for (n in 1:N) {
    int d1 = bernoulli_rng(p1[n]);
    if (d1 == 0) y_rep[n] = 0;
    else {
      int d2 = bernoulli_rng(p2[n]);
      if (d2 == 1) y_rep[n] = 1;
      else {
        real a = fmax(mu[n]*phi[n], 1e-4);
        real b = fmax((1 - mu[n])*phi[n], 1e-4);
        y_rep[n] = beta_rng(a, b);
      }
    }
  }

  // ---- 逐点 log-lik ----
  {
    vector[N] u;
    vector[N] u_left;
    vector[N] u_right;

    // 先算每个 t 的边际部分 log p(y_t | θ) 并存起来
    for (n in 1:N) {
      if (y[n] == 0) {
        log_lik[n] = bernoulli_lpmf(0 | p1[n]);
        u[n]       = clamp01(1 - p1[n]);
        u_left[n]  = 0;
        u_right[n] = clamp01(1 - p1[n]);
      } else if (y[n] == 1) {
        log_lik[n] = bernoulli_lpmf(1 | p1[n]) + bernoulli_lpmf(1 | p2[n]);
        u[n]       = 1;
        u_left[n]  = clamp01(1 - p1[n]*p2[n]);
        u_right[n] = 1;
      } else {
        real a = fmax(mu[n]*phi[n], 1e-4);
        real b = fmax((1 - mu[n])*phi[n], 1e-4);
        log_lik[n] = bernoulli_lpmf(1 | p1[n]) +
                     bernoulli_lpmf(0 | p2[n]) +
                     beta_lpdf(y[n] | a, b);

        real Fbeta = beta_cdf(y[n] | a, b);
        u[n]       = clamp01(1 - p1[n] + p1[n]*(1 - p2[n]) * Fbeta);
        u_left[n]  = 0;
        u_right[n] = 0;
      }
    }

    // 加上 pair-copula 项
    real rho_eff = fmin(0.9999, fmax(-0.9999, tanh(rho_raw)));
    for (n in 2:N) {
      if (y[n] > 0 && y[n] < 1 && y[n-1] > 0 && y[n-1] < 1) {
        log_lik[n] += log( gauss_c(u[n], u[n-1], rho_eff) );
      } else if ((y[n] == 0 || y[n] == 1) && (y[n-1] == 0 || y[n-1] == 1)) {
        real mass = gauss_rect_mass_dd(u_left[n],  u_right[n],
                                       u_left[n-1],u_right[n-1], rho_eff);
        log_lik[n] += log(mass);
      } else if (y[n] > 0 && y[n] < 1) {
        real part = gauss_dC_dv(u[n], u_right[n-1], rho_eff)
                  - gauss_dC_dv(u[n], u_left[n-1],  rho_eff);
        log_lik[n] += log(fmax(part, 1e-300));
      } else {
        real part = gauss_dC_du(u_right[n], u[n-1], rho_eff)
                  - gauss_dC_du(u_left[n],  u[n-1], rho_eff);
        log_lik[n] += log(fmax(part, 1e-300));
      }
    }
  }
}

