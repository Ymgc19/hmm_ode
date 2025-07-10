data {
  int T;                       // 時系列の長さ
  real delta_t;               // Δt 
  vector[T] y;                // 大学進学率（目的変数）
  vector[T] x1;               // 高卒と大卒の所得差
  vector[T] x2;               // 有効求人倍率
  vector[T] x3;               // 完全失業率
}

parameters {
  vector[T] b0;               // 切片の係数（時変）
  vector[T] b1;               // 説明変数1の係数（時変）
  vector[T] b2;               // 説明変数2の係数（時変）
  vector[T] b3;               // 説明変数3の係数（時変）
  real<lower=0> t0;           // 切片の変動幅
  real<lower=0> t1;           // 係数1の変動幅
  real<lower=0> t2;           // 係数2の変動幅
  real<lower=0> t3;           // 係数3の変動幅
  real<lower=0> v;            // 観測誤差の標準偏差
}

transformed parameters {
  vector[T] alpha;  // 離散化された状態の予測値
  alpha[1] = y[1];  // 初期値（観測値から与える）

  for (i in 2:T) {
    real r = b0[i-1] + b1[i-1]*x1[i-1] + b2[i-1]*x2[i-1] + b3[i-1]*x3[i-1];
    alpha[i] = y[i-1] + delta_t * y[i-1] * (r - y[i-1]*r);
  }
}

model {
  t0 ~ exponential(1);
  t1 ~ exponential(1);
  t2 ~ exponential(1);
  t3 ~ exponential(1);
  v ~ exponential(1);

  for (i in 2:T) {
    b0[i] ~ normal(b0[i-1], t0);
    b1[i] ~ normal(b1[i-1], t1);
    b2[i] ~ normal(b2[i-1], t2);
    b3[i] ~ normal(b3[i-1], t3);
  }

  for (i in 2:T) {
    y[i] ~ normal(alpha[i], v);
  }
}
