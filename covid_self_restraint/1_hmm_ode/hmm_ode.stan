data {
  int t;               // 時系列の長さ
  int K;               // 潜在変数の数
  real Dt;                      // 時間幅
  vector[t] y;                  // 観測された自粛率（0〜1）
  vector[t] D;                  // 緊急事態宣言ダミー（0/1）
  vector[t] tau;                // 宣言解除後の経過日数
  vector[t] infectious;         // 感染者数
}

parameters {
  simplex[K] pi;                        // 初期状態分布
  array[K] simplex[K] A;                // 遷移確率行列
  ordered[K] phi;
  vector[K] intercept;
  vector[K] rho;                         // 外出の期待利得
  array[K] real<lower=0, upper=1> delta; // 割引因子
  array[K] real<lower=10e-6> sigma;          // 観測ノイズの標準偏差
}


model {
  // 事前分布
  phi ~ normal(0,.5);
  rho ~ normal(0,.5);
  sigma ~ normal(0, 0.1);
  
  matrix[t-1, K] log_alpha;

  // 初期ステップ（t=2を予測）
  for (k in 1:K) {
    real decay = pow(delta[k], tau[1]);
    real drift = Dt * y[1] * (1 - y[1]) * (intercept[k] + rho[k]*infectious[1] - phi[k]* D[1] * decay);
    log_alpha[1, k] = log(pi[k]) + normal_lpdf(y[2] | y[1] - drift, sigma[k]);
  }

  // 再帰ステップ
  for (time in 2:(t-1)) {
    for (k in 1:K) {
      vector[K] temp;
      // 次時点の平均値を計算
      real decay = pow(delta[k], tau[time]);
      real drift = Dt * y[time] * (1 - y[time]) * (intercept[k] + rho[k]*infectious[time] - phi[k]* D[time] * decay);
      // t時点に状態jにいる確率を計算（対数）
      for (j in 1:K) {
        temp[j] = log_alpha[time - 1, j] + log(A[j, k]);
      }
      log_alpha[time, k] = log_sum_exp(temp) + normal_lpdf(y[time+1] | y[time] - drift, sigma[k]);
    }
  }

  target += log_sum_exp(log_alpha[t-1]);
}

