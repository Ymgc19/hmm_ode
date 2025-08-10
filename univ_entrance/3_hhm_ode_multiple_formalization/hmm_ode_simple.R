# これは，隠れマルコフ微分方程式をを用いて，大学進学率を扱うファイル
# 説明変数を使話ないバージョン

# データの読み込みと前処理
source("univ_entrance/data_preprocessing.R")

# ===== MCMCのためのデータ ===== #
data = list(
  t = nrow(df),
  K = 3,
  delta_t = 1,
  y = df$enter_univ_rate,
  salary_diff = df$salary_diff,
  jobs_to_applicants_ratio = df$jobs_to_applicants_ratio,
  unemployment_rate = df$unemployment_rate
)




# ========== モデルを読み込む ========== #
model <- cmdstan_model("univ_entrance/3_hhm_ode_multiple_formalization/hmm_ode_multiple_formalization.stan")


# ========== サンプリング ========== #
fit <- model$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  thin = 5,
  refresh = 1000
)

# 結果を見る
options(scipen = 4)
fit$cmdstan_summary()


# お絵描き
#draws <- fit$draws(
#  variables = c("H", "sigma"),
#  format = "draws_array"
#)

# traceplot
#color_scheme_set("mix-blue-red")  # 色設定（任意）
#mcmc_trace(draws)

# posterior
#mcmc_dens_overlay(draws)




# ===== ログサム指数関数 =====
log_sum_exp <- function(x) {
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}

# ===== 潜在状態の事後分布の計算 =====
compute_posterior_states <- function(y, delta_t, pi, A, i, b, lambda, sigma, 
                                     salary_diff, jobs_to_applicants_ratio, unemployment_rate) {
  T_graph <- length(y) - 1
  K <- length(pi)
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)

  # forward
  for (k in 1:K) {
    mu <- switch(k,
      y[1] + delta_t * (i[1] + b[1] * salary_diff[1]) * y[1] * (1 - y[1]),
      y[1] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[1]) * y[1] * (1 - y[1]),
      y[1] + delta_t * (i[3] + b[3] * unemployment_rate[1]) * y[1] * (1 - y[1]),
      y[1] + delta_t * lambda * y[1] * (1 - y[1])
    )
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }

  for (t_r in 2:T_graph) {
    for (k in 1:K) {
      mu <- switch(k,
        y[t_r] + delta_t * (i[1] + b[1] * salary_diff[t_r]) * y[t_r] * (1 - y[t_r]),
        y[t_r] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[t_r]) * y[t_r] * (1 - y[t_r]),
        y[t_r] + delta_t * (i[3] + b[3] * unemployment_rate[t_r]) * y[t_r] * (1 - y[t_r]),
        y[t_r] + delta_t * lambda * y[t_r] * (1 - y[t_r])
      )
      temp <- sapply(1:K, function(j) {
        log_alpha[t_r - 1, j] + log(A[j, k])
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mu, sigma[k], log = TRUE)
    }
  }

  # backward
  for (k in 1:K) {
    mu <- switch(k,
      y[T_graph] + delta_t * (i[1] + b[1] * salary_diff[T_graph]) * y[T_graph] * (1 - y[T_graph]),
      y[T_graph] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[T_graph]) * y[T_graph] * (1 - y[T_graph]),
      y[T_graph] + delta_t * (i[3] + b[3] * unemployment_rate[T_graph]) * y[T_graph] * (1 - y[T_graph]),
      y[T_graph] + delta_t * lambda * y[T_graph] * (1 - y[T_graph])
    )
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mu, sigma[k], log = TRUE)
  }

  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
        mu <- switch(j,
          y[t_r] + delta_t * (i[1] + b[1] * salary_diff[t_r]) * y[t_r] * (1 - y[t_r]),
          y[t_r] + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio[t_r]) * y[t_r] * (1 - y[t_r]),
          y[t_r] + delta_t * (i[3] + b[3] * unemployment_rate[t_r]) * y[t_r] * (1 - y[t_r]),
          y[t_r] + delta_t * lambda * y[t_r] * (1 - y[t_r])
        )
        log(A[k, j]) + dnorm(y[t_r + 1], mu, sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }

  # posterior distribution of states
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t)
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  return(gamma)
}





















y <- df$enter_univ_rate

# パラメータ抽出
pi <- fit$draws("pi", format = "draws_matrix") |> colMeans()
A <- fit$draws("A", format = "draws_matrix") |> colMeans() |> matrix(nrow = 4, byrow = TRUE)
i <- fit$draws("i", format = "draws_matrix") |> colMeans()
b <- fit$draws("b", format = "draws_matrix") |> colMeans()
lambda <- fit$draws("lambda", format = "draws_matrix") |> colMeans()
sigma <- fit$draws("sigma", format = "draws_matrix") |> colMeans()

gamma <- compute_posterior_states(
  y = y,
  delta_t = 1,
  pi = pi,
  A = A,
  i = i,
  b = b,
  lambda = lambda,
  sigma = sigma,
  salary_diff = df$salary_diff,
  jobs_to_applicants_ratio = df$jobs_to_applicants_ratio,
  unemployment_rate = df$unemployment_rate
)
gamma























# 潜在状態の決定
s <- apply(gamma, 1, which.max)
s
s %>% table


# 状態ごとの微分方程式
de <- function(y, k, i, b, lambda, salary_diff, jobs_to_applicants_ratio, unemployment_rate, delta_t = 1) {
  mu <- switch(k,
    y + delta_t * (i[1] + b[1] * salary_diff) * y * (1 - y),
    y + delta_t * (i[2] + b[2] * jobs_to_applicants_ratio) * y * (1 - y),
    y + delta_t * (i[3] + b[3] * unemployment_rate) * y * (1 - y),
    y + delta_t * lambda * y * (1 - y)
  )
  return(mu)
}

# 予測値計算
y_pred <- c(y[1])
y_low <- c()
y_high <- c()

for (j in 1:length(s)) {
  mu <- de(
    y[j], s[j], i, b, lambda,
    df$salary_diff[j],
    df$jobs_to_applicants_ratio[j],
    df$unemployment_rate[j]
  )
  y_pred <- c(y_pred, mu)
  y_low <- c(y_low, mu - 2 * sigma[s[j]])
  y_high <- c(y_high, mu + 2 * sigma[s[j]])
}

# 可視化
library(ggplot2)
ggplot() +
  geom_ribbon(
    aes(x = df$year[2:nrow(df)],
        ymin = y_low, ymax = y_high),
    fill = "black", alpha = 0.125
  ) +
  geom_line(aes(x = df$year, y = df$enter_univ_rate), color = "black") +
  geom_point(aes(x = df$year, y = df$enter_univ_rate), color = "black", shape = 15) +
  geom_line(aes(x = df$year[1:length(y_pred)], y = y_pred, color = c(s, NA))) +
  geom_point(aes(x = df$year[1:length(y_pred)], y = y_pred, color = c(s, NA)), shape = 3) +
  labs(x = "Year", y = "University Enrollment Rate") +
  theme_minimal() +
  scale_color_gradientn(colours = c("purple", "turquoise", "tomato")) +
  theme(legend.position = "none")

s
b
i
