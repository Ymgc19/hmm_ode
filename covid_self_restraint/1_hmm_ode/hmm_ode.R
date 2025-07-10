# これは，コロナの自粛率のモデルです
# データの下処理
source("covid_self_restraint/data_preprocessing.R")






# ====================================================================================================== #
# =========================================== MCMCを回すぜ！ =========================================== #
# ====================================================================================================== #
# モデルの読み込み
model <- cmdstan_model("covid_self_restraint/1_hmm_ode/hmm_ode.stan")




# ===== データの調整 ===== #
df_for_model <- df_pref %>% 
  filter(!is.na(M30_smooth)) %>% 
  # 便宜的に自粛率が0を下回っているところを0.001に置換
  mutate(M30_smooth = if_else(M30_smooth < 0, 0.01, M30_smooth)) %>% 
  filter(date >= "2020-01-16", #date <= "2020-07-10", 
         )

df_for_model %>% glimpse()
df_for_model$M30_smooth %>% plot
df_for_model$Tokyo %>% plot


#df_for_model %>% view()
# モデルに渡すデータの整備
y <- df_for_model$M30_smooth # 自粛率
#y <- df_for_model$M30 # 自粛率
D <- df_for_model$emergency # 緊急事態宣言ダミー
tau <- df_for_model$pass_days # 宣言後の経過日数
infectious <- df_for_model$Tokyo # 新規感染者数


data <- list(
  t = length(D),
  K = 2,
  Dt = 1,
  y = y,
  D = D,
  tau = tau,
  infectious = infectious
)
data


fit <- model$sample(
  data = data,
  chains = 1,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 10000,
  refresh = 1
)

# 結果を見る
fit$cmdstan_summary()




# ======================= 推定されたパラメタの下での予測値 =========================== #
compute_posterior_states <- function(y, D, tau, pi, A, intercept, rho, phi, delta, sigma, Dt = 1) {
  T_graph <- length(y) - 1 
  K <- length(pi) 
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K) 
  log_beta <- matrix(NA, nrow = T_graph, ncol = K) 
  # forward
  for (k in 1:K) {
    mu <- y[1] + Dt*y[1]*(1-y[1])*(intercept[k] + rho[k] - phi[k]*D[1]*(delta[k]^tau[1]))*(-1)
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  for (t_r in 2:T_graph) { 
    for (k in 1:K) { 
      mu <- y[t_r] + Dt*y[t_r]*(1-y[t_r])*(intercept[k] + rho[k] - phi[k]*D[t_r]*(delta[k]^tau[t_r]))*(-1)
      temp <- sapply(1:K, function(j) { 
        log_alpha[t_r - 1, j] + log(A[j, k]) 
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mu, sigma[k], log = TRUE)
    }
  }
  # backward
  for (k in 1:K) {
    mu_final_beta <- y[T_graph] + Dt*y[T_graph]*(1-y[T_graph])*(intercept[k] + rho[k] - phi[k]*D[1]*(delta[k]^tau[1]))*(-1)
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mu_final_beta, sigma[k], log = TRUE)
  }
  for (t_r in (T_graph - 1):1) { 
    for (k in 1:K) { 
      temp <- sapply(1:K, function(j) { 
        mu_inner <- y[t_r] + Dt*y[t_r]*(1-y[t_r])*(intercept[k] + rho[k] - phi[k]*D[t_r]*(delta[k]^tau[t_r]))*(-1)
        log(A[k, j]) + dnorm(y[t_r + 1], mu_inner, sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }
  # posterior copmputatation
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t) # 正規化
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  return(gamma)
}

log_sum_exp <- function(x) {
  m <- max(x)
  return(m + log(sum(exp(x - m))))
}












pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
intercept <- fit$draws(variables = c("intercept"), format = "draws_matrix") %>% colMeans() %>% as.vector()
rho <- fit$draws(variables = c("rho"), format = "draws_matrix") %>% colMeans() %>% as.vector()
phi <- fit$draws(variables = c("phi"), format = "draws_matrix") %>% colMeans() %>% as.vector()
delta <- fit$draws(variables = c("delta"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
phi; intercept; rho; delta; sigma


gamma <- compute_posterior_states(y, D, tau, pi, A, intercept, rho, phi, delta, sigma, Dt = 1)

gamma

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table



# 予測値を入れるリスト
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_low <- c(NA)
y_high <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, phi, intercept, rho, delta, Dt = 1, D, tau) {
  return(y + Dt*y*(1-y)*(intercept + rho - phi*D*(delta^tau))*(-1))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], phi[s[i]], intercept[s[i]], rho[s[i]], delta[s[i]], Dt = 1, D[i], tau[i]))
  y_pred_mu <- c(y_pred_mu, de(y_pred_mu[i], phi[s[i]], intercept[s[i]], rho[s[i]], delta[s[i]], Dt = 1, D[i], tau[i]))
  y_low <- c(y_low, y_pred_mu[i+1] - 2*sigma[s[i]])
  y_high <- c(y_high, y_pred_mu[i+1] + 2*sigma[s[i]])
}




ggplot() +
  # 実測値
  geom_line(
    aes(x = 0:(length(y)-1), y = y),
    color = "black"
  ) +
  geom_point(
    aes(x = 0:(length(y)-1), y = y),
    color = "black", shape = 15
  ) +
  # 微分方程式
  geom_ribbon(
    aes(x = 0:(length(y)-1),
        ymin = y_low, ymax = y_high),
    fill = "black", alpha = .125
  ) +
  geom_line(
    aes(
      x = 0:(length(y)-1), y_pred_mu,
      color = c(s, NA)
    )
  ) +
  geom_point(
    aes(
      x = 0:(length(y)-1), y_pred_mu,
      color = c(s, NA)
    ), shape = 3
  ) +
  labs(x = "days", y = "jishuku") +
  scale_colour_gradientn(colours = c("tomato", "orange", "turquoise", "royalblue", "purple")) +
  theme_minimal()



# 予測値と実測値の差
ggplot() +
  geom_point(
    aes(x = 0:(length(y)-1), y = y - y_pred),
    color = "black", shape = 3
  ) +
  geom_abline(intercept = 0, slope = 0) +
  theme_minimal()


































