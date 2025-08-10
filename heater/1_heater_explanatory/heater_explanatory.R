library(rstan)
library(cmdstanr)

# データの読み込み
source("heater/heater_data_handling.R")

# モデルの読み込み
model <- cmdstan_model("heater/1_heater_explanatory/heater_explanatory.stan")




# ======================================================== #
# ===================== K = 1の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 1,
  delta_t = 1,
  y = df$heater_preveil*.01,
  kero = df$kerosene_price*.01
)
data

# MCMC回す
fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)


# 結果を見る
fit$cmdstan_summary()




# 潜在状態を計算する
compute_posterior_states <- function(y, kero, delta_t, pi, A, b0, b1, sigma) {
  T_graph <- length(y) - 1  # y[1] to y[T+1]
  K <- length(pi)
  
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)
  
  # ----- Forward pass -----
  for (k in 1:K) {
#    mu <- y[1] + delta_t * (b1[k]*kero[1]) * y[1] * (1-y[1])
    mu <- y[1] + delta_t * (b0[k] + b1[k]*kero[1]) * y[1] * (1-y[1])
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in 2:T_graph) {
    for (k in 1:K) {
#      mu <- y[t_r] + delta_t * (b1[k]*kero[t_r]) * y[t_r] * (1-y[t_r])
      mu <- y[t_r] + delta_t * (b0[k] + b1[k]*kero[t_r]) * y[t_r] * (1-y[t_r])
      temp <- sapply(1:K, function(j) {
        log_alpha[t_r - 1, j] + log(A[j, k])
      })
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mean = mu, sd = sigma[k], log = TRUE)
    }
  }
  
  # ----- Backward pass -----
  for (k in 1:K) {
#    mu <- y[T_graph] + delta_t * b[k]
    mu <- y[T_graph] + delta_t * (b0[k] + b1[k]*kero[T_graph]) * y[T_graph] * (1-y[T_graph])
    log_beta[T_graph, k] <- dnorm(y[T_graph + 1], mean = mu, sd = sigma[k], log = TRUE)
  }
  
  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
#        mu <- y[t_r + 1] + delta_t * b[j]
        mu <- y[t_r + 1] + delta_t * (b0[k] + b1[k]*kero[t_r + 1]) * y[t_r + 1] * (1-y[t_r + 1])
        log(A[k, j]) + dnorm(y[t_r + 2], mean = mu, sd = sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }
  
  # ----- Posterior state probabilities -----
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t)
    gamma[t_r, ] <- exp(log_gamma_t)
  }
  
  return(gamma)
}

# 安定な log-sum-exp
log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b0
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, df$kerosene_price*.01, delta_t = 1, pi, A, b0, b1, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table



# 予測値を入れるリスト
y <- df$heater_preveil*.01
kero <- df$kerosene_price*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_lower <- c(NA)
y_upper <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, kero, b0, b1, delta_t = 1) {
  return(y + y*(1-y)*delta_t*(b0+b1*kero))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], kero[i], b0[s[i]], b1[s[i]]))
  y_lower <- c(y_lower, y_pred[i+1] - 2*sigma[s[i]])
  y_upper <- c(y_upper, y_pred[i+1] + 2*sigma[s[i]])
}

plot(df$year, y_pred)

# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Time") +
  scale_fill_manual(values = c("royalblue3", "yellow", "turquoise"), 
                    labels = c("0")
  ) +
  geom_ribbon(
    aes(x = df$year, ymin = c(y_lower), ymax = c(y_upper)), 
    alpha = .25
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

























# ======================================================== #
# ===================== K = 2の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 2,
  delta_t = 1,
  y = df$heater_preveil*.01,
  kero = df$kerosene_price*.01
)
data

# sMCMC回す
fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)


# 結果を見る
fit$cmdstan_summary()


# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b0
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, df$kerosene_price*.01, delta_t = 1, pi, A, b0, b1, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table


# 予測値を入れるリスト
y <- df$heater_preveil*.01
kero <- df$kerosene_price*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_lower <- c(NA)
y_upper <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, kero, b0, b1, delta_t = 1) {
  return(y + y*(1-y)*delta_t*(b0+b1*kero))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], kero[i], b0[s[i]], b1[s[i]]))
  y_lower <- c(y_lower, y_pred[i+1] - 2*sigma[s[i]])
  y_upper <- c(y_upper, y_pred[i+1] + 2*sigma[s[i]])
}

plot(df$year, y_pred)


# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Year") +
  scale_fill_manual(values = c("yellow", "turquoise"), 
                    labels = c("0", "1")
  ) +
  geom_ribbon(
    aes(x = df$year, ymin = c(y_lower), ymax = c(y_upper)), 
    alpha = .25
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

















# ======================================================== #
# ===================== K = 3の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 3,
  delta_t = 1,
  y = df$heater_preveil*.01,
  kero = df$kerosene_price*.01
)
data

# sMCMC回す
fit <- model$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)


# 結果を見る
fit$cmdstan_summary()


# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b0
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, df$kerosene_price*.01, delta_t = 1, pi, A, b0, b1, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table



# 予測値を入れるリスト
y <- df$heater_preveil*.01
kero <- df$kerosene_price*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_lower <- c(NA)
y_upper <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, kero, b0, b1, delta_t = 1) {
  return(y + y*(1-y)*delta_t*(b0+b1*kero))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], kero[i], b0[s[i]], b1[s[i]]))
  y_lower <- c(y_lower, y_pred[i+1] - 2*sigma[s[i]])
  y_upper <- c(y_upper, y_pred[i+1] + 2*sigma[s[i]])
}


# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Year") +
  scale_fill_manual(values = c("royalblue3", "yellow", "turquoise"), 
                    labels = c("0", "1", "2")
  ) +
  geom_ribbon(
    aes(x = df$year, ymin = c(y_lower), ymax = c(y_upper)), 
    alpha = .25
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

b0; b1
s

# ===== それぞれの潜在状態の切片と係数に基づいた直線の描画 ===== #
p <- seq(df$kerosene_price %>% min, df$kerosene_price %>% max)
s1 <- b0[1] + b1[1]*p*.01
s2 <- b0[2] + b1[2]*p*.01
s3 <- b0[3] + b1[3]*p*.01

# 点を打つためのコード
y <- c()
for (i in 1:length(s)){
  if (s[i] == 1){
    y <- c(y, b0[1] + b1[1]*df$kerosene_price[i]*.01)
  }
  else if (s[i] == 2){
    y <- c(y, b0[2] + b1[2]*df$kerosene_price[i]*.01)
  }
  else{
    y <- c(y, b0[3] + b1[3]*df$kerosene_price[i]*.01)
  }
}
ggplot() +
  geom_line(aes(p, s1), color = "royalblue3", lwd = 3) +
  geom_line(aes(p, s2), color = "yellow", lwd = 3) +
  geom_line(aes(p, s3), color = "turquoise", lwd = 3) +
  geom_line(aes(p, s1), lwd = .25) +
  geom_line(aes(p, s2), lwd = .25) +
  geom_line(aes(p, s3),  lwd = .25) +
  geom_point(aes(df$kerosene_price[1:39], y), size = 1.5) +
  theme_minimal() +
  geom_hline(yintercept = 0, color = "black", linetype = "dotted", lwd = .75) +
  labs(x = "Kerosene price", y = "Utility fan heater possession (lambda)")









# ======================================================== #
# ===================== K = 4の場合 ======================= #
# ======================================================== #
# データの用意
data <- list(
  t = nrow(df),
  K = 4,
  delta_t = 1,
  y = df$heater_preveil*.01,
  kero = df$kerosene_price*.01
)
data

# sMCMC回す
fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)


# 結果を見る
fit$cmdstan_summary()


# パラメタを取り出す
pi_draws <- fit$draws(variables = "pi", format = "draws_matrix")
pi <- colMeans(pi_draws) %>% as.vector()
pi

library(posterior)
A_mat <- as_draws_matrix(fit$draws("A"))
A_names <- colnames(A_mat)
K <- sqrt(length(grep("^A\\[", A_names)))
A <- matrix(colMeans(A_mat), nrow = K, byrow = TRUE)
A

b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b0
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1

sigma <- fit$draws(variables = c("sigma"), format = "draws_matrix") %>% colMeans() %>% as.vector()
sigma

gamma <- compute_posterior_states(df$heater_preveil*.01, df$kerosene_price*.01, delta_t = 1, pi, A, b0, b1, sigma)

# 潜在状態の決定
s <- c()
for (i in 1:nrow(gamma)){
  s <- c(s, gamma[i,] %>% which.max)
}
s
s %>% table



# 予測値を入れるリスト
y <- df$heater_preveil*.01
kero <- df$kerosene_price*.01
y_pred <- c(y[1])
y_pred_mu <- c(y[1])
y_lower <- c(NA)
y_upper <- c(NA)
# 予測値を計算
# 基礎方程式
de <- function(y, kero, b0, b1, delta_t = 1) {
  return(y + y*(1-y)*delta_t*(b0+b1*kero))
}
# 平均ではなく中央値とか？
for (i in 1:length(s)){
  y_pred <- c(y_pred, de(y[i], kero[i], b0[s[i]], b1[s[i]]))
  y_lower <- c(y_lower, y_pred[i+1] - 2*sigma[s[i]])
  y_upper <- c(y_upper, y_pred[i+1] + 2*sigma[s[i]])
}


# 例: gamma_matrix が T行4列の数値行列
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2", "state3")
gamma_df$time <- 1985:2023

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2", "state3")))

ggplot() +
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  labs(y = "Probability / Penetration rate of fan heater", x = "Year") +
  scale_fill_manual(values = c("springgreen", "royalblue3", "yellow", "turquoise"), 
                    labels = c("0", "1", "2", "3")
  ) +
  geom_ribbon(
    aes(x = df$year, ymin = c(y_lower), ymax = c(y_upper)), 
    alpha = .25
  ) +
  geom_line(
    data = df, color = "black", lwd = .5,
    aes(year, heater_preveil*.01, color = c(s, NA))
  ) +
  geom_point(
    aes(df$year, y_pred),
    shape = 4
  ) +
  geom_point(
    data = df, 
    aes(year, heater_preveil*.01)
  ) +
  theme_minimal() +
  theme(legend.position = "none")

b0; b1
s














# 事後分布の描画
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")
draws <- as_draws_df(fit$draws())
mcmc_combo(draws, pars = c("b0[1]", "b0[2]", "b0[3]", "b1[1]", "b1[2]", "b1[3]"))
mcmc_intervals(draws, pars = c("b0[1]", "b0[2]", "b0[3]", "b1[1]", "b1[2]", "b1[3]"))
mcmc_areas(draws, pars = c("b0[1]", "b0[2]", "b0[3]", "b1[1]", "b1[2]", "b1[3]"), prob = 0.9)
