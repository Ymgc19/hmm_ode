source("recipients/data_handling.R")
library(rstan)
library(cmdstanr)

model <- cmdstan_model("recipients/2_hmm_ode/hmm_ode_umemployment_model.stan")

K <- 3
N <- 5

# データの用意
data <- list(
  t = nrow(df),
  K = K,
  delta_t = 1,
  y = df$recipients_density,
  p = df$umemploy,
  N = N
)
data


fit <- model$sample(
  data = data,
  chains = 2,
  parallel_chains = 4,
  iter_warmup = 3000,
  iter_sampling = 3000,
  refresh = 100
)
fit$cmdstan_summary()


# モデル保存
#saveRDS(fit, file = "recipients/2_hmm_ode/fit.rds")
# モデル読み込み
#fit <- readRDS("recipients/2_hmm_ode/fit.rds")

1+1










# ==================== 事後分布の計算 ====================== #
calc_mu <- function(k, y_prev, p_prev, N, lambda, u_recipients, delta_t) {
  if (k == 1) {
    q <- 1 - p_prev^(1.0 / N)
    u <- 0
    for (m in 0:N) {
      bin_prob <- exp(dbinom(m, N, q, log = TRUE))
      u <- u + ((lambda[k]+1)^m - 1) * bin_prob
    }
    y_prev + delta_t * y_prev * (1 - y_prev) *
      (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
  } else if (k == 2) {
    q <- 1 - p_prev^(1.0 / N)
    u <- 0
    for (m in 0:N) {
      bin_prob <- exp(dbinom(m, N, q, log = TRUE))
      u <- u + (lambda[k] * log(m + 1)) * bin_prob
    }
    y_prev + delta_t * y_prev * (1 - y_prev) *
      (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
  } else if (k == 3) {
    q <- 1 - p_prev^(1.0 / N)
    u <- 0
    for (m in 0:N) {
      bin_prob <- exp(dbinom(m, N, q, log = TRUE))
      u <- u + lambda[k] * m * bin_prob
    }
    y_prev + delta_t * y_prev * (1 - y_prev) *
      (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
  } else {
    y_prev + 100 * delta_t * y_prev * (1 - y_prev)
  }
}


log_sum_exp <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

compute_posterior_states <- function(y, p, delta_t, pi, A, lambda, u_recipients, sigma, N) {
  T_graph <- length(y) - 1
  K <- length(pi)
  log_alpha <- matrix(NA, nrow = T_graph, ncol = K)
  log_beta <- matrix(NA, nrow = T_graph, ncol = K)

  # ---- μ計算関数 ----
  calc_mu <- function(k, y_prev, p_prev, N, lambda, u_recipients, delta_t) {
    if (k == 1) {
      q <- 1 - p_prev^(1.0 / N)
      u <- 0
      for (m in 0:N) {
        bin_prob <- exp(dbinom(m, N, q, log = TRUE))
        u <- u + ((lambda[k]+1)^m - 1) * bin_prob
      }
      y_prev + delta_t * y_prev * (1 - y_prev) * 
        (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
    } else if (k == 2) {
      q <- 1 - p_prev^(1.0 / N)
      u <- 0
      for (m in 0:N) {
        bin_prob <- exp(dbinom(m, N, q, log = TRUE))
        u <- u + lambda[k] * log(m + 1) * bin_prob
      }
      y_prev + delta_t * y_prev * (1 - y_prev) * 
        (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
    } else if (k == 3) {
      q <- 1 - p_prev^(1.0 / N)
      u <- 0
      for (m in 0:N) {
        bin_prob <- exp(dbinom(m, N, q, log = TRUE))
        u <- u + lambda[k] * m * bin_prob
      }
      y_prev + delta_t * y_prev * (1 - y_prev) * 
        (u_recipients[k] - (y_prev * u_recipients[k] + (1 - y_prev) * u))
    } else {
      y_prev + 100 * delta_t * y_prev * (1 - y_prev)
    }
  }

  # ---- forward ----
  for (k in 1:K) {
    mu <- calc_mu(k, y[1], p[1], N, lambda, u_recipients, delta_t)
    log_alpha[1, k] <- log(pi[k]) + dnorm(y[2], mean = mu, sd = sigma[k], log = TRUE)
  }

  for (t_r in 2:T_graph) {
    for (k in 1:K) {
      mu <- calc_mu(k, y[t_r], p[t_r], N, lambda, u_recipients, delta_t)
      temp <- sapply(1:K, function(j) log_alpha[t_r - 1, j] + log(A[j, k]))
      log_alpha[t_r, k] <- log_sum_exp(temp) + dnorm(y[t_r + 1], mu, sigma[k], log = TRUE)
    }
  }

  # ---- backward ----
  log_beta[T_graph, ] <- 0  # 最終時刻の対数βは0（確率1）

  for (t_r in (T_graph - 1):1) {
    for (k in 1:K) {
      temp <- sapply(1:K, function(j) {
        mu <- calc_mu(j, y[t_r + 1], p[t_r + 1], N, lambda, u_recipients, delta_t)
        log(A[k, j]) + dnorm(y[t_r + 2 - 1], mu, sigma[j], log = TRUE) + log_beta[t_r + 1, j]
      })
      log_beta[t_r, k] <- log_sum_exp(temp)
    }
  }

  # ---- posterior γ ----
  gamma <- matrix(NA, nrow = T_graph, ncol = K)
  for (t_r in 1:T_graph) {
    log_gamma_t <- log_alpha[t_r, ] + log_beta[t_r, ]
    log_gamma_t <- log_gamma_t - log_sum_exp(log_gamma_t)
    gamma[t_r, ] <- exp(log_gamma_t)
  }

  gamma
}




y <- df$recipients_density
p <- df$umemploy

# パラメータ抽出
pi <- fit$draws("pi", format = "draws_matrix") |> colMeans()
A <- fit$draws("A", format = "draws_matrix") |> colMeans() |> matrix(nrow = K, byrow = TRUE)
lambda <- fit$draws("lambda", format = "draws_matrix") |> colMeans()
u_recipients <- fit$draws("u_recipients", format = "draws_matrix") |> colMeans()
sigma <- fit$draws("sigma", format = "draws_matrix") |> colMeans()


lambda
u_recipients


gamma <- compute_posterior_states(y, p, delta_t = 1, pi, A, lambda, u_recipients, sigma, N = N)
gamma

# 潜在状態の決定
s <- apply(gamma, 1, which.max)
s
s %>% table








# ======================= 可視化など ======================== #
# ==== 潜在状態の決定 ====
s <- apply(gamma, 1, which.max)
table(s)
delta_t <- 1



# ==== 予測値計算 ====
y_pred <- c(y[1])
y_low <- c(NA)
y_high <- c(NA)

for (j in 1:length(s)) {
  mu <- calc_mu(
    k = s[j],
    y_prev = y[j],
    p_prev = p[j],
    N = N,
    lambda = lambda,
    u_recipients = u_recipients,
    delta_t = delta_t
  )
  y_pred <- c(y_pred, mu)
  y_low <- c(y_low, mu - 2 * sigma[s[j]])
  y_high <- c(y_high, mu + 2 * sigma[s[j]])
}







# ==== 可視化 ====
# 成功確率の計算
df <- df %>% 
  mutate(
    q = 1 - umemploy^(1/N)
  )

# 潜在変数の事後確率
gamma_df <- as.data.frame(gamma)
colnames(gamma_df) <- c("state0", "state1", "state2")
gamma_df$time <- min(df$year):(max(df$year)-1)

# long形式に変換
gamma_long <- gamma_df %>%
  pivot_longer(cols = starts_with("state"), names_to = "state", values_to = "probability") %>%
  mutate(state = factor(state, levels = c("state0", "state1", "state2")))

# スケーリング関数
scale_factor <- 0.019 - 0.005  # 0.014
offset <- 0.005

ggplot() +
  # 事後確率
  geom_area(
    data = gamma_long, 
    aes(x = time, y = probability, fill = state), 
    alpha = .5
  ) +
  # 実測値
  geom_point(
    aes(x = df$year, y = (df$recipients_density - offset) / scale_factor),
    color = "black"
  ) + 
  geom_line(
    aes(x = df$year, y = (df$recipients_density - offset) / scale_factor),
    color = "black"
  ) + 
  # 成功確率
  geom_line(
    aes(df$year, df$q, color = df$q)
  ) +
  geom_point(
    aes(df$year, df$q, color = df$q)
  ) +
  annotate(
    "text", 
    x = 1995, y = .575,         # 文字を置く位置
    label = "↓求職成功確率 (q_t)",  # 表示する文字
    size = 3.5,              # フォントサイズ
    color = "black",
    fontface = "bold"
  ) +
  # 予測値
  geom_point(
    aes(x = df$year, y = (y_pred - offset) / scale_factor),
    shape = 4
  ) + 
  geom_ribbon(
    aes(
      x = df$year, 
      ymin = (y_low - offset) / scale_factor,
      ymax = (y_high - offset) / scale_factor
    ), 
    alpha = .25
  ) +
  # 雑多な設定
  scale_fill_manual(values = c("royalblue3", "turquoise", "yellow"), 
                    labels = c("0", "1", "2"),
                    guide = "none") +
  scale_color_gradient2(
    high = "hotpink", mid = "royalblue", low = "royalblue", midpoint = .5
  ) +
  scale_y_continuous(
    name = "Probability",                     
    sec.axis = sec_axis(~ . * scale_factor + offset, 
                        name = "Recipients density",
                        breaks = seq(0.005, 0.019, by = 0.002))
  ) +
  coord_cartesian(ylim = c(0, 1)) +  # 左軸は0～1
  theme_minimal() +
  labs(x = "Year", color = "q_t") +
  theme(legend.position = "none")

1 + 1






# ========== 効用関数の可視化 ========== #
m <- 0:N
# ===== 指数型の関数 ===== #
u1 <- ((lambda[1]+1)^m - 1)

# ===== 対数型の関数 ===== #
u2 <- (lambda[2] * log(m + 1))

# ===== 比例型の関数 ===== #
u3 <- lambda[3] * m

ggplot() +
  # u1
  geom_point(aes(m, u1), color = "royalblue3", size = 2) +
  geom_line(aes(m, u1), color = "royalblue3", lwd = 1, linetype = "dotted") +
  geom_line(aes(m, u1), lwd = .2) +
  # u2
  geom_point(aes(m, u2), color = "turquoise", size = 4) +
  geom_line(aes(m, u2), color = "turquoise", lwd = 3) +
  geom_line(aes(m, u2), lwd = .2) +
  # u3
  geom_point(aes(m, u3), color = "yellow", size = 4) +
  geom_line(aes(m, u3), color = "yellow", lwd = 3) +
  geom_line(aes(m, u3), lwd = .2) +
  # 雑多な設定
  theme_minimal() +
  labs(x = "number of success", y = "utility")






