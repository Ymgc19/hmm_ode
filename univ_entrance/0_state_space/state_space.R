# これは，状態空間モデルの時変係数モデルを用いて，大学進学率を扱うファイル
# データの読み込みと前処理
source("univ_entrance/data_preprocessing.R")

# ===== MCMCのためのデータ ===== #
data = list(
  T = nrow(df),
  delta_t = 1,
  y = df$enter_univ_rate,
  x1 = df$salary_diff, # 初任給差
  x2 = df$jobs_to_applicants_ratio, # 有効求人倍率
  x3 = df$unemployment_rate # 失業率
)

data


# ========== モデルを読み込む ========== #
model <- cmdstan_model("univ_entrance/0_state_space/state_space.stan")

# ========== サンプリング ========== #
fit <- model$sample(
  data = data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 5000,
  refresh = 1000
)

# 結果を見る
options(scipen = 4)
fit$cmdstan_summary()





# ========== 係数の確認 ========== #
# 時変係数
b0 <- fit$draws(variables = c("b0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b1 <- fit$draws(variables = c("b1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b2 <- fit$draws(variables = c("b2"), format = "draws_matrix") %>% colMeans() %>% as.vector()
b3 <- fit$draws(variables = c("b3"), format = "draws_matrix") %>% colMeans() %>% as.vector()
# 分散
t0 <- fit$draws(variables = c("t0"), format = "draws_matrix") %>% colMeans() %>% as.vector()
t1 <- fit$draws(variables = c("t1"), format = "draws_matrix") %>% colMeans() %>% as.vector()
t2 <- fit$draws(variables = c("t2"), format = "draws_matrix") %>% colMeans() %>% as.vector()
t3 <- fit$draws(variables = c("t3"), format = "draws_matrix") %>% colMeans() %>% as.vector()
# 観測誤差
v <- fit$draws(variables = c("v"), format = "draws_matrix") %>% colMeans() %>% as.vector()
# 状態
alpha <- fit$draws(variables = c("alpha"), format = "draws_matrix") %>% colMeans() %>% as.vector()


# ========== 実測値とトレンド（alpha）の比較 ========== #

ggplot() +
  # 実測値
  geom_line(
    aes(df$year, df$enter_univ_rate), colour = "tomato"
  ) +
  geom_point(
    aes(df$year, df$enter_univ_rate), colour = "tomato"
  ) +
  # 状態空間のトレンド
  geom_line(
    aes(df$year, alpha), colour = "royalblue"
  ) +
  geom_point(
    aes(df$year, alpha), colour = "royalblue"
  ) +
  theme_minimal()






# ==================== 係数の詳細な可視化 ==================== #
# 切片
p0 <- ggplot() +
  geom_ribbon(
    aes(x = df$year, 
         ymin = b0 - 2*t0,
         ymax = b0 + 2*t0),
    fill = "cyan", alpha = .5
  ) +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b0 - t0,
        ymax = b0 + t0),
    fill = "royalblue", alpha = .5
  ) +
  geom_point(
    aes(df$year, b0)
  ) +
  geom_line(
    aes(df$year, b0)
  ) +
  theme_minimal() +
  labs(x = "year", y = "b0, intercept")
# 初任給差の係数
p1 <- ggplot() +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b1 - 2*t1,
        ymax = b1 + 2*t1),
    fill = "cyan", alpha = .5
  ) +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b1 - t1,
        ymax = b1 + t1),
    fill = "royalblue", alpha = .5
  ) +
  geom_point(
    aes(df$year, b1)
  ) +
  geom_line(
    aes(df$year, b1)
  ) +
  theme_minimal() +
  labs(x = "year", y = "b1, differences in starting salary")
# 有効求人倍率の係数
p2 <- ggplot() +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b2 - 2*t2,
        ymax = b2 + 2*t2),
    fill = "cyan", alpha = .5
  ) +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b2 - t2,
        ymax = b2 + t2),
    fill = "royalblue", alpha = .5
  ) +
  geom_point(
    aes(df$year, b2)
  ) +
  geom_line(
    aes(df$year, b2)
  ) +
  theme_minimal() +
  labs(x = "year", y = "b2, jobs to applicants ratio")
# 完全失業率の係数
p3 <- ggplot() +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b3 - 2*t3,
        ymax = b3 + 2*t3),
    fill = "cyan", alpha = .5
  ) +
  geom_ribbon(
    aes(x = df$year, 
        ymin = b3 - t3,
        ymax = b3 + t3),
    fill = "royalblue", alpha = .5
  ) +
  geom_point(
    aes(df$year, b3)
  ) +
  geom_line(
    aes(df$year, b3)
  ) +
  theme_minimal() +
  labs(x = "year", y = "b3, unemployment rate")

p0 + p1 + p2 + p3
