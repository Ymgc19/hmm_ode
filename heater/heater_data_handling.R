library(tidyverse)

df <- read_csv("heater/heater.csv")
df %>% glimpse()

plot(df$year, df$heater_preveil)
plot(df$year, df$kerosene_price)




# ggplotによる描画
df %>% 
  ggplot() +
  # ファンヒーター普及率
#  geom_line(
#    aes(year, heater_preveil), 
#    color = "tomato", alpha = .35
#  ) +
#  geom_point(
#    aes(year, heater_preveil), 
#    color = "tomato", alpha = .35
#  ) +
  # 灯油価格指数
  geom_line(
    aes(year, kerosene_price)
  ) +
  geom_point(
    aes(year, kerosene_price)
  ) +
  theme_minimal() +
  labs(x = "Year", y = "Kerosene price index")
