library(tidyverse)

df <- read_csv("heater/heater.csv")
df %>% glimpse()

plot(df$year, df$heater_preveil)
plot(df$year, df$kerosene_price)

