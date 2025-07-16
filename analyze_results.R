library(here)
library(data.table)
library(ggplot2)

df <- lapply(list.files("outputs_unif/", full.names = TRUE), function(file){
  dt <- readRDS(file)
  return(dt$est)
})

dt <- rbindlist(df)

ggplot(dt, aes(X1, X2)) + 
  geom_segment(aes(col = method)) + 
  coord_equal()

ggplot(dt, aes(x = 0, y = 0, xend = X1, yend = X2)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), aes(col = method), alpha = 0.6) +
  coord_equal() +
  theme_minimal() +
  labs(title = "Estimated Vectors and True Vector",
       x = "X", y = "Y") +
  theme(legend.title = element_blank())
