library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(purrr)

######## C-SCNR ########
replicate = 20
all_results <- lapply(1:5, function(x) rep(list(NULL), replicate))

for (i in 1:5) {
  for (rep in 1:replicate) {
  #### C-SCNR
  test.cscnr = fread(paste0("prediction_interval/C-SCNR/PI_test_",i,"_rep_",rep,".txt"),header = T)
  all_results[[i]][[rep]] <- test.cscnr
  }
}

final_summarized_list <- map(all_results, function(subset_reps) {
  summarized_bounds <- subset_reps %>%
    map_dfr(~ as.data.frame(.x)[, c("IID", "lower", "upper")]) %>%
    group_by(IID) %>%
    summarise(
      mean_lower = mean(lower, na.rm = TRUE),
      mean_upper = mean(upper, na.rm = TRUE),
      .groups = "drop"
    )
  
  template <- as.data.frame(subset_reps[[1]]) %>%
    select(-any_of(c("Cluster_ID", "lower", "upper")))

  final_df <- template %>%
    left_join(summarized_bounds, by = "IID") %>%
    rename(lower = mean_lower, upper = mean_upper)
  
  return(final_df)
})

summarize <- function(df, method) {
  df %>%
    group_by(Base_Group) %>%
    summarise(
      Coverage = mean(HDL >= lower & HDL <= upper),
      Mean_Width = mean(upper - lower)
    ) %>%
    mutate(Method = method)
}

A <- final_summarized_list %>%
  map_dfr(~ summarize(.x, method = "C-SCNR"), .id = "Subset")

######## SCNR ########
results_scnr = rep(list(NULL), 5)
  
for (i in 1:5) {
  #### SCNR
  test.scnr = fread(paste0("prediction_interval/SCNR/PI_test_",i,".txt"),header = T)
  test.scnr$Base_Group = final_summarized_list[[i]]$Base_Group
  results_scnr[[i]] = test.scnr
}
  
B <- results_scnr %>%
  map_dfr(~ summarize(.x, method = "SCNR"), .id = "Subset")



######## Plot by Base Group ########
my_colors <- c(
  "C-SCNR" = "#EEAD0E",
  "SCNR" = "#E64B35"
)

my_theme <- theme_bw() + 
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    legend.title = element_blank(),
    legend.position = "top",
    strip.background = element_rect(fill = "grey95"),
  )

plot_data = rbind(A,B)
plot_data$Method = factor(plot_data$Method, levels = c("C-SCNR","SCNR"))

#### Prediction Coverage Rate 
p1 <- ggplot(plot_data, aes(x = Base_Group, y = Coverage, color = Method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.9), 
             size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "text",
               aes(label = after_stat(sprintf("%.3f", y))), 
               position = position_dodge(width = 0.85), vjust = -4, size = 2, fontface = "bold", show.legend = FALSE) +
  geom_hline(yintercept = 0.95, color = "orange", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0.9, 1), breaks = seq(0.9, 1, 0.025)) +
  labs(title = "HDL",
       x = "Base Group",
       y = "Coverage Rate") +
  my_theme 
p1


#### Mean Width
p2 <- ggplot(plot_data, aes(x = Base_Group, y = Mean_Width, color = Method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, jitter.height = 0, dodge.width = 0.9), 
             size = 2, alpha = 0.6) +
  stat_summary(fun = mean, geom = "text",
               aes(label = after_stat(sprintf("%.2f", y))), 
               position = position_dodge(width = 0.85), vjust = -3, size = 2, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values = my_colors) +
  coord_cartesian(ylim = c(min(plot_data$Mean_Width) - 0.5, 
                           max(plot_data$Mean_Width) + 0.5)) +
  labs(title = "HDL",
       x = "Base Group",
       y = "Mean Interval Width") +
  my_theme
p2


