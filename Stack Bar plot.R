
# ------------------------- Stack Bar plot -----------------------
library(dplyr)
library(ggplot2)

# Create a metadata dataframe 

meta_df <- data.frame(
  cluster   = Idents(cluster1) %>% as.character(),
  condition = cluster1@meta.data$orig.ident %>% as.character()
)


# Count number of cells per cluster per condition, and calculate their proportion

counts_df <- meta_df %>%
  group_by(cluster, condition) %>%
  summarize(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

# Set the order of condition levels
counts_df <- counts_df %>%
  mutate(condition = factor(condition, levels = c("Control", "Experiment")))



# Set color
my_cluster_colors <- c(
  "0" = "#D8D8D8",
  "1" = "#FFA726",
  "2" = "#BFCFE7",
  "3" = "#26C6DA",
  "4" = "#E3DFFF",
  "5" = "#CFE0E8"
)


# 1
ggplot(counts_df, aes(x = cluster, y = n, fill = condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Cluster", y = "Cell Counts", fill = "Condition") +
  theme_classic()

# 2
ggplot(counts_df, aes(x = condition, y = prop, fill = cluster)) +
  geom_bar(stat = "identity",width = 0.6, position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(x = "Condition", y = "Proportion", fill = "cluster") +
  theme_minimal()


# 3
ggplot(counts_df, aes(x = condition, y = prop, fill = cluster)) +
  geom_bar(stat = "identity", width = 0.4) +
  geom_text(aes(label = percent(prop)),
            position = position_stack(vjust = 0.5),
            size = 3) +
  scale_fill_manual(values = FAP_cluster_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Condition", y = "Proportion", fill = "Cluster") +
  theme_classic() +   
  theme(
    panel.background      = element_blank(),  
    panel.border          = element_blank(),  
    panel.grid.major      = element_blank(),  
    panel.grid.minor      = element_blank()   
  ) 







