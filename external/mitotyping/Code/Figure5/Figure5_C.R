library(tidyverse)
library(broom)

color_celltypes =  c(
  "unspecific" = "grey",
  "Inhibitory" = "#FF7F00",# 
  "Excitatory" = "#007200",
  "OPC" = "#A6CEE3", 
  "Oligodendrocytes" = "#1F78B4",
  "Astrocytes" = "#E15A57",
  "Microglia" ="#FFFF99",
  "VLMCs" =  "#6A3D9A",
  "Pericytes" = "#D0BEDA",
  "Endothelial" = "#B789DA",
  "Brain"= "#1DB100"
  
)

mitoPPS_MBM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "MBM_mitoPPS.csv")) %>%
  filter(!Voxel %in% "WM") %>%
  group_by(CellType, Pathway1) %>%
  summarize(mitoPPS = mean(mitoPPS, na.rm = TRUE), .groups = "drop") 

mitoPPS_RM <- read_csv(here::here("Data", "MBM_ROSMAP", "ProcessedData", "ROSMAP_mitoPPS.csv")) %>%
  separate(Group, into = c("CellType", "ID"), sep = "_", extra = "merge") %>%
  mutate(CellType = case_when(
        CellType == "Ast" ~ "Astrocytes",
        CellType == "End" ~ "Endothelial",
        CellType == "Exc" ~ "Excitatory",
        CellType == "Inh" ~ "Inhibitory",
        CellType == "Mic" ~ "Microglia",
        CellType == "Oli" ~ "Oligodendrocytes",
        CellType == "OPC" ~ "OPC",
        CellType == "Peri" ~ "Pericytes",
        TRUE ~ CellType)) 

RM_av <- mitoPPS_RM %>%
  group_by(CellType, Pathway1) %>%
  mutate(av = median(mitoPPS, na.rm = T)) %>%
  select(CellType, Pathway1, av) %>%
  unique()

combined <- RM_av %>%
  full_join(mitoPPS_MBM, by = c("CellType", "Pathway1")) %>%
  dplyr::rename(RM = av, MBM = mitoPPS)
cor.test(combined$RM, combined$MBM, method = "pearson")


test <- combined %>%
  group_by(CellType) %>%
  nest(-CellType) %>%
  mutate(cor=map(data,~cor.test(.x$RM, .x$MBM, method = "pearson"))) %>%
  mutate(tidied = map(cor, tidy)) %>%
  unnest(tidied, .drop = T) %>%
  ungroup() %>%
  mutate(padj = p.adjust(p.value, method = "BH"))
mean(test$estimate)
sd(test$estimate)

top <- combined %>%
  pivot_longer(cols = c("RM", "MBM"), names_to = "Dataset", values_to = "mitoPPS") %>%
  group_by(CellType, Dataset) %>%
  arrange(desc(mitoPPS)) %>%
  dplyr::slice(1)
combined %>%
  ggplot(aes(x = RM, y = MBM, color= CellType)) +
  geom_point(size = 0.5, alpha = 0.8) +
  theme_bw() +
  scale_color_manual(values = color_celltypes) +
  geom_smooth(method = "lm", linewidth = 0.3,
              color = "black", linetype ="dotted", se = F) +
  facet_wrap(~fct_relevel(CellType, c("Astrocytes", "Oligodendrocytes", 
                                      "Excitatory", "Inhibitory", "OPC", 
                                      "Microglia", "Endothelial", "Pericytes")), scales = "free", ncol = 4) +
  theme( legend.position = "none",
        axis.text =  element_text(size = 6),
        axis.title =  element_text(size = 7),
        strip.text = element_text(size = 6, face = "bold"))

ggsave(here::here("Figures", "Figure5",  "Figure_5C_Corr_MBM_RM.png"), width = 4.3, height = 2.2, unit = "in", dpi = 1200)
