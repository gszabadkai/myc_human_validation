library(tidyverse)

mitoPPS <- readRDS(here::here("Data/Guo_et_al", "ProcessedData", "mitoPPS_bins.rds"))
meta_binned <- read.csv(here::here("Data/Guo_et_al", "ProcessedData", "meta_binned.csv"))
mitoPPS %>%
  group_by(Pathway1) %>%
  inner_join(meta_binned, by = "bin") %>%
  filter(grepl("assembly|subunits", Pathway1)) %>%
  filter(!Pathway1 %in% c("OXPHOS subunits", "Mitochondrial ribosome assembly",
                          "OXPHOS assembly factors")) %>%
  mutate(Category = case_when(
    grepl("assembly" , Pathway1) ~ "Assembly factors",
    grepl("subunits" , Pathway1) ~ "Subunits"
  )) %>%
  separate(Pathway1, into = c("Complex", "rem1", "rem2"), sep = " ", remove = F) %>%
  
  ggplot(aes(x = bin, y = mitoPPS, color = Category)) +
  geom_point(aes(y = mitoPPS), alpha = 0.3, size = 0.1) + 
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs", k = 15), linewidth = 0.6) +
  scale_color_manual(values = c(
    "Assembly factors" = "#4A4C85",
    "Subunits" = "#B24A5C"
  )) +
  theme_classic() +
  facet_wrap(~Complex, scales = "free", ncol = 6) +
  theme(legend.position = "none",
        axis.title = element_text(size =7),
        axis.text = element_text(size = 6),
        strip.text = element_text(size = 7),
        axis.line = element_line(size = 0.3),
        strip.background = element_blank()) 
ggsave(here::here("Figures", "Figure6", "Figure6D_assembly_vs_subunits.png"),
       width = 8, height =1.6, dpi = 1200)
