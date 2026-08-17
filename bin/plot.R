#setwd("~/stingele_del/analysis")
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)

files <- c(
 "AAVS1_24_U.del",
 "AAVS1_31_U.del",
 "AAVS1_no_U.del",
 "HMCES_KO_24_U.del",
 "HMCES_KO_31_U.del",
 "HMCES_KO_no_U.del"
)

make_del_plot <- function(file,
                           show_x = TRUE,
                           show_y = TRUE,
                           show_legend = FALSE) {
  df <- read.csv(file, stringsAsFactors = FALSE) %>%
    mutate(
      focal_del_start = as.numeric(fdel_start),
      focal_del_end   = as.numeric(fdel_end),
      focal_del_size  = as.numeric(fdel_size)
    )

  plot_df <- df %>%
    filter(!is.na(focal_del_start), focal_del_size > 0,
           focal_del_end >= 70, focal_del_end <= 115) %>%
    mutate(
      size_group = cut(
        focal_del_size,
        breaks = c(2,12,22,32,42,100),
        labels = c("2–11","12–21","22–31","32–41", ">42"),
        include.lowest = TRUE,
        right = TRUE
      )
    ) %>%
    count(focal_del_end, size_group, .drop = FALSE) %>%
    mutate(n_thousands = n / 1000)

  p <- ggplot(plot_df, aes(x = focal_del_end, y = n_thousands, fill = size_group)) +
    geom_col(width = 1, colour = "white", linewidth = 0.15) +
    scale_x_continuous(limits = c(70,115), breaks = seq(70,115,5),
                        minor_breaks = NULL, expand = c(0,0)) +
    scale_y_continuous(limits = c(0,45), breaks = c(0,10,20,30,40),
                        expand = expansion(mult = c(0,0.02))) +
    scale_fill_brewer(palette = "Blues", direction = 1,
                       name = "Deletion\nsize (bp)", drop = FALSE) +
    labs(x = "", y = "Read counts (thousands)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title  = element_text(face = "bold"),
      axis.line   = element_line(colour = "black"),
      axis.ticks.length = unit(2, "mm"),
      plot.margin = margin(t = 8, r = 4, b = 8, l = 4)
    )

  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                               axis.ticks.x = element_blank(),
                               axis.title.x = element_blank())
  if (!show_y) p <- p + theme(axis.title.y = element_blank())
  if (!show_legend) p <- p + theme(legend.position = "none")

  p
}

plots <- list(
  make_del_plot(files[1], show_x = FALSE, show_y = FALSE),
  make_del_plot(files[2], show_x = FALSE, show_y = FALSE),
  make_del_plot(files[3], show_x = TRUE,  show_y = FALSE),
  make_del_plot(files[4], show_x = FALSE, show_y = FALSE),
  make_del_plot(files[5], show_x = FALSE, show_y = FALSE),
  make_del_plot(files[6], show_x = TRUE,  show_y = FALSE)
)
names(plots) <- files

plots[[1]] <- plots[[1]] + ggtitle("AAVS1")    + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
plots[[4]] <- plots[[4]] + ggtitle("HMCES KO") + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

final_plot <- (plots[[1]] | plots[[4]]) /
              (plots[[2]] | plots[[5]]) /
              (plots[[3]] | plots[[6]]) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right", plot.margin = margin(t = 8, r = 4, b = 4, l = 4))

y_label <- ggdraw() + draw_label("Read counts (thousands)", angle = 90, size = 14)
x_label <- ggdraw() + draw_label("Deletion end position (bp)", size = 14) +
  theme(plot.margin = unit(c(0.1,0.1,1.2,0.1), "cm"))

row1 <- plot_grid(y_label, final_plot, ncol = 2, rel_widths = c(0.04, 1))
final_plot_labeled <- plot_grid(row1, x_label, ncol = 1, rel_heights = c(1, 0.05))

ggsave(final_plot_labeled, 
       filename = "deletion_distribution.pdf",
       device = "pdf",
       height = 6, width = 5, units = "in")

