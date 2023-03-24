
path <- "data/CEATTLE/"

biomass <- rbind(cbind.data.frame(read.csv(paste0(path, "2020/ceattle_intrasp_biomass.csv")),
                                  assessment = "2020"),
                 cbind.data.frame(read.csv(paste0(path, "2022/ceattle_intrasp_biomass.csv")),
                                  assessment = "2022"))

# Add bounds for error & set 0 as minimum for plotting
biomass$min <- biomass$value - (2 * biomass$error)
biomass$min[biomass$min < 0] <- 0
biomass$max <- biomass$value + (2 * biomass$error)

popdy_plot <- ggplot(biomass, aes(x=year, y=value, color = assessment, fill = assessment)) +
    geom_line() +
    geom_ribbon(aes(ymin=min, ymax=max), alpha = 0.2, color = NA) + 
    scale_color_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) +  
    scale_fill_viridis(discrete = TRUE, direction = -1, begin = 0.1, end = 0.9) + 
    # geom_vline(xintercept = , linetype = 2, colour = "gray") +  # Add line at end of hindcast
    ylim(0, NA) +
    ylab(" ") +
    labs(color = "Assessment") +
    facet_wrap(~type, ncol = 2, scales = "free_y") 
popdy_plot

ggsave(filename="plots/CEATTLE/cannibalism/compare_assessments.png", popdy_plot, width=220, height=70, units="mm", dpi=300)
