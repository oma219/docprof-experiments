library(ggplot2)
library(tidyr)

data_file <- "/Users/omarahmed/Downloads/output.csv"
df <- read.csv(data_file, header=TRUE)
new_df <- pivot_longer(df, 
                      c("heuristic", "noheuristic", "twopasstotal", 
                        "firstpass", "secondpass"),
                      names_to="build",
                      values_to="time")

new_df$build <- factor(new_df$build, levels=c("noheuristic", "heuristic", "twopasstotal", "firstpass", "secondpass"))
plot <- ggplot(new_df, aes(x=num)) + 
        geom_point(aes(group=build, y=time, color=build)) +
        geom_line(aes(group=build, y=time, color=build)) +
        theme_bw() +
        theme(axis.title.x=element_text(size=14),
              axis.text.x=element_text(size=12),
              axis.title.y=element_text(size=14),
              axis.text.y=element_text(size=12)) +
        labs(x="Number of Documents",
             y="Time (sec)",
             title="") +
        scale_color_discrete(name="Build Mode",
                             labels=c("heuristic"="Heuristic",
                                    "noheuristic"="No Heuristic",
                                    "twopasstotal"="Both Passes",
                                    "firstpass"="1st Pass",
                                    "secondpass"="2nd Pass"))
plot
ggsave("/Users/omarahmed/Downloads/output.png", plot=plot, device="png", dpi=500, units="in")