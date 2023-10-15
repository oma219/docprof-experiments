library(ggplot2)
library(tidyr)
library(ggpubr)

#############################################
# Time Plot
#############################################
df <- read.csv("/Users/omarahmed/Downloads/output.csv", header=TRUE)
new_df <- pivot_longer(df, 
                       c("time", "pass1", "pass2"),
                       names_to="build",
                       values_to="time")

new_df$build <- factor(new_df$build, levels=c("time", "pass1", "pass2"))
group.colors <- c("time" = "#F8766D", 
                  "pass1" = "#A3A500", 
                  "pass2" ="#00BF7D")

plot1 <- ggplot(new_df, aes(x=num)) + 
        geom_point(aes(group=build, y=time, color=build)) +
        geom_line(aes(group=build, y=time, color=build)) +
        theme_bw() +
        theme(axis.title.x=element_text(size=14),
              axis.text.x=element_text(size=12),
              axis.title.y=element_text(size=14),
              axis.text.y=element_text(size=12)) +
        labs(x="Number of Human Genomes",
             y="Time (sec)",
             title="") +
        scale_x_continuous(breaks=seq(0, 10, 1)) +
        scale_y_continuous(breaks=seq(0, 40000, 5000)) +
        scale_color_manual(name="Build Mode",
                           labels=c("time"="Total Time",
                                    "pass1"="1st Pass",
                                    "pass2"="2nd Pass"), 
                           values=group.colors)
plot1

ggsave("/Users/omarahmed/Downloads/output_time_plot.png", 
       plot=plot1, 
       device="png", 
       dpi=400,
       width=7,
       height=5,
       units=c("in"))

#############################################
# Memory Plot
#############################################

df <- read.csv("/Users/omarahmed/Downloads/output.csv", header=TRUE)
new_df <- pivot_longer(df, 
                       c("maxrss", "tmpmem"),
                       names_to="metric",
                       values_to="memory")

new_df$build <- factor(new_df$metric, levels=c("maxrss", "tmpmem"))
group.colors <- c("maxrss" = "#F8766D", 
                  "tmpmem" = "#A3A500")

plot2 <- ggplot(new_df, aes(x=n)) + 
          geom_point(aes(group=metric, y=memory, color=metric)) +
          geom_line(aes(group=metric, y=memory, color=metric)) +
          theme_bw() +
          theme(axis.title.x=element_text(size=14),
                axis.text.x=element_text(size=12),
                axis.title.y=element_text(size=14),
                axis.text.y=element_text(size=12)) +
          labs(x="Size of Input File (GB)",
               y="Memory (GB)",
               title="") +
          # scale_x_continuous(breaks=seq(0, 10, 1)) +
          # scale_y_continuous(breaks=seq(0, 40000, 5000)) +
          scale_color_manual(name="Metric:",
                             labels=c("maxrss"="Max RSS",
                                      "tmpmem"="Temp. Mem Used"), 
                             values=group.colors)

plot2

ggsave("/Users/omarahmed/Downloads/output_mem_plot.png", 
       plot=plot2, 
       device="png", 
       dpi=400,
       width=7,
       height=5,
       units=c("in"))