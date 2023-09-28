library(ggplot2)
library(tidyr)


# Time results:
data_file <- "/Users/omarahmed/Downloads/output_time.csv"
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
        scale_x_continuous(breaks=seq(0, 40, 5)) +
        scale_y_continuous(breaks=seq(0, 1750, 250)) +
        scale_color_discrete(name="Build Mode",
                             labels=c("heuristic"="Heuristic",
                                      "noheuristic"="No Heuristic",
                                      "twopasstotal"="Both Passes",
                                      "firstpass"="1st Pass",
                                      "secondpass"="2nd Pass"))
plot
ggsave("/Users/omarahmed/Downloads/output_time_plot.png", 
       plot=plot, 
       device="png", 
       dpi=400,
       width=7,
       height=5,
       units=c("in"))


# Memory results:
data_file <- "/Users/omarahmed/Downloads/output_mem.csv"
df <- read.csv(data_file, header=TRUE)

df['heurmaxrss'] <- df['heurmaxrss']/1048576
df['noheurmaxrss'] <- df['noheurmaxrss']/1048576
df['tpassmaxrss'] <- df['tpassmaxrss']/1048576

df['heurwriteio'] <- df['heurwriteio']/1073741824
df['heurreadio'] <- df['heurreadio']/1073741824

df['noheurwriteio'] <- df['noheurwriteio']/1073741824
df['noheurreadio'] <- df['noheurreadio']/1073741824

df['tpasswriteio'] <- df['tpasswriteio']/1073741824
df['tpassreadio'] <- df['tpassreadio']/1073741824

new_df <- pivot_longer(df, 
                       c("heurmaxrss", "heurwriteio", "heurreadio", 
                         "noheurmaxrss", "noheurwriteio", "noheurreadio",
                         "tpassmaxrss", "tpasswriteio", "tpassreadio"),
                       names_to="type",
                       values_to="value")
new_df["metric"] <- ""

# Go through and split it 
for (row in 1:nrow(new_df)) {
  name <- new_df[row, "type"]
  if (name == "heurmaxrss" || name == "noheurmaxrss" || name == "tpassmaxrss") {
    new_df[row, 'metric'] <- "maxrss"
  }
  if (name == "heurwriteio" || name == "noheurwriteio" || name == "tpasswriteio") {
    new_df[row, 'metric'] <- "write"
  }
  if (name == "heurreadio" || name == "noheurreadio" || name == "tpassreadio") {
    new_df[row, 'metric'] <- "read"
  }
  
  if (name == "heurmaxrss" || name == "heurreadio" || name == "heurwriteio") {
    new_df[row, 'type'] <- "heuristic"
  }
  if (name == "noheurmaxrss" || name == "noheurreadio" || name == "noheurwriteio") {
    new_df[row, 'type'] <- "no_heuristic"
  }
  if (name == "tpassmaxrss" || name == "tpassreadio" || name == "tpasswriteio") {
    new_df[row, 'type'] <- "two_pass"
  }
}

facet_names <- list(
  'maxrss'="Maximum RSS",
  'read'="Reads",
  'write'="Writes"
)

facet_labeller <- function(variable,value){
  return(facet_names[value])
}

pd <- position_dodge(3)
plot <- ggplot(new_df, aes(x=num)) + 
        geom_point(aes(group=type, y=value, color=type), position=pd) +
        geom_line(aes(group=type, y=value, color=type), position=pd) +
        facet_wrap(~metric, labeller=facet_labeller) +
        theme_bw() +
        theme(axis.title.x=element_text(size=14),
              axis.text.x=element_text(size=12),
              axis.title.y=element_text(size=14),
              axis.text.y=element_text(size=12)) +
        labs(x="Number of Documents",
             y="Memory (GB)") +
        #scale_x_continuous(breaks=seq(0, 40, 5)) +
        scale_y_continuous(breaks=seq(0, 5, 0.5)) +
        scale_color_discrete(name="Build",
                             labels=c("heuristic"="Heuristic",
                                      "no_heuristic"="No Heuristic",
                                      "two_pass"="Two Passes"))

plot
ggsave("/Users/omarahmed/Downloads/output_mem_plot.png", 
       plot=plot, 
       device="png", 
       dpi=400,
       width=9,
       height=3,
       units=c("in"))