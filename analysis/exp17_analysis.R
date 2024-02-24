
library(ggplot2)
library(viridis)

df <- read.csv("/Users/omarahmed/Downloads/work_dir/docarray_exps/exp_17/output.csv",
               header=FALSE,
               col.names=c("optimization", "error", "length", "newtime", "oldtime"))
df["speedup"] <- round(df["oldtime"]/df["newtime"], 2)
df["optimization"] <- factor(df$optimization, levels=c("optimized", "ftab", "minimizer", "minimizer+ftab"))

facet_names <- list(
  'optimized'="optimized",
  'ftab'="optimized + ftab",
  'minimizer'="optimized + minimizer",
  'minimizer+ftab'="optimized + minimizer + ftab"
)

facet_labeller <- function(variable,value){
  return(facet_names[value])
}

plot <- ggplot(df[df$length < 1500,]) +
        geom_point(aes(x=error, y=length, size=speedup, color=speedup)) +
        geom_text(aes(label=speedup, x=error, y=length+30), hjust=0.5, vjust=0) +
        facet_wrap(~optimization, labeller=facet_labeller)+
        theme_bw() +
        theme(axis.title.x=element_text(size=13),
              axis.title.y=element_text(size=13),
              strip.text.x=element_text(size=13),
              legend.title=element_text(size=13),
              axis.text.x=element_text(size=10),
              axis.text.y=element_text(size=10)) +
        labs(x="Error Rate", y="Read Length") +
        #scale_color_continuous(name="Speed-Up") +
        scale_size_continuous(name="Speed-Up") +
        scale_color_viridis(name="Speed-Up", option = "D") +
        scale_x_continuous(limit=c(0.01,0.10), breaks=seq(0.01, 0.10, 0.01)) +
        scale_y_continuous(limit=c(100, 1050), breaks=seq(0, 1000, 100))
plot

ggsave("/Users/omarahmed/Downloads/work_dir/docarray_exps/exp_17/speedup.png", 
       device="png",
       dpi=800,
       width=10,
       height=10,
       units="in")
