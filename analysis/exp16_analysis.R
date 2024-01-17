
# Name: exp16_analysis.R
# Description: Analyze the speed-up achieved when using the optimized
#              query algorithm and the ftab querying strategy
# Author: Omar Ahmed

library(ggplot2)
library(viridis)

df <- read.csv("/Users/omarahmed/Downloads/output.csv",
               header=FALSE,
               col.names=c("optimization", "error", "length", "newtime", "oldtime"))
df["speedup"] <- df["oldtime"]/df["newtime"]


facet_names <- list(
  'ftab'="Ftab Query + Optimized Algorithm",
  'optimized'="Optimized Algorithm"
)

facet_labeller <- function(variable,value){
  return(facet_names[value])
}

plot <- ggplot(df) +
        geom_point(aes(x=error, y=length, size=speedup, color=speedup)) +
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
        scale_x_continuous(limit=c(0.01,0.10), breaks=seq(0.01, 0.10, 0.01))
plot

ggsave("/Users/omarahmed/Downloads/speedup.png", 
       device="png",
       dpi=800,
       width=9,