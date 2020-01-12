library(extrafont)
#loadfonts(device = "win")
library(ggplot2)
library(cowplot)
###
TextSize = 13 #Presentation(s)
#TextSize = 10 #Presentation(s)
###
geom.text.size =  TextSize / ggplot2::.pt #geom_text uses "mm", ggplot theme uses point
geom.text.size.small = geom.text.size / 2
#theme.size = ggplot2::.pt * geom.text.size
theme.size = TextSize
theme.size.small = theme.size / 2
theme_kat <- theme_cowplot(line_size = 1) + theme(axis.ticks.length=unit(-0.2, "cm"), axis.text.x = element_text(margin=margin(8,0,0,0,"pt"), size = theme.size),
                                                  axis.title.y = element_text(size=theme.size),
                                                  axis.text.y = element_text(margin=margin(0,10,0,0,"pt"), size = theme.size),
                                                  axis.title.x = element_text(size=theme.size),
                                                  legend.text = element_text(size=theme.size), legend.title = element_text(size=theme.size), legend.text.align = 1,
                                                  strip.text = element_text(size=theme.size), line = element_line(size=unit(0.7, "cm")),
                                                  text = element_text(family = "Arial", size = TextSize))

theme_set(theme_kat)
