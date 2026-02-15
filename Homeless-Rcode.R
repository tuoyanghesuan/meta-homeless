
library(meta)


data <- read.csv("SUDs.csv",header = T)


shapiro.test(data$p) 
qqnorm(data$arcsin.size)


data <- transform(data, p = event/n, log = log(event/n), 
                  logit = log((event/n)/(1 - event/n)), 
                  arcsin.size = asin(sqrt(event/(n + 1))), 
                  darcsin = 0.5 * (asin(sqrt(event/(n + 1)) )) + 
                    asin((sqrt(event + 1)/(n + 1))))

shapiro.test(data$arcsin.size)


library(dplyr)
data <- data %>% arrange(desc(data$p))
meta1 <- metaprop(event,n,data=data,
                  studlab = paste (data$study,data$year,seq=" "), 
                  sm="PFT",
                  )


forest(meta1, 
       digits=3, digits.I2 = 2, digits.tau2 = 3, digits.pval =5
      )
meta1


forest(meta1, layout = "Revman",
       
       leftcols = c("studlab", "event", "n","effect.ci"),
       
       rightcols = F,
       
       digits=3, digits.I2 = 2, digits.tau2 = 3, digits.pval =3,
       
       family="sans",
       
       fontsize=9.5,lwd=2,
       
       col.diamond.fixed="lightslategray",
       
       col.diamond.lines.fixed="lightslategray", 
       
       col.diamond.random="#669999",
       
       col.diamond.lines.random="maroon",
       
       col.square="#336699",
       
       col.square.lines = "white",
       
       col.study="#006666",
       
       col.inside = "#FF0033",
       
       col.label.right = "#FF0033", 
       
       col.random = "#660099",
       
       col.subgroup = "#003399",
       
       lty.fixed=4,
       
       plotwidth="8cm",
       
       colgap.forest.left="0.5cm",
       
       colgap.forest.right="1cm",
       
       just.forest="right",
       
       colgap.left="0.5cm",
       
       colgap.right="0.5cm",
       
       common = F, 
       
       xlab = "disease", 
       
       xlim = c(0,1),
       
       ref = 0,
       
       squaresize = 0.6, w.random = T,
       
       smlab = "Proportion",
       
       sep.subgroup =":"
      
       )



#leave-one-study-out sensitivity,
meta2 <- metainf(meta1,pooled = "random")
meta2
forest(meta2, rightcols = "effect.ci", "I2")

funnel(meta1, yaxis = "size")

metabias(meta1, k.min = 5)

tf <- trimfill(meta1,comb.random = TRUE)

summary(tf)

funnel(tf)

metareg(meta1,~bias + n + year + MA + SR) 

meta3 <- metaprop(event,n,data=data,
                  studlab = paste (data$study,data$year,seq=" "), 
                  sm="PFT",
                  subgroup = Instrument,
                  method.tau = "SJ",
                  hakn = T,
                  overall = F)
meta3
forest(meta3,
       digits=3, digits.I2 = 2, digits.tau2 = 3, digits.pval =3, common = F,
       leftcols = c("studlab", "event", "n","effect.ci"),
       rightcols = F
       )



