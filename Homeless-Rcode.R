#加载Meta包
library(meta)

#读取文件
data <- read.csv("SUDs.csv",header = T)

#检验正态性(p值大于0.05表示正态分布)
shapiro.test(data$p) 
qqnorm(data$arcsin.size)

#转换
data <- transform(data, p = event/n, log = log(event/n), 
                  logit = log((event/n)/(1 - event/n)), 
                  arcsin.size = asin(sqrt(event/(n + 1))), 
                  darcsin = 0.5 * (asin(sqrt(event/(n + 1)) )) + 
                    asin((sqrt(event + 1)/(n + 1))))

shapiro.test(data$arcsin.size)

#Meta分析
library(dplyr)
data <- data %>% arrange(desc(data$p))
meta1 <- metaprop(event,n,data=data,
                  studlab = paste (data$study,data$year,seq=" "), 
                  sm="PFT",
                  )

#森林图
forest(meta1, 
       digits=3, digits.I2 = 2, digits.tau2 = 3, digits.pval =5
      )
meta1

#优化的森林图
forest(meta1, layout = "Revman",
       
       leftcols = c("studlab", "event", "n","effect.ci"),
       
       rightcols = F,
       
       digits=3, digits.I2 = 2, digits.tau2 = 3, digits.pval =3,
       
       family="sans",
       
       fontsize=9.5,lwd=2,
       
       col.diamond.fixed="lightslategray",#固定效应菱形颜色
       
       col.diamond.lines.fixed="lightslategray", #固定效应菱形线条颜色
       
       col.diamond.random="#669999",#随机效应菱形颜色
       
       col.diamond.lines.random="maroon",#随机效应菱形线条颜色
       
       col.square="#336699",# 反映权重的方块颜色
       
       col.square.lines = "white",# 方块框颜色
       
       col.study="#006666",# 横线置信区间颜色
       
       col.inside = "#FF0033",#如果置信限完全在方块内，个体研究结果和置信限的颜色
       
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
       
       common = F, #是否统计固定效应模型
       
       xlab = "disease", #x轴标签
       
       xlim = c(0,1),
       
       ref = 0,# 参考线的横坐标
       
       squaresize = 0.6, w.random = T,
       
       smlab = "Proportion",
       
       sep.subgroup =":"
      
       )



#leave-one-study-out sensitivity,敏感性分析
meta2 <- metainf(meta1,pooled = "random")
meta2
forest(meta2, rightcols = "effect.ci", "I2")


#漏斗图-定性
funnel(meta1, yaxis = "size")


#Egger's-定量（p值小于0.05表示有发表偏倚）
metabias(meta1, k.min = 5)


#剪补法
tf <- trimfill(meta1,comb.random = TRUE)

summary(tf)

funnel(tf)


#Meta回归
metareg(meta1,~bias + n + year + MA + SR) 

metareg(meta1,~bias)
metareg(meta1,~n)
metareg(meta1,~year)
metareg(meta1,~MA)
metareg(meta1,~SR)


#亚组分析
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

subgroup = Risk.of.bias
subgroup = Sample.Size
subgroup = Continent
subgroup = Degree.of.Development
subgroup = Publication.Year
subgroup = Mean.Age
subgroup = Sex.Ratio
subgroup = Instrument


