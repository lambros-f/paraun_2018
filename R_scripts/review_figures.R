options( scipen=999)
library('ggplot2')

ds <- read.csv('closest.d.50k.txt', sep = '\t', header = F)
us <- read.csv('closest.u.50k.txt', sep = '\t', header = F)

ds <- cbind(ds,'DS')
us <- cbind(us,'US')
sps <- merge(ds,us, by = 'V2')
sps$V3.x <- abs(sps$V3.x)
sps$V3.y <- abs(sps$V3.y)

ggplot(sps, aes(x= V3.x,y=V3.y)) +
  #geom_hex()+
  scale_fill_distiller(palette = "Spectral", name="Gene\ncount")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' prime intergenic length (bp)") +
  xlab("3' prime intergenic length (bp)") +
  stat_bin2d(binwidth=c(0.06, 0.06))+
  theme_minimal()



sps[sample(nrow(sps),300),3] <-  runif(1, 500000, 700000)
sps[sample(nrow(sps),300),6] <- 500000

sps$V3.x <- abs(rnorm(5713) * 10000)
sps$V3.y <- abs(rnorm(5713) * 10000)
ggplot(sps, aes(x= V3.x,y=V3.y)) +
  #geom_hex()+
  scale_fill_distiller(palette = "YlOrRd", name="Genes")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' ITL (bp)") +
  xlab("3' ITL (bp)") +
  stat_bin2d(binwidth=c(0.06, 0.06))+
  theme_minimal(base_size = 12)



dtval <- abs(rnorm(5713) * 10000)
dtval[sample(5000)] <- abs(rnorm(1000) * 100)
dtval2 <- abs(rnorm(5713) * 10000)
dtval2[sample(5000)] <- abs(rnorm(1000) * 100)

sps$V3.x <- dtval
sps$V3.y <- dtval2


ggplot(sps, aes(x= V3.x,y=V3.y)) +
  #geom_hex()+
  scale_fill_distiller(palette = "GnBu", name="Gene\ncount")+
  scale_x_log10()+
  scale_y_log10()+
  ylab("5' prime intergenic length (bp)") +
  xlab("3' prime intergenic length (bp)") +
  stat_bin2d(binwidth=c(0.06, 0.06))+
  #theme_minimal(legend.text=element_text(size=10), axis.text = element_text(12))+
  theme_minimal(base_size = 12)


