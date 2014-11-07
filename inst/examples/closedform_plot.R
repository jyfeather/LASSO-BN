library("ggplot2")

df <- data.frame()
ggplot(df) + geom_point() + 
  geom_hline(yintercept=c(-1, 1)) +
  geom_vline(xintercept=c(-1, 1)) +
  annotate("text", x=-1.65, y=1.65, label="1\nb1 < 0\nb2 > 0") +
  annotate("text", x=0, y=1.65, label="2\nb1 = 0\nb2 > 0", colour="red") +
  annotate("text", x=1.65, y=1.65, label="3\nb1 > 0\nb2 > 0") +
  annotate("text", x=-1.65, y=0, label="4\nb1 < 0\nb2 = 0", colour="red") +
  annotate("text", x=0, y=0, label="5\nb1 = 0\nb2 = 0") +
  annotate("text", x=1.65, y=0, label="6\nb1 > 0\nb2 = 0", colour="red") +
  annotate("text", x=-1.65, y=-1.65, label="7\nb1 < 0\nb2 < 0") +
  annotate("text", x=0, y=-1.65, label="8\nb1 = 0\nb2 < 0", colour="red") +
  annotate("text", x=1.65, y=-1.65, label="9\nb1 > 0\nb2 < 0") +
  xlim(-2, 2) + xlab("zt1") +
  ylim(-2, 2) + ylab("zt2") +
  ggtitle("LASSO-BN Solution Region") + 
  theme_bw() + 
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) 

a = 0.25
ggplot(df) + geom_point() + 
  geom_segment(aes(x=a-1, xend=a+1, y=1, yend=1)) +
  geom_segment(aes(x=-a-1, xend=1-a, y=-1, yend=-1)) +
  geom_segment(aes(x=-a-1, xend=a-1, y=-1, yend=1)) +
  geom_segment(aes(x=1-a, xend=1+a, y=-1, yend=1)) +
  geom_vline(xintercept = c(a-1, 1-a, -1-a, a+1)) +
  geom_abline(intercept=a^2+a+1, slope=-a) +
  geom_abline(intercept=-a^2+a-1, slope=-a) +
  geom_abline(intercept=a^2-a+1, slope=-a) +
  geom_abline(intercept=-a^2-a-1, slope=-a) +
  annotate("text", x=0, y=0, label="5\nb1 = 0\nb2 = 0") +
  xlim(-2, 2) + xlab("zt1") +
  ylim(-2, 2) + ylab("zt2") +
  ggtitle("VS-MSPC Solution Region") +
  theme_bw() + 
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
  ) 