library(ggplot2)
library(dplyr)
library(plotly)

algs <- c(
  'nelder-mead'
  ,'annealing'
  ,'pso'
  # ,'bayesian'
  )

dff <- data.frame()
for(n in seq(9)) {
  for(alg in algs) {
    df <- read.csv(paste0('output/', alg, '_', n, '_trace.csv'))
    df$alg <- alg
    df$n_matrices <- as.character(n)
    dff <- bind_rows(dff, df)
    
    df2 <- df
    df2$alg <- paste0(alg, '_smoothed')
    df2$error <- cummin(df2$error)
    dff <- bind_rows(dff, df2)
    cat(paste0(alg, ' n=', n, '\n'))
  }
  
  df <- read.csv(paste0('output/gpu/bayesian_', n, '_trace.csv'))
  df$alg <- 'bayesian'
  df$n_matrices <- as.character(n)
  df$type <- 'gpu'
  dff <- bind_rows(dff, df)
  
  df2 <- df
  df2$alg <- paste0('bayesian_smoothed')
  df2$error <- cummin(df2$error)
  dff <- bind_rows(dff, df2)
  
  df <- read.csv(paste0('output/cpu/bayesian_', n, '_trace.csv'))
  df$alg <- 'bayesian'
  df$n_matrices <- as.character(n)
  df$type <- 'cpu'
  dff <- bind_rows(dff, df)
  
  df2 <- df
  df2$alg <- paste0('bayesian_smoothed')
  df2$error <- cummin(df2$error)
  dff <- bind_rows(dff, df2)
  
  cat(paste0('bayesian n=', n, '\n'))
}

# Error plots by algorithm
for(alg in algs) {
  png(paste0('output/alg_', alg, '.png'), width=1000, height=600)
  plt <- ggplot(dff[dff$alg==alg,], aes(x=index, y=error, color=n_matrices)) + 
    geom_line() +
    scale_x_continuous(labels=function(x) format(x, scientific=F, big.mark=',')) +
    scale_y_continuous(limits=c(0,10)) +
    xlab('Evaluations') +
    ylab('Error') +
    ggtitle(paste0(alg, ' by n_matrices'))
  # ggsave(plot=plt, filename=paste0('output/alg_', alg, '.png'), width = 10, height=10)
  print(plt)
  dev.off()
  
  png(paste0('output/alg_', alg, '_smoothed.png'), width=1000, height=600)
  plt <- ggplot(dff[dff$alg==paste0(alg, '_smoothed'),], aes(x=index, y=error, color=n_matrices)) + 
    geom_line() +
    scale_x_continuous(labels=function(x) format(x, scientific=F, big.mark=',')) +
    scale_y_continuous(limits=c(0,10)) +
    xlab('Evaluations') +
    ylab('Error') +
    ggtitle(paste0(alg, '(smoothed) by n_matrices'))
  # ggsave(plot=plt, filename=paste0('output/alg_', alg, '_smoothed.png'), width = 10, height=10)
  print(plt)
  dev.off()
}

# Plots by n_matrices
for(n in seq(9)) {
  png(paste0('output/n_', n, '.png'), width=1000, height=600)
  df.p <- dff[dff$n_matrices==n,]
  df.p[is.na(df.p$type), 'type'] <- 'cpu'
  plt <- ggplot(df.p[!endsWith(df.p$alg, '_smoothed') & df.p$type=='cpu',], aes(x=index, y=error, color=alg)) + 
    geom_line() +
    scale_x_continuous(labels=function(x) format(x, scientific=F, big.mark=',')) +
    scale_y_continuous(limits=c(0,max(df.p$error))) +
    xlab('Evaluation') +
    ylab('Error') +
    ggtitle(paste0('n_matrices = ', n, ' by method'))
  # ggsave(plot=plt, filename=paste0('output/alg_', alg, '.png'), width = 10, height=10)
  print(plt)
  dev.off()
  
  pdf(paste0('../../output/n_', n, '_smoothed.pdf'), width=9, height=4)
  plt <- ggplot(df.p[endsWith(df.p$alg, '_smoothed') & df.p$type == 'cpu',], aes(x=log10(index+1), y=error, color=alg)) + 
    geom_line(size=1) +
    scale_x_continuous(limits=c(0,4),labels=function(x) format(x, scientific=F, big.mark=',')) +
    scale_y_continuous(limits=c(0.5,1.5)) +
    scale_color_manual(name='Method', 
                         breaks=c('nelder-mead_smoothed', 
                                  'annealing_smoothed', 
                                  'pso_smoothed',
                                  'bayesian_smoothed'), 
                         labels=c('Nelder-Mead',
                                  'Simulated Annealing',
                                  'PSO', 
                                  'BO-SE'),
                         values=c('#a3a500',
                                  '#e76bf3',
                                  '#00bf7d',
                                  '#f8766d')) +
    xlab('Number of evaluations (log scale)') +
    ylab('Cumulative minimum error') +
    theme_minimal() +
    theme(legend.position = 'bottom')
    # ggtitle(paste0('n_matrices = ', n, ' by method (smoothed)'))
  # ggsave(plot=plt, filename=paste0('output/alg_', alg, '.png'), width = 10, height=10)
  print(plt)
  dev.off()
}

### Bayesian outputs

dfb <- dff[dff$alg=='bayesian',]
dfb$label <- paste0(dfb$n_matrices, ' (', dfb$type, ')')

dfbs <- dfb.cpu %>% group_by(n_matrices) %>% summarise(avg.time=mean(time), median.time=median(time))
dfbs$n_matrices <- as.numeric(dfbs$n_matrices)

dfb.cpu <- dfb[dfb$type=='cpu',]
dfbs.cpu <- dfb.cpu %>% group_by(n_matrices) %>% summarise(avg.time=mean(time), median.time=median(time))
dfbs.cpu$n_matrices <- as.numeric(dfbs.cpu$n_matrices)

dfb.gpu <- dfb[dfb$type=='gpu',]
dfbs.gpu <- dfb.gpu %>% group_by(n_matrices) %>% summarise(avg.time=mean(time), median.time=median(time))
dfbs.gpu$n_matrices <- as.numeric(dfbs.gpu$n_matrices)

# m.mean.exp <- lm(log(avg.time)~n_matrices, dfbs)
# dfbs$expected.avg.time.exp <- exp(predict(m.mean.exp, list(n_matrices=dfbs$n_matrices)))
# 
# m.median <- lm(median.time~exp(n_matrices), dfbs)
# dfbs$expected.median.time <- exp(predict(m.median, list(n_matrices=dfbs$n_matrices)))
# 
# m.mean.squared <- lm(avg.time^(1/2)~n_matrices, dfbs)
# dfbs$expected.avg.time.squared <- exp(predict(m.mean.squared, list(n_matrices=dfbs$n_matrices)))
# 
# m.mean.cubic <- lm(avg.time^(1/3)~n_matrices, dfbs)
# dfbs$expected.avg.time.cubic <- exp(predict(m.mean.cubic, list(n_matrices=dfbs$n_matrices)))

# BO time plots by n_matrices
plt.time <- ggplot(dfb, aes(x=index, y=time, color=n_matrices, linetype=type)) + 
  geom_line() +
  xlab('Iteration') +
  ylab('Time (s)') +
  ggtitle('BO iteration time by n_matrices')
ggsave('output/figures/bo_time.png', plt.time, width=10, height=8)
print(plt.time)

# BO time plot for many evaluations (CPU + 1 matrix)
df.long <- read.csv('output/bayesian_1_trace_2000iter_cpu.csv')
df.long$index <- 1:nrow(df.long)
df.long$col <- 'red'
plt.long.cpu <- ggplot(df.long, aes(x=index, y=time)) + 
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), mapping=aes(color='red')) +
  geom_smooth(method='lm', formula=y~poly(x,3), mapping=aes(color='blue')) +
  xlab('Iteration') +
  ylab('Time (s)') +
  scale_color_manual(values=c('red', 'blue'), breaks=c('red', 'blue'), labels=c('Quadratic', 'Cubic'), name='Fit') +
  ggtitle('BO iteration time (CPU, 1 matrix)')
ggsave('output/figures/bo_time_cpu_long.png', plt.long.cpu, width=10, height=8)
print(plt.long.cpu)

# BO time plot for many evaluations (GPU + 1 matrix)
df.long <- read.csv('output/bayesian_1_trace_2000iter_gpu.csv')
df.long$index <- 1:nrow(df.long)
plt.long.gpu <- ggplot(df.long, aes(x=index, y=time)) + 
  geom_point() +
  geom_smooth(method='lm', formula=y~poly(x,2), mapping=aes(color='red')) +
  geom_smooth(method='lm', formula=y~poly(x,3), mapping=aes(color='blue')) +
  xlab('Iteration') +
  ylab('Time (s)') +
  scale_color_manual(values=c('red', 'blue'), breaks=c('red', 'blue'), labels=c('Quadratic', 'Cubic'), name='Fit') +
  ggtitle('BO iteration time (GPU, 1 matrix)')
ggsave('output/figures/bo_time_gpu_long.png', plt.long.gpu, width=10, height=8)
print(plt.long.gpu)

# BO time boxplots (CPU & GPU)
plt.box <- ggplot(dfb, aes(x=n_matrices, y=time, fill=type)) + 
  geom_boxplot() +
  xlab('Number of matrices') +
  ylab('Time (s)') +
  ggtitle('BO (CPU) iteration time boxplot by n_matrices')
ggsave('output/figures/bo_boxplot.png', plt.box, width=10, height=8)
print(plt.box)


# BO time boxplots (CPU)
plt.box.cpu <- ggplot(dfbs.cpu, aes(x=n_matrices, y=avg.time)) + 
  geom_boxplot(data=dfb.cpu, mapping = aes(x=n_matrices, y=time)) +
  geom_point(aes(x=n_matrices, y=avg.time), color='orange', pch=18, size=4) +
  geom_smooth(method='lm', formula=y~exp(x), mapping=aes(color='red')) +
  geom_smooth(method='lm', formula=y~poly(x,2), mapping=aes(color='blue')) +
  geom_smooth(method='lm', formula=y~poly(x,3), mapping=aes(color='green')) +
  xlab('Number of matrices') +
  ylab('Time (s)') +
  scale_color_manual(name='Fit', 
                     values=c('red', 'blue', 'green'), 
                     labels=c('Quadratic', 'Cubic', 'Exponential')) +
  ggtitle('BO (CPU) iteration time boxplot by n_matrices')
ggsave('output/figures/bo_boxplot_cpu.png', plt.box.cpu, width=10, height=8)
print(plt.box.cpu)


# BO time boxplots (GPU)
plt.box.gpu <- ggplot(dfbs.gpu, aes(x=n_matrices, y=avg.time)) + 
  geom_boxplot(data=dfb.gpu, mapping = aes(x=n_matrices, y=time)) +
  geom_point(aes(x=n_matrices, y=avg.time), color='orange', pch=18, size=4) +
  geom_smooth(method='lm', formula=y~exp(x), mapping=aes(color='red')) +
  geom_smooth(method='lm', formula=y~poly(x,2), mapping=aes(color='blue')) +
  geom_smooth(method='lm', formula=y~poly(x,3), mapping=aes(color='green')) +
  xlab('Number of matrices') +
  ylab('Time (s)') +
  scale_color_manual(name='Fit', 
                     values=c('red', 'blue', 'green'), 
                     labels=c('Quadratic', 'Cubic', 'Exponential')) +
  ggtitle('BO (GPU Quadro P620) iteration time boxplot by n_matrices')
ggsave('output/figures/bo_boxplot_gpu.png', plt.box.gpu, width=10, height=8)
print(plt.box.gpu)




# Cumulative errors on 1 matrix
df.nm <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_nelder-mead_1mat.txt', col.names = 'x')
df.nm$method <- 'nelder-mead'
df.nm$iter <- 1:nrow(df.nm)
df.nm$xm <- cummin(df.nm$x)

df.annealing <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_annealing_1mat.txt', col.names = 'x')
df.annealing$method <- 'annealing'
df.annealing$iter <- 1:nrow(df.annealing)
df.annealing$xm <- cummin(df.annealing$x)

df.pso <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_pso_1mat.txt', col.names = 'x')
df.pso$method <- 'pso'
df.pso$iter <- 1:nrow(df.pso)
df.pso$xm <- cummin(df.pso$x)

df.bo <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_bayesian_1mat.txt', col.names = 'x')
df.bo$method <- 'bo'
df.bo$iter <- 1:nrow(df.bo)
df.bo$xm <- cummin(df.bo$x)

df.bo.oak <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_bayesian_oak_1mat.txt', col.names = 'x')
df.bo.oak$method <- 'bo.oak.20'
df.bo.oak$iter <- 1:nrow(df.bo.oak)
df.bo.oak$xm <- cummin(df.bo.oak$x)

dff <- rbind(df.nm, df.annealing)
dff <- rbind(dff, df.pso)
dff <- rbind(dff, df.bo)
dff <- rbind(dff, df.bo.oak)

breaks <- c('initial', 'nelder-mead', 'annealing', 'pso', 'bo', 'bo.oak.20', 'bo.oak.25')
method.names <- c('initial', 'Nelder-Mead', 'Simulated Annealing', 'PSO', 'BO-SE', 'BO-OAK', 'BO-OAK (25 iterations)')
colors <- c('black', '#619cff', '#00ba38', '#f752df', '#f8766d', '#a8261d', 'yellow')

ggplot(dff, aes(x=log10(iter), y=xm, color=method)) + 
  geom_line() +
  coord_cartesian(ylim=c(0,3)) +
  theme_minimal() +
  xlab('Number of evaluations (log scale)') +
  ylab('Cumulative minimum error') +
  scale_color_manual(breaks=breaks, values=colors, labels=method.names, name='Method') +
  theme(plot.title=element_text(hjust=.5), legend.position = 'bottom')



# Comparison BO methods for 9 matrices

df.se <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_bo_se_delay0.txt')
df.se$error.m <- cummin(df.se$error)
df.se$method <- 'bo'
df.oak <- read.csv('../../output/stepwise_calibration/conventional_unconstrained_bo_oak_delay0.txt')
df.oak$error.m <- cummin(df.oak$error)
df.oak$method <- 'bo.oak.20'
df.both.bo <- rbind(df.se, df.oak)
ggplot(df.both.bo, aes(x=iter,y=error,color=method)) +
  geom_line() +
  coord_cartesian(ylim=c(0,6000)) +
  theme_minimal() +
  xlab('Number of evaluations') +
  ylab('Cumulative minimum error') +
  scale_color_manual(breaks=breaks, values=colors, labels=method.names, name='Method') +
  theme(plot.title=element_text(hjust=.5), legend.position = 'bottom')


# Table of conventional vs stepwise by method, 9 age groups
