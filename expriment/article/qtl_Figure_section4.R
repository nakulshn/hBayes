#'
#' Run qtl_simulation_figs_section4.R first to get the MSE_sim.rds data
#'
#'
library(xtable)
library(beeswarm)
library(tidyverse)
library(latex2exp )
library(ggplot2)
library(tidyverse)
library(NPBayes)
save.fig=T
list.data <- readRDS("MSE_sim.rds")
list.data2 <- readRDS("MSE_sim2.rds")
list.data3 <- readRDS("MSE_sim3.rds")
list2env(list.data,globalenv())
list2env(list.data2,globalenv())
#MSE_beta_fb <- MSE_beta_fb.
Table <- rbind(c(mean(MSE_XB_lasso),mean(MSE_XB_ridge),mean(MSE_XB_fb),mean(MSE_XB_npbayes), mean(MSE_XB_HS),mean(MSE_XB_oracle)),
               c(mean(MSE_beta_lasso),mean(MSE_beta_ridge),mean(MSE_beta_fb),mean(MSE_beta_npbayes), mean(MSE_beta_HS), mean(MSE_beta_oracle)),
               c(mean(MSE_beta_smooth.lasso),mean(MSE_beta_smooth.ridge),mean(MSE_beta_smooth.fb),mean(MSE_beta_smooth.bayes), mean(MSE_beta_smooth_HS), mean(MSE_beta_smooth.oracle)))
rownames(Table) <- c("RMSE of $ {\\bf X} \\beta $","RMSE of $ \\beta $","RMSE of $ \\beta $ smooth")
colnames(Table) <- c("lasso",  "ridge", "mix","hBeta","horseshoe","oracle")

print.xtable(xtable(Table, digits=2),type="latex",sanitize.text.function = function(x) x)
data.RMSEbeta <- data.frame(Lasso = MSE_beta_lasso,
                            Ridge = MSE_beta_ridge,
                            hBayes =MSE_beta_npbayes,
                            Mix     =MSE_beta_fb,
                            Mr.Ash  = MSE_beta_ash,
                            EMVS    = MSE_beta_EMVS,
                            horseshoe = MSE_beta_HS,
                            oracle = MSE_beta_oracle)
data.RMSEbeta <- pivot_longer(data.RMSEbeta, everything(),names_to = "method", values_to = "RSN")
#if(save.fig)
#    pdf('beta_rmse_bee.pdf')
# beeswarm(RSN ~ method, data=data.RMSEbeta, col = 1:4,pch = 16, main=TeX("RMSE of $\\beta$"), cex.main=2, cex.lab=1.5)
#if(save.fig)
#    dev.off()

beta_violon <- ggplot(data.RMSEbeta, aes(x=method, y=RSN,fill=method)) +ggtitle(TeX("$\\beta$"))+
    geom_violin(trim=T) + stat_summary(fun=mean, geom="point", shape=23, size=2)+
    theme(legend.position="none",
         plot.title = element_text(size=30,face="bold",hjust = 0.5),
         axis.text.x  = element_text(size=20),
         axis.text.y  = element_text(size=20),
         axis.title.y = element_text(size=22,face="bold"),
         axis.title.x = element_blank())
         #axis.title.x = element_text(size=20,face="bold"))
print(beta_violon)
if(save.fig)
    ggsave("beta_violon.pdf",beta_violon, height=6, width=10)



data.RMSEbeta_smooth <- data.frame(Lasso   = MSE_beta_smooth.lasso,
                            Ridge   = MSE_beta_smooth.ridge,
                            HBayes  = MSE_beta_smooth.bayes,
                            Mix     = MSE_beta_smooth.fb,
                            Mr.ash  = MSE_beta_smooth_ash,
                            EMVS    = MSE_beta_smooth_EMVS,
                            horseshoe = MSE_beta_smooth_HS,
                            oracle = MSE_beta_smooth.oracle)
data.RMSEbeta_smooth <- pivot_longer(data.RMSEbeta_smooth, everything(),names_to = "method", values_to = "RSN")
#if(save.fig)
#    pdf('beta_rmse_smooth_bee.pdf')
#beeswarm(RSN ~ method, data=data.RMSEbeta, col = 1:4,pch = 16, main=TeX("RSN of $\\beta$ smoothed"), cex.main=2, cex.lab=1.5)
#if(save.fig)
#    dev.off()

beta_violon <- ggplot(data.RMSEbeta_smooth, aes(x=method, y=RSN,fill=method)) +ggtitle(TeX("$\\beta$ smoothed"))+
    geom_violin(trim=T) + stat_summary(fun=mean, geom="point", shape=23, size=2)+
    theme(legend.position="none",
          plot.title = element_text(size=30,face="bold",hjust = 0.5),
          axis.text.x  = element_text(size=20),
          axis.text.y  = element_text(size=20),
          axis.title.y = element_text(size=22,face="bold"),
          axis.title.x = element_blank())
print(beta_violon)
if(save.fig)
    ggsave("beta_smooth_violon.pdf",beta_violon, height=6, width=10)


data.XBbeta <- data.frame(Lasso = MSE_XB_lasso,
                             Ridge   = MSE_XB_ridge,
                             HBayes  = MSE_XB_npbayes,
                             Mix     = MSE_XB_fb,
                             Mr.ash  = MSE_XB_ash,
                             EMVS    = MSE_XB_EMVS,
                          horseshoe= MSE_XB_HS,
                          oracle = MSE_XB_oracle)
data.XBbeta <- pivot_longer(data.XBbeta, everything(),names_to = "method", values_to = "RSN")
#if(save.fig)
#    pdf('Xbeta_rmse_bee.pdf')
#beeswarm(RSN ~ method, data=data.XBbeta, col = 1:4,pch = 16, main=TeX("RSN of $X\\beta$"), cex.main=2, cex.lab=1.5)
#if(save.fig)
#    dev.off()
xbeta_violon <- ggplot(data.XBbeta, aes(x=method, y=RSN,fill=method)) +ggtitle(TeX("$X\\beta$"))+
    geom_violin(trim=T) + stat_summary(fun=mean, geom="point", shape=23, size=2)+  theme(legend.position="none",
                                                                                         plot.title = element_text(size=30,face="bold",hjust = 0.5),
                                                                                         axis.text.x  = element_text(size=20),
                                                                                         axis.text.y  = element_text(size=20),
                                                                                         axis.title.y = element_text(size=22,face="bold"),
                                                                                         axis.title.x = element_blank())

print(xbeta_violon)
if(save.fig)
    ggsave("xbeta_violon.pdf",xbeta_violon, height=6, width=10)


beta_true <- c(0.2,0.2,-0.2)
Table <- matrix(0, nrow=5, ncol = 3)
Table_smooth <- matrix(0, nrow=5, ncol = 3)
qtl_list <- list(beta.qtl.lasso,
                 beta.qtl.ridge,
                 beta.qtl.fb,
                 beta.qtl.bayes,
                 beta.qtl.oracle)
qtl_list_smooth <- list(beta.smooth.qtl.lasso,
                 beta.smooth.qtl.ridge,
                 beta.smooth.qtl.fb,
                 beta.smooth.qtl.bayes,
                 beta.smooth.qtl.oracle)
for(i in 1:5){

    for(ii in 1:3){
        Table[i, ii] <- sqrt(mean((qtl_list[[i]][,ii]-beta_true[ii])^2))
        Table_smooth[i, ii] <- sqrt(mean((qtl_list_smooth[[i]][,ii]-beta.smooth.qtl.true[,ii])^2))
    }
}

colnames(Table) <- c("$ {\\bf \\beta}_1$","$ {\\bf \\beta}_2$","$ {\\bf \\beta}_3$")
rownames(Table) <- c("lasso",  "ridge", "mix","hBeta","oracle")
colnames(Table_smooth) <- c("$ {\\bf \\beta}_1$","$ {\\bf \\beta}_2$","$ {\\bf \\beta}_3$")
rownames(Table_smooth) <- c("lasso",  "ridge", "mix","hBeta","oracle")
print.xtable(xtable(Table, digits=2),type="latex",sanitize.text.function = function(x) x)
print.xtable(xtable(Table_smooth, digits=2),type="latex",sanitize.text.function = function(x) x)

