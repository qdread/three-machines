mlatent <- readRDS("three_machine_fit_latent_gammaprior.rds")
pdf('traceplots_rescor.pdf')
plot(mlatent, ask = FALSE)
dev.off()

library(dagitty)

g <- dagitty('dag {
    GAC [pos="0,1"]
    Vol [pos="1,1"]
    Perten [pos="2,1"]
    TWtrue [unobserved,pos="1,0"]
    Plot_Trial_Year [pos="1,-1"]

    Plot_Trial_Year -> TWtrue
    TWtrue -> GAC
    TWtrue -> Vol
    TWtrue -> Perten
}')
plot(g)
