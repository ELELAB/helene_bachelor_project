#Load DEG files
upreg_hugo <- readRDS(file = "../Upreg_hugo.RDS")
downreg_hugo <- readRDS(file = "../Downreg_hugo.RDS" )

upreg_hugo_oldGCn <- readRDS(file = "../Upreg_hugo_oldGCn.RDS" )
downreg_hugo_oldGCn <- readRDS(file = "../Downreg_hugo_oldGCn.RDS" )


#Find ID that are unique for oldGCn and find position
dif_upreg_oldGCn <- setdiff(upreg_hugo_oldGCn$ID ,upreg_hugo$ID )
pos_dif_upreg_oldGCn <- match(dif_upreg_oldGCn,upreg_hugo_oldGCn$ID )
saveRDS(dif_upreg_oldGCn, "dif_upreg_oldGCn.RDS")

dif_downreg_oldGCn <- setdiff(downreg_hugo_oldGCn$ID ,downreg_hugo$ID )
pos_dif_downreg_oldGCn <- match(dif_downreg_oldGCn,downreg_hugo_oldGCn$ID )
saveRDS(dif_downreg_oldGCn, "dif_downreg_oldGCn.RDS")

#Find ID that are unique for new method and find position
dif_upreg <- setdiff(upreg_hugo$ID ,upreg_hugo_oldGCn$ID )
pos_dif_upreg <- match(dif_upreg,upreg_hugo$ID )
saveRDS(dif_upreg, "dif_upreg.RDS")

dif_downreg <- setdiff(downreg_hugo$ID ,downreg_hugo_oldGCn$ID )
pos_dif_downreg <- match(dif_downreg,downreg_hugo$ID )
saveRDS(dif_downreg, "dif_downreg.RDS")

#extract logFC for unique oldGCn IDs
logFC_dif_upreg_oldGCn<-data.frame(ID = dif_upreg_oldGCn,
                             logFC_upreg = upreg_hugo_oldGCn$logFC[pos_dif_upreg_oldGCn])
logFC_dif_downreg_oldGCn <- data.frame(ID = dif_downreg_oldGCn,
                                       logFC_downreg = downreg_hugo_oldGCn$logFC[pos_dif_downreg_oldGCn])

#extract logFC for unique IDs
logFC_dif_upreg<-data.frame(ID = dif_upreg,
                                   logFC_upreg = upreg_hugo$logFC[pos_dif_upreg])
logFC_dif_downreg <- data.frame(ID = dif_downreg,
                                       logFC_downreg = downreg_hugo$logFC[pos_dif_downreg])


#plot logFC for unique oldGCn IDs

png('logFC_unique_oldGCn.png', width = 20,height =8,res = 300, units = "cm")
par(mfrow=c(1,2))
plot(x=seq(1:length(logFC_dif_upreg_oldGCn$ID )),y = logFC_dif_upreg_oldGCn$logFC_upreg,
     xlab = "unique IDs", ylab = "logFC", main = "Upregulated unique oldGCn"  )
abline(h = 0.5, col = "red", lwd = 3)
abline(h = mean(logFC_dif_upreg_oldGCn$logFC_upreg ), col = "blue", lwd = 3)

plot(x=seq(1:length(logFC_dif_downreg_oldGCn$ID )),y = logFC_dif_downreg_oldGCn$logFC_downreg,
     xlab = "unique IDs", ylab = "logFC", main = "Downregulated unique oldGCn"  )
abline(h = -0.5, col = "red", lwd = 3)
abline(h = mean(logFC_dif_downreg_oldGCn$logFC_downreg ), col = "blue", lwd = 3)

dev.off()

#plot logFC for unique IDs

png('logFC_unique_new.png', width = 20,height =8,res = 300, units = "cm")
par(mfrow=c(1,2))
plot(x=seq(1:length(logFC_dif_upreg$ID )),y = logFC_dif_upreg$logFC_upreg,
     xlab = "unique IDs", ylab = "logFC", main = "Upregulated unique new"  )
abline(h = 0.5, col = "red", lwd = 3)
abline(h = mean(logFC_dif_upreg$logFC_upreg ), col = "blue", lwd = 3)

plot(x=seq(1:length(logFC_dif_downreg$ID )),y = logFC_dif_downreg$logFC_downreg,
     xlab = "unique IDs", ylab = "logFC", main = "Downregulated unique new"  )
abline(h = -0.5, col = "red", lwd = 3)
abline(h = mean(logFC_dif_downreg$logFC_downreg ), col = "blue", lwd = 3)

dev.off()