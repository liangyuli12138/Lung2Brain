
# deal with E-MTAB-8659
# 2021-9-8
#===========================================================================================================================================
samplelist <- dir("./")
tmp.res <- matrix(NA,nrow=48107)
for(i in 1:length(samplelist)){
    tmp <- read.table(paste0("./",samplelist[i]),header=T)
    tmp.res <- cbind(tmp.res,tmp[,c(1,2)])
}

tmp.res <- tmp.res[,-1]
apply(tmp.res,2,function(x){all(as.vector(tmp.res[,1])==as.vector(x))}) # to check probe ID names 

res <- tmp.res[,grep("mean",colnames(tmp.res))]
res$probe <- tmp.res[,1]

tmp.ann <- data.table::fread("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/A-GEOD-10558_comments.txt")
ann.f <- as.data.frame(tmp.ann[,c("Comment[Symbol]","Comment[Probe_Id]")])
colnames(ann.f) <- c("Symbol","Probe")

# filter some gene 
res.f <- res[which(res$probe%in%ann.f$Probe),]
rs.dat <- merge(res.f,ann.f,by.x="probe",by.y="Probe")
rs.dat$probe <- NULL

rs.f <- aggregate(.~Symbol,data=rs.dat,FUN=median)

# trans date to gene name 
tmp.date <- read.table("/public/workspace/lily/software/date_gene.txt")
for(i in 1:nrow(tmp.date)){
    rs.f[grep(paste0("^",tmp.date[i,2],"$"),rs.f$Symbol),"Symbol"] <- as.vector(tmp.date[i,1])
}

rs.ff <- rs.f[-1,]
rownames(rs.ff) <- as.vector(rs.ff$Symbol)
rs.ff$Symbol <- NULL

saveRDS(rs.ff,file="/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")




############################################################################################################################################
# 2021-9-8 
# analysis in a new data ,what need to do 
# 1. MDM and MG in sample
# 2. BMS Lung and Brain like tumor score 

dat <- readRDS("/public/workspace/lily/metastasis/data/verify/E-MTAB-8659/LCBM_S63_data.RDS")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod <- mod.analyze2(as.matrix(dat),c("BMS_update","Brain_gene","Lung_gene","BMDM_marker","MG_marker","RES","SEN","immunecheck"),"/public/workspace/lily/MOD_file/",permN=0)
mod <- as.data.frame(mod)




























































