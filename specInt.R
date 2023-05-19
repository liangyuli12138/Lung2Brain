# cell interaction
celldb = read.table('/public/workspace/lily/Lung2Brain/Version5/PairsLigRec.txt',header=T,sep='\t')

# e.g., specInt(inte123.expr,inte123$celltype2,celldb[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')],'3','myeloid')
specInt <- function(mat,label,intdb,L,R){
  rs = as.data.frame(t(apply(intdb,1,function(x){
    if(length(which(rownames(mat) %in% x)) == 2){
      tmp1 = tapply(unlist(mat[x[1],]),label,mean);
      tmp2 = tapply(unlist(mat[x[2],]),label,mean);
      idx1 = which(names(tmp1) %in% L);
      idx2 = which(names(tmp2) %in% R);
      sec1max = max(tmp1[-idx1],na.rm=T);
      sec2max = max(tmp2[-idx2],na.rm=T);
      sig1 = wilcox.test(unlist(mat[x[1],which(label %in% L)]),unlist(mat[x[1],-which(label %in% L)]),alternative='greater')$p.value
      sig2 = wilcox.test(unlist(mat[x[2],which(label %in% R)]),unlist(mat[x[2],-which(label %in% R)]),alternative='greater')$p.value
      c(mean(tmp1[idx1])/sec1max,mean(tmp2[idx2])/sec2max,mean(tmp1[idx1]),mean(tmp2[idx2]),mean(tmp1[-idx1],na.rm=T),mean(tmp2[-idx2],na.rm=T),sig1,sig2)
    }else{
      rep(NA,8)
    }
  })))
  colnames(rs) = c('L.Ratio','R.Ratio','L.mean','R.mean','nL.mean','nR.mean','L.pval','R.pval')
  rs[order(rs$L.Ratio*rs$R.Ratio,decreasing=T),]
}
#
#
# rs = specInt(inte123.expr[,which(inte123$celltype!='unknown')],inte123$seurat_clusters[which(inte123$celltype!='unknown')],celldb[,c('Ligand.ApprovedSymbol','Receptor.ApprovedSymbol')],c('4','9'),c('1','10','12','15','17'))
# rs.f = rs[which(rs$L.mean>0.1 & rs$R.mean>0.1),]
# rs.f$L.fc = rs.f$L.mean/rs.f$nL.mean
# rs.f$R.fc = rs.f$R.mean/rs.f$nR.mean
# head(rs.f[order(apply(rs.f[,c('L.fc','R.fc')],1,min),decreasing=T),],10)




