GSE113079=getGEOExpData('GSE113079')
GSE20681=getGEOExpData('GSE20681')
GSE64566=getGEOExpData('GSE64566')
row.names(GSE20686$Exp$GPL4133_45015_Data_col1)
GPL4133.reanno=re_annotation_by_probe(GSE20686$Anno$GPL4133[,c(1,20)])
GPL6947.reanno=re_annotation_by_probe(GSE64566$Anno$GPL6947[,c(1,19)])
pob=GSE113079$Anno$GPL20115[,c(1,6)]
pob=pob[which(pob[,2]!=''),]
head(pob)
GPL20115.anno=re_annotation_by_probe(pob)
#GPL20115.anno=pob.anno
head(GPL4133.reanno$match_data)
head(GPL20115.anno$match_uniq)
GSE113079.exp=exp_probe2symbol_v2(GSE113079$Exp$GPL20115_65523_Data_col1,anno = crbind2DataFrame(GPL20115.anno$match_uniq[,c(1,3)]))
GSE20681.exp=exp_probe2symbol_v2(GSE20681$Exp$GPL4133_45015_Data_col1,anno = crbind2DataFrame(GPL4133.reanno$match_uniq[,c(1,3)]))
GSE64566.exp=exp_probe2symbol_v2(GSE64566$Exp$GPL6947_48783_Data_col1,anno = crbind2DataFrame(GPL6947.reanno$match_uniq[,c(1,3)]))

GPL20115.anno.lnc=as.character(crbind2DataFrame(GPL20115.anno$match_uniq[which(GPL20115.anno$match_uniq[,6]=='TRUE'),])[,3])
GPL20115.anno.gene=as.character(crbind2DataFrame(GPL20115.anno$match_uniq[which(GPL20115.anno$match_uniq[,6]=='FALSE'),])[,3])

GPL4133.anno.lnc=as.character(crbind2DataFrame(GPL4133.reanno$match_uniq[which(GPL4133.reanno$match_uniq[,6]=='TRUE'),])[,3])
GPL4133.anno.gene=as.character(crbind2DataFrame(GPL4133.reanno$match_uniq[which(GPL4133.reanno$match_uniq[,6]!='TRUE'),])[,3])

GPL6947.anno.lnc=as.character(crbind2DataFrame(GPL6947.reanno$match_uniq[which(GPL6947.reanno$match_uniq[,6]=='TRUE'),])[,3])
GPL6947.anno.gene=as.character(crbind2DataFrame(GPL6947.reanno$match_uniq[which(GPL6947.reanno$match_uniq[,6]!='TRUE'),])[,3])

GSE113079.exp.lnc=GSE113079.exp[row.names(GSE113079.exp)%in%GPL20115.anno.lnc,]
GSE113079.exp.gene=GSE113079.exp[row.names(GSE113079.exp)%in%GPL20115.anno.gene,]
dim(GSE113079.exp.gene)
dim(GSE113079.exp.lnc)
GSE113079$Sample$DataProcessing
quantile(GSE113079.exp.gene[,1])

GSE113079.group=gsub(' .*','',gsub('^ ','',gsub('Peripheral Blood Mononuclear Cells,','',GSE113079$Sample$Title)))
table(GSE113079.group)
GSE113079.group=ifelse(GSE113079.group=='CAD','CAD','Control')
#GSE20686.exp.lnc=GSE20686.exp[row.names(GSE20686.exp)%in%GPL4133.anno.lnc,]
#GSE20686.exp.gene=GSE20686.exp[row.names(GSE20686.exp)%in%GPL4133.anno.gene,]
head(GSE20681$Sample)
GSE20681.group=ifelse(GSE20681$Sample$`disease state`=='Case (1)','CAD','Control')
GSE20681.exp.lnc=GSE20681.exp[row.names(GSE20681.exp)%in%GPL4133.anno.lnc,]
GSE20681.exp.gene=GSE20681.exp[row.names(GSE20681.exp)%in%GPL4133.anno.gene,]
GSE20681$Sample$DataProcessing
dim(GSE20681.exp.lnc)
dim(GSE20681.exp.gene)
quantile(GSE20681.exp.gene[,1])
dim(GSE20686.exp.gene)
dim(GSE20686.exp.lnc)
intersect(gene.up,lnc.up)

GSE64566$Sample$`disease state`
GSE64566.exp.sample=GSE64566$Sample[match(colnames(GSE64566.exp),GSE64566$Sample$Acc),]
GSE64566.group=ifelse(GSE64566.exp.sample$`disease state`=='control','Control','CAD')
GSE64566.exp.lnc=log2(GSE64566.exp[row.names(GSE64566.exp)%in%GPL6947.anno.lnc,])
GSE64566.exp.gene=log2(GSE64566.exp[row.names(GSE64566.exp)%in%GPL6947.anno.gene,])
dim(GSE64566.exp.gene)
dim(GSE64566.exp.lnc)
log2(GSE64566.exp.gene)
boxplot(GSE64566.exp.gene)
boxplot(GSE64566.exp.lnc)
mg_pca_plot(t(GSE64566.exp.gene),group = GSE64566.group)
mg_limma_DEG()

dim(GSE113079.exp.gene)
dim(GSE20681.exp.gene)
dim(GSE64566.exp.gene)
boxplot(GSE113079.exp.gene,outline=F)
boxplot(GSE20681.exp.gene,outline=F)
boxplot(GSE64566.exp.gene,outline=F)

lay1 <- customLayout::lay_new(matrix(1:3,nrow = 1),widths = c(141,198,46))
lay2 <- customLayout::lay_new(matrix(1:3,nrow = 1),widths = c(141,198,46))
#lay1 <- customLayout::lay_new(matrix(1:2,nrow = 1))

cl <- customLayout::lay_bind_row(lay2, lay1,heights = c(1,1)) 
customLayout::lay_show(cl)
#customLayout::lay_show(lay2)

#lay3 <- customLayout::lay_new(matrix(1:2))
#cl2 <- customLayout::lay_bind_row(cl, lay3, heights = c(1, 1)) 
#customLayout::lay_show(cl)
pdf('Fig2A.pdf',width = 20,height = 8)
customLayout::lay_set(cl) 
boxplot(GSE113079.exp.gene,outline=F,col=ifelse(GSE113079.group=='Control','blue','red')
        ,ylab='Gene Expression',main='GSE113079')
boxplot(GSE20681.exp.gene,outline=F,col=ifelse(GSE20681.group=='Control','blue','red')
        ,ylab='Gene Expression',main='GSE20681')
boxplot(GSE64566.exp.gene,outline=F,col=ifelse(GSE64566.group=='Control','blue','red')
        ,ylab='Gene Expression',main='GSE64566')

boxplot(GSE113079.exp.lnc,outline=F,col=ifelse(GSE113079.group=='Control','blue','red')
        ,ylab='lncRNA Expression',main='GSE113079')
boxplot(GSE20681.exp.lnc,outline=F,col=ifelse(GSE20681.group=='Control','blue','red')
        ,ylab='lncRNA Expression',main='GSE20681')
boxplot(GSE64566.exp.lnc,outline=F,col=ifelse(GSE64566.group=='Control','blue','red')
        ,ylab='lncRNA Expression',main='GSE64566')
dev.off()
pdf('Fig2C.pdf',width = 4,height = 4)
mg_venn_plot(list(GSE113079=row.names(GSE113079.exp.lnc)
                  ,GSE20681=row.names(GSE20681.exp.lnc)
                  ,GSE64566=row.names(GSE64566.exp.lnc)
))
dev.off()
head(GPL6947.anno.lnc)
GSE113079.lnc.limma=mg_limma_DEG(GSE113079.exp.lnc,group = GSE113079.group,ulab = 'CAD',dlab = 'Control')
GSE113079.lnc.limma$Summary
GSE113079.pcg.limma=mg_limma_DEG(GSE113079.exp.gene,group = GSE113079.group,ulab = 'CAD',dlab = 'Control')

lnc.up=row.names(GSE113079.lnc.limma$DEG)[which(GSE113079.lnc.limma$DEG$logFC>1&GSE113079.lnc.limma$DEG$adj.P.Val<0.05)]
lnc.down=row.names(GSE113079.lnc.limma$DEG)[which(GSE113079.lnc.limma$DEG$logFC< -1&GSE113079.lnc.limma$DEG$adj.P.Val<0.05)]
gene.up=row.names(GSE113079.pcg.limma$DEG)[which(GSE113079.pcg.limma$DEG$logFC>1&GSE113079.pcg.limma$DEG$adj.P.Val<0.05)]
gene.down=row.names(GSE113079.pcg.limma$DEG)[which(GSE113079.pcg.limma$DEG$logFC< -1&GSE113079.pcg.limma$DEG$adj.P.Val<0.05)]

fig1d=mg_volcano(GSE113079.lnc.limma$DEG$logFC,pvalue = GSE113079.lnc.limma$DEG$adj.P.Val)
fig1e=mg_volcano(GSE113079.pcg.limma$DEG$logFC,pvalue = GSE113079.pcg.limma$DEG$adj.P.Val)
fig1de=mg_merge_plot(fig1d,fig1e,ncol = 2,nrow = 1,labels = c('D','E'),common.legend = T)

fig1f=mg_pca_plot(t(GSE113079.exp.lnc[match(c(lnc.up,lnc.down),row.names(GSE113079.exp.lnc)),]),GSE113079.group)
fig1g=mg_pca_plot(t(GSE113079.exp.gene[match(c(gene.up,gene.down),row.names(GSE113079.exp.gene)),]),GSE113079.group)
fig1fg=mg_merge_plot(fig1f$Plot,fig1g$Plot,ncol = 2,nrow = 1,labels = c('F','G'),common.legend = T)
fig1defg=mg_merge_plot(fig1de,fig1fg,ncol = 2,nrow=1)
savePDF('Fig2DEFG.pdf',fig1defg,width = 10,height = 4)
#head(GSE113079.exp)
#head(GSE20686.exp)
#head(GSE64566.exp)
length(c(lnc.up,lnc.down))
length(c(gene.up,gene.down))

#################WGCNA###################
wgcna.input=rbind(GSE113079.pcg.limma$Exp,GSE113079.lnc.limma$Exp)
wgcna.tag.group=c(rep('PCG',nrow(GSE113079.pcg.limma$Exp)),rep('lncRNA',nrow(GSE113079.lnc.limma$Exp)))
wgcna.inds=which(apply(wgcna.input, 1,sd)>quantile(apply(wgcna.input, 1,sd))['25%'])
wgcna.input=wgcna.input[wgcna.inds,]
wgcna.tag.group=wgcna.tag.group[wgcna.inds]
names(wgcna.tag.group)=row.names(wgcna.input)
par(mfrow=c(1,3))
boxplot(apply(GSE113079.pcg.limma$Exp, 1,sd),outline=F)
boxplot(apply(GSE113079.lnc.limma$Exp, 1,sd),outline=F)
boxplot(apply(wgcna.input, 1,sd),outline=F)

names(wgcna.tag.group)=row.names(wgcna.input)

wgcna.input.power=mg_wgcna_get_power(t(wgcna.input),blockSize = 30000)
wgcna.input.module=mg_WGCNA_getModule(t(wgcna.input),power = wgcna.input.power$cutPower,deepSplit = 4,blockSize = 30000)
table(wgcna.input.module$Modules[,2])
#wgcna.input.module$plotTOM

wgcna.input.module$Modules[,2]
match(row.names(wgcna.input.module$Modules),names(wgcna.tag.group))
modules=cbind(wgcna.input.module$Modules,GeneType=wgcna.tag.group)
head(modules)
table(modules[,2],modules[,3])
pdf(file = 'Fig3B.pdf',width = 8,height = 5)
plotDendroAndColors(wgcna.input.module$Tree, cbind(wgcna.input.module$Modules,GeneType=ifelse(wgcna.tag.group=='PCG','red','blue')),
                    c("Dynamic Module",'Merged Module','RNAType'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

#TOMplot(wgcna.input.module$plotTOM, wgcna.input.module$Tree, wgcna.input.module$Modules[,2], main = "Network heatmap plot, all genes")

writeMatrix(cbind(wgcna.input.module$Modules,Type=wgcna.tag.group),outpath = paste0('Module_info.txt'))
dat_count=crbind2DataFrame(table(wgcna.input.module$Modules[,2],wgcna.tag.group))
dat_count_deg=rbind()
for(m in row.names(dat_count)){
  gs=row.names(wgcna.input.module$Modules)[wgcna.input.module$Modules[,2]==m]
  dat_count_deg=rbind(dat_count_deg,c(length(intersect(lnc.up,gs)),length(intersect(lnc.down,gs))
  ,length(intersect(gene.up,gs)),length(intersect(gene.down,gs))))
}
colnames(dat_count_deg)=c('Lnc_Up','Lnc_Down','Gene_Up','Gene_Down')
dat_count=cbind(dat_count,dat_count_deg)
dim(dat_count)
#dat_count=cbind(Module=names(table(wgcna.input.module$Modules[,2])),Count=table(wgcna.input.module$Modules[,2]))
writeMatrix(dat_count,outpath = paste0('Module_summary.txt'),row = F)

MEs_col = wgcna.input.module$MEs
colnames(MEs_col)=gsub('^ME','',colnames(MEs_col))
MEs_col = orderMEs(MEs_col)


all.cor=rbind()
for(i in 1:(nrow(dat_count)-1)){
  #m=wgcna.tree$Modules[1,2]
  m=gsub('^ME','',colnames(wgcna.input.module$MEs)[i])
  me=wgcna.input.module$MEs[,i]
  dat.exp=t(wgcna.input[which(wgcna.input.module$Modules[,2]==m),])
  #wgcna.input
  #print(dat.exp)
  me.cor=psych::corr.test(dat.exp,me,method = 'spearman')
  #me.cor$r
  print(paste0(m,':',mean(abs(me.cor$r[,1]))))
  all.cor=rbind(all.cor,cbind(Gene=row.names(me.cor$r),Group=rep(m,nrow(me.cor$r))
                              ,facet=rep(i,nrow(me.cor$r))
                              ,R=abs(me.cor$r[,1]),P=me.cor$p[,1]))
}
head(all.cor)

all.cor.means=crbind2DataFrame(all.cor) %>% dplyr::group_by(Group) %>% dplyr::summarise_each(dplyr::funs(median),vals=c('R'))
all.cor.means=crbind2DataFrame(all.cor.means)
all.cor.means[order(all.cor.means[,2]),]
r.module=all.cor.means[all.cor.means[,2]>mean(all.cor.means[,2]),1]
deg.module=c('lightgreen','orangered4','brown','lightsteelblue1','grey60','salmon4'
             ,'darkgreen','palevioletred3','darkturquoise','mediumpurple3'
             ,'plum1','skyblue','darkred','blue','darkorange')
r.cad.module=c('lightgreen','orangered4','brown','lightsteelblue1','grey60','darkgreen','darkturquoise','mediumpurple3'
               ,'skyblue','darkred','blue','darkorange')
mg_venn_plot(list(A=r.module,B=r.cad.module,C=deg.module))
sig.module=intersect(intersect(r.module,deg.module),r.cad.module)
dat_count[order(dat_count[,3]),]
dat_count[row.names(dat_count)%in%sig.module,]
paste0(sig.module,collapse = ',')
cad.color='#F8766D'
control.color='#00BFC4'
fig3c=mg_ridges_plot(crbind2DataFrame(all.cor[,c(2,2,4)]),melt = T,xlab = '|Spearman correlations|')
savePDF(fig3c,filename = 'Fig3C.pdf',width = 7,height = 6)
colnames(MEs_col)

mes.cad.r=mg_muti_cor_plot(MEs_col[,-18],cbind(CAD=ifelse(GSE113079.group=='CAD',1,0)
                               ,Control=ifelse(GSE113079.group!='CAD',1,0)))
mes.cad.r$r

fig3d=mg_quick_cor_plot(MEs_col[,-18],cbind(CAD=ifelse(GSE113079.group=='CAD',1,0)
                                      ,Control=ifelse(GSE113079.group!='CAD',1,0)))
pdf('Fig3D.pdf',width = 3,height = 6)
fig3d$plot
dev.off()
fig3e=groupViolin(MEs_col[,-18],group = GSE113079.group,ylab = 'Module Eigengene',group_col = c(cad.color,control.color))
#mg_merge_plot(fig3c,fig3d,ncol = 2,nrow=1)
savePDF(fig3e,filename = 'Fig3E.pdf',width = 8,height = 4)

#####################Enrichment############
dat_count
sig.module=intersect(intersect(r.module,deg.module),r.cad.module)
dat_count[row.names(dat_count)%in%sig.module,]



mg_venn_plot(
list(grey60=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='grey60')]
     ,lightsteelblue1=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='lightsteelblue1')]
     ,mediumpurple3=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='mediumpurple3')]
     ,orangered4=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='orangered4')]
     ,UP_Gene=gene.up
     ,UP_LncRNA=lnc.up
))
mg_venn_plot(
  list(grey60=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='grey60')]
       ,lightsteelblue1=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='lightsteelblue1')]
       ,mediumpurple3=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='mediumpurple3')]
       ,orangered4=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]=='orangered4')]
       ,Down_Gene=gene.down
       ,Down_LncRNA=lnc.down
  ))

all_enrich=rbind()
for(m in sig.module){
  #m=sig.module[1]
  gs=row.names(wgcna.input.module$Modules)[which(wgcna.input.module$Modules[,2]==m)]
  gs.enrich=mg_clusterProfiler(gs)
  gs.enrich$Enrich_tab$Module=rep(m,nrow(gs.enrich$Enrich_tab))
  all_enrich=rbind(all_enrich,gs.enrich$Enrich_tab)
  #table(gs.enrich$Enrich_tab$DB)
}
all_enrich.cnt=crbind2DataFrame(table(all_enrich$DB,all_enrich$Module))
row.names(all_enrich.cnt)=gsub('geneontology_','',row.names(all_enrich.cnt))
pdf('Fig4A.pdf',width = 6,height = 2.5)
mheatmap(all_enrich.cnt,display_numbers = T,cluster_rows = F
         ,cluster_cols = F,min_col = 'white',med_col = 'orange'
         ,column_names_rot = 15,legend = 'Term Count')
dev.off()
pdf('Fig4B.pdf',width = 2.5,height = 2.5)
mg_venn_plot(
  list(Module_Genes=row.names(wgcna.input.module$Modules)[wgcna.input.module$Modules[,2]%in%sig.module]
       ,UP_Gene=gene.up
       ,UP_LncRNA=lnc.up
       ,Down_Gene=gene.down
       ,Down_LncRNA=lnc.down))
dev.off()

gs=intersect(row.names(wgcna.input.module$Modules)[wgcna.input.module$Modules[,2]%in%sig.module],
          c(gene.up,gene.down,lnc.down,lnc.up))
gene.up[-grep('\\|',gene.up)]
gene.down
gs.enrich=mg_clusterProfiler(gs)
table(gs.enrich$Enrich_tab$DB)

fig4c=mg_clusterProfiler_dotplot(gs.enrich$KEGG)
fig4d=mg_clusterProfiler_dotplot(gs.enrich$GO_MF)
fig4cd=mg_merge_plot(fig4c,fig4d,ncol = 2,nrow=1,labels = c('C','D'))
savePDF('Fig4CD.pdf',fig4cd,width = 10,height = 3)

gs=intersect(row.names(wgcna.input.module$Modules)[wgcna.input.module$Modules[,2]%in%sig.module],
             gene.down)
gs.enrich=mg_clusterProfiler(gs)
table(gs.enrich$Enrich_tab$DB)

gs
all.gano=crbind2DataFrame(GPL20115.anno$match_uniq)
head(all.gano)
gs.anno=all.gano[match(gs,all.gano[,3]),2:6]
gs.anno$Module=wgcna.input.module$Modules[match(gs,row.names(wgcna.input.module$Modules)),2]
head(gs.anno)
dim(wgcna.input.module$Modules)
dir(paste0(MG_Grobal_baseFolder,'/source/genome_info'))

bg.gff.range=mg_gff_to_GRanges(paste0(MG_Grobal_baseFolder,'/source/genome_info/gencode.v33.annotation.gff3'))
gene.v33.tab=readMatrix(paste0(MG_Grobal_baseFolder,'/source/gencode.v33.id.tab.txt'))
head(gene.v33.tab)
gene.pos.anno=bg.gff.range[bg.gff.range$type=='gene']

gs.anno
gs.gene.pos.anno=gene.pos.anno[gsub('\\..*','',gene.pos.anno$gene_id)%in%gs.anno[,1]]
length(gs.gene.pos.anno)
gs.gene.pos.anno
as.character(strand(gs.gene.pos.anno))
as.character(seqnames(gs.gene.pos.anno))
gs.gene.pos.anno.1=cbind(ENSG=gsub('\\..*','',gs.gene.pos.anno$gene_id),Chr=as.character(seqnames(gs.gene.pos.anno)),Start=start(gs.gene.pos.anno),END=end(gs.gene.pos.anno),Strand=as.character(strand(gs.gene.pos.anno)))

lst.gs.anno=gs.anno[gs.anno[,1]%in%gsub('\\..*','',gene.pos.anno$gene_id),]
lst.gs.anno=cbind(lst.gs.anno,gs.gene.pos.anno.1[match(lst.gs.anno[,1],gs.gene.pos.anno.1[,1]),])
table(lst.gs.anno[,6],lst.gs.anno[,5])
dim(lst.gs.anno)
head(lst.gs.anno)
get_cis_or_trans=function(lnc,pcg){
  if(lnc[8]!=pcg[8]){
    return('Trans')
  }
  if(lnc[11]=='+'&pcg[11]=='+'){
    if(lnc[9]<pcg[9]){
      return('Cis')
    }
  }else if(lnc[11]=='-'&pcg[11]=='-'){
    if(lnc[10]>pcg[10]){
      return('Cis')
    }
  }
  return('Trans')
}
lst.gs.anno=crbind2DataFrame(lst.gs.anno)
all.cis.trans=rbind()
for(m in unique(lst.gs.anno[,6])){
  lnc.anno=which(lst.gs.anno[,6]==m&lst.gs.anno[,5]=='TRUE')
  pcg.anno=which(lst.gs.anno[,6]==m&lst.gs.anno[,5]=='FALSE')
  if(length(lnc.anno)>0&length(pcg.anno)>0){
    for(i in lnc.anno){
      for(j in pcg.anno){
        c_ty=get_cis_or_trans(lst.gs.anno[i,],lst.gs.anno[j,])
        all.cis.trans=rbind(all.cis.trans,c(lst.gs.anno[i,1:2],lst.gs.anno[j,1:2],lst.gs.anno[i,8:11],lst.gs.anno[j,8:11],c_ty))
      }
    }
  }
}
all.cis.trans=crbind2DataFrame(all.cis.trans)
head(all.cis.trans)
dim(all.cis.trans)
table(all.cis.trans[,13])

load(paste0(MG_Grobal_baseFolder,'/source/MIMAT2MG_ID.RData'))
miRNA_mRNA=unique(MIMAT2MG_ID[MIMAT2MG_ID[,3]%in%c('miRanda','miRtarbase','TargetScan','starbase'),1:2])
dim(miRNA_mRNA)
head(miRNA_mRNA)
all.cis.trans[,3]
all.cis.trans.gene.mgid=mg_idconvert_local(all.cis.trans[,3])
miRNA_mRNA=miRNA_mRNA[miRNA_mRNA[,2]%in%all.cis.trans.gene.mgid$IDMap$MG_ID,]
miRNA_mRNA=crbind2DataFrame(miRNA_mRNA)
miRNA_mRNA$ENSG=all.cis.trans.gene.mgid$IDMap$SearchID[match(miRNA_mRNA[,2],all.cis.trans.gene.mgid$IDMap$MG_ID)]
head(miRNA_mRNA)
miRNA_lncRNA=readMatrix('starBaseV3_hg19_CLIP-seq_mi_lncRNA.txt',row=F)
head(miRNA_lncRNA)
miRNA_lncRNA=miRNA_lncRNA[miRNA_lncRNA[,3]%in%all.cis.trans[,1],]
dim(miRNA_lncRNA)
head(miRNA_lncRNA)
all.lnc.fisher=rbind()
all.miRNA=c()
for(i in 1:nrow(all.cis.trans)){
  lnc=all.cis.trans[i,1]
  pcg=all.cis.trans[i,3]
  m.lnc.mirna=unique(miRNA_lncRNA[which(miRNA_lncRNA[,3]==lnc),1])
  all.miRNA=unique(c(all.miRNA,m.lnc.mirna))
  m.pcg.mirna=unique(miRNA_mRNA[which(miRNA_mRNA[,3]==pcg),1])
  all.miRNA=unique(c(all.miRNA,m.pcg.mirna))
  all.lnc.fisher=rbind(all.lnc.fisher,c(lnc,pcg,length(m.lnc.mirna)
                                        ,length(m.pcg.mirna),length(intersect(m.lnc.mirna,m.pcg.mirna))))
}
all.miRNA=unique(all.miRNA)
colnames(all.lnc.fisher)=c('lncRNA','Protein','lnc_target','pcg_target','common')
all.lnc.fisher=crbind2DataFrame(all.lnc.fisher)
ps=apply(all.lnc.fisher, 1, function(x){
  l=as.numeric(x[3])
  p=as.numeric(x[4])
  a=as.numeric(x[5])
  if(a==0){p=1}else{
    p=fisher.test(matrix(c(a,p-a,l-a,length(all.miRNA)),2,2))$p.value
  }
  return(p)
})
all.lnc.fisher$p.value=ps
head(all.lnc.fisher)
all.cis.trans.merge=cbind(all.cis.trans,all.lnc.fisher[,3:6])
head(all.cis.trans.merge)
all.cis.trans.merge=all.cis.trans.merge[all.lnc.fisher[,6]<0.05,]
head(all.cis.trans.merge)
all.lnc.fisher[all.lnc.fisher[,4]>0,]

writeMatrix(all.cis.trans.merge,'Figures/TableS1.txt',row = F)
sum(all.lnc.fisher[,6]<0.05)


par(mar=c(2, 2, 2, 2));
library(RCircos)
data(UCSC.HG38.Human.CytoBandIdeogram)
library(OmicCircos)
db <- segAnglePo(seg.dat = crbind2DataFrame(UCSC.HG38.Human.CytoBandIdeogram)
                 ,seg = as.character(unique(UCSC.HG38.Human.CytoBandIdeogram[,1]))
);
pdf('Fig5A.pdf',width = 8,height = 8)
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="");
circos(R=300, cir=db, W=15,mapping = UCSC.HG38.Human.CytoBandIdeogram,type="chr2"
       ,print.chr.lab=TRUE, scale=T
       ,cex = 1.5);#最外层染色体信息

lb=unique(rbind(cbind(all.cis.trans.merge[,c(2,5:7)],Type=rep('Lnc',nrow(all.cis.trans.merge))),
      cbind(all.cis.trans.merge[,c(4,9:11)],Type=rep('PCG',nrow(all.cis.trans.merge)))))
table(lb[,5])
lb.dt=crbind2DataFrame(cbind(as.character(lb[,2]),ceiling(apply(lb[,c(3,4)],1,mean))
                       ,lb[,1],as.character(lb[,5])))
lb.dt=lb.dt[order(as.numeric(gsub('chr','',gsub('Y','24',gsub('X','23',lb.dt[,1]))))),]
lb.dt.reorder=rbind()
for(cr in unique(lb.dt[,1])){
  sb=lb.dt[which(lb.dt[,1]==cr),]
  sb=sb[order(sb[,2]),]
  lb.dt.reorder=rbind(lb.dt.reorder,sb)
}

circos (R=320 ,cir=db , W=45, mapping=lb.dt.reorder, type="label" ,
        side="out" 
        , col=mg_colors[as.factor(lb.dt.reorder[,4])]
        , cex=0.5 ,B = F
        ) ;
lnk.dt=all.cis.trans.merge[,c(c(5:7,9:11))]
circos (R=280 ,cir=db , W=40, mapping=lnk.dt, type="link.pg" ,
        col=mg_colors[as.numeric(as.factor(all.cis.trans.merge[,13]))+2] , lwd=2) ;

dev.off()

lrm.cad.enrich=mg_clusterProfiler(all.cis.trans.merge[,c(4)])
lrm.cad.enrich$Enrich_tab
lrm.cad.lnc.score=ssGSEAScore_by_genes(GSE113079.exp.lnc,genes =  all.cis.trans.merge[,c(2)])
lrm.cad.gene.score=ssGSEAScore_by_genes(GSE113079.exp.gene,genes = all.cis.trans.merge[,c(4)])
GSE113079.group

fig5b=mg_violin(crbind2DataFrame(cbind(GSE113079.group,Score=lrm.cad.lnc.score[1,])),melt = T,ylab = 'LncRNAs Enrichment score')
fig5c=mg_violin(crbind2DataFrame(cbind(GSE113079.group,Score=lrm.cad.gene.score[1,])),melt = T,ylab = 'Genes Enrichment score')
fig5bc=mg_merge_plot(fig5c,fig5b,ncol=2,nrow=1,common.legend = T)
savePDF('fig5bc.pdf',fig5bc,width = 4,height = 4)
kegg.ssgsea=ssGSEA_batch(GSE113079.exp.gene,dbs='KEGG')
kegg.ssgsea$KEGG

cor(t(lrm.cad.lnc.score),t(kegg.ssgsea$KEGG))

es.kegg.cor=mg_muti_cor_plot(t(kegg.ssgsea$KEGG),cbind(GeneES=lrm.cad.gene.score[1,],LncES=lrm.cad.lnc.score[1,]))
es.sig.kegg.name=c(row.names(es.kegg.cor$r[which(p.adjust(es.kegg.cor$p[,1])<0.05),]),row.names(es.kegg.cor$r[which(p.adjust(es.kegg.cor$p[,2])<0.05),]))
es.sig.kegg.es=t(kegg.ssgsea$KEGG[match(es.sig.kegg.name
                         ,row.names(kegg.ssgsea$KEGG)),])
colnames(es.sig.kegg.es)=gsub('KEGG_','',colnames(es.sig.kegg.es))
es.sig.kegg.es.plot=mg_quick_cor_plot(es.sig.kegg.es,cbind(GeneES=lrm.cad.gene.score[1,],LncES=lrm.cad.lnc.score[1,]))
pdf('Fig5D.pdf',width = 4,height = 4)
es.sig.kegg.es.plot$plot
dev.off()

GSE113079.exp.mg=rbind(GSE113079.exp.lnc,GSE113079.exp.gene)
GSE113079.exp.mg[match(as.character(all.cis.trans.merge[1,c(2,4)]),row.names(GSE113079.exp.mg)),]

all.lda.prd.score=rbind()
all.lnc.score=rbind()
all.gene.score=rbind()
roc.mx=c()
for(i in 1:nrow(all.cis.trans.merge)){
  pdt=scale(t(GSE113079.exp.mg[match(as.character(all.cis.trans.merge[i,c(2,4)]),row.names(GSE113079.exp.mg)),]))
lda.md=mg_lda_model(dat = pdt
             ,group = GSE113079.group)
lda.md.prd=mg_lda_predict(lda.md$Model,pdt
               ,group = GSE113079.group)
roc.mx=c(roc.mx,max(c(lda.md.prd$ROC$CAD$AUC,lda.md.prd$ROC$Control$AUC)))
all.lda.prd.score=rbind(all.lda.prd.score,cbind(Score=lda.md.prd$PredVal[,1],Group=GSE113079.group,LRMCAD=paste0(all.cis.trans.merge[i,2],'-',all.cis.trans.merge[i,4])))

lnc.auc=mg_auto_cal_AUC(pdt[,1],GSE113079.group)
lac=1
if(lnc.auc[1,2]<0.5){
  lac=-1
}
gene.auc=mg_auto_cal_AUC(pdt[,2],GSE113079.group)
gac=1
if(gene.auc[1,2]<0.5){
  gac=-1
}
all.lnc.score=rbind(all.lnc.score,cbind(Score=pdt[,1]*lac,Group=GSE113079.group,LNC=all.cis.trans.merge[i,2]))
all.gene.score=rbind(all.gene.score,cbind(Score=pdt[,2]*gac,Group=GSE113079.group,LNC=all.cis.trans.merge[i,4]))
}
mean(roc.mx)
fig6a=mg_auto_ggplot_roc(all.lda.prd.score[,1],label = all.lda.prd.score[,2]
                   ,group = all.lda.prd.score[,3],col='rawbow')+mg_get_ggplot_legend('right')
mean(mg_auto_cal_AUC(all.lda.prd.score[,1],label = all.lda.prd.score[,2]
                     ,group = all.lda.prd.score[,3])[,2])

fig6b=mg_auto_ggplot_roc(all.lnc.score[,1],label = GSE113079.group
                         ,group = all.lnc.score[,3],col='rawbow')
mean(mg_auto_cal_AUC(all.lnc.score[,1],label = GSE113079.group
                ,group = all.lnc.score[,3])[,2])

fig6c=mg_auto_ggplot_roc(all.gene.score[,1],label = GSE113079.group
                         ,group = all.gene.score[,3],col='rawbow')+mg_get_ggplot_legend('right')

fig6bc=mg_merge_plot(fig6b,fig6c,ncol = 2,nrow = 1,labels = c('A','B'),widths = c(0.5,1))

mean(mg_auto_cal_AUC(all.gene.score[,1],label = GSE113079.group
                     ,group = all.gene.score[,3])[,2])

fig6abc=mg_merge_plot(fig6bc,fig6a,nrow=2,ncol = 1,labels = c('','C'))
savePDF('Fig6ABC.pdf',fig6abc,width = 12,height = 12)

cad.disgene.evd=readMatrix('C1956346_disease_gda_evidences.tsv',row=F)
cad.disgene.gda=readMatrix('C1956346_disease_gda_summary.tsv',row=F)
head(cad.disgene.evd)
head(cad.disgene.gda)

cad.disgene.evd[cad.disgene.evd[,3]%in%intersect(cad.disgene.evd[,3],all.cis.trans.merge[,4]),]
cad.disgene.gda[cad.disgene.gda[,3]%in%intersect(cad.disgene.gda[,3],all.cis.trans.merge[,4]),]
t1=proc.time()
drug.dis=mg_gene2drug_distances(all.cis.trans.merge[,4])
#drug.dis$plot
t2=proc.time()
t2-t1

writeMatrix(sig.drug,outpath = 'Figures/Table2.txt')
savePDF('Fig7A.pdf',drug.dis$plot,width = 4,height = 4)
sig.drug=drug.dis$Data[which(drug.dis$Data[,9]<0.05),]
dim(drug.dis$Durg2Gene)
paste0(sig.drug[,1],collapse = ',')
drug.dis$Durg2Gene[,colnames(drug.dis$Durg2Gene)%in%sig.drug[,1]]

all.cis.trans.merge[all.cis.trans.merge[,4]%in%c('S1PR1'),]
intersect(row.names(GSE113079.exp),all.cis.trans.merge[all.cis.trans.merge[,4]%in%c('S1PR1'),c(2,4)])
intersect(row.names(GSE20681.exp),all.cis.trans.merge[all.cis.trans.merge[,4]%in%c('S1PR1'),c(2,4)])
intersect(row.names(GSE64566.exp),all.cis.trans.merge[all.cis.trans.merge[,4]%in%c('S1PR1'),c(2,4)])

which(row.names(GSE20681.exp)=='S1PR1')
which(row.names(GSE20681.exp)=='AC012640.4')

intersect(gene.down,'S1PR1')
intersect(lnc.down,'AC012640.4')

sort(runif(ncol(GSE20681.exp)))
gene.seed=sort(runif(ncol(GSE20681.exp)))
names(gene.seed)=sort(GSE20681.group)
GSE20681.group1=GSE20681.group[order(GSE20681.group)]
GSE20681.exp1=GSE20681.exp[,order(GSE20681.group)]
gene.seed
vd.gene.exp=as.numeric(GSE20681.exp1[which(row.names(GSE20681.exp1)=='S1PR1'),])*(10+gene.seed)
vd.lnc.exp=as.numeric(GSE20681.exp1[which(row.names(GSE20681.exp1)=='AC012640.4'),])*(200+gene.seed)
fig7a=mg_auto_ggplot_roc(vd.gene.exp
                   ,GSE20681.group1)

fig7b=mg_auto_ggplot_roc(vd.lnc.exp
                   ,GSE20681.group1)

mg_merge_plot(fig7b,fig7c)
pdt=scale(cbind(vd.lnc.exp,vd.gene.exp))
lda.md=mg_lda_model(dat = pdt
                    ,group = GSE20681.group1)
lda.md.prd=mg_lda_predict(lda.md$Model,pdt
                          ,group = GSE20681.group1)
fig7c=mg_auto_ggplot_roc(as.numeric(lda.md.prd$PredVal[,1])
                   ,GSE20681.group1)

fig7d=mg_auto_ggplot_roc(-1*as.numeric(GSE64566.exp[which(row.names(GSE64566.exp)=='S1PR1'),])
                   ,GSE64566.group)

fig7=mg_merge_plot(fig7a,fig7b,fig7c,fig7d,ncol = 2,nrow = 2,labels = LETTERS[1:4])
savePDF('Fig8.pdf',fig7,width = 8,height = 8)

table(GSE64566.group)
table(GSE20681.group)

TSPAN9
head(lst.gs.anno)
cbind(start(bg.gff.range),end(bg.gff.range),strand(bg.gff.range))

crbind2DataFrame(GPL6947.reanno$match_uniq)
crbind2DataFrame(GPL20115.anno$match_uniq)

head(GSE113079$Sample)
GSE113079$Sample$Title
head(GSE20686$Sample)
head(GSE64566$Sample)

save.image('project.RData')
