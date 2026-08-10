cellphoneDB_Dotplot <- function(pvals.data,
                                means.data,
				key,
                                target.cells_1,
                                target.cells_2 = NA,
                                gene_a = NA,
                                gene_b = NA,
                                p.cutoff = 0.05,
                                xlab=NULL,
                                ylab=NULL,
                                title=NULL,
                                text.size = 8,
                                text.angle = 90,
                                text.hjust = 0,
                                legend.position = "right",
                                ...
                                ){
  library(stringr)
  library(tidyverse)
  mytheme <- theme(plot.title = element_text(size = text.size+2,color="black",hjust = 0.5),
                   axis.ticks = element_line(color = "black"),
                   axis.title = element_text(size = text.size,color ="black"),
                   axis.text = element_text(size=text.size,color = "black"),
                   axis.text.x = element_text(angle = text.angle, hjust = text.hjust ), #,vjust = 0.5
                   #axis.line = element_line(color = "black"),
                   #axis.ticks = element_line(color = "black"),
                   #panel.grid.minor.y = element_blank(),
                   #panel.grid.minor.x = element_blank(),
                   panel.grid=element_blank(), # 去网格线
                   legend.position = legend.position,
                   legend.text = element_text(size= text.size),
                   legend.title= element_text(size= text.size),
                   panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
                   strip.background = element_rect(color="black",size= 1, linetype="solid") # fill="#FC4E07",
  )
  colnames(pvals.data) = str_replace(string = colnames(pvals.data),
                                  pattern = "\\.",
                                  replacement = "\\_")

  colnames(means.data) = str_replace(string = colnames(means.data),
                                pattern = "\\.",
                                replacement = "\\_")

  if(is.na(target.cells_2[1])){
    kp = grepl(pattern = target.cells_1[1], colnames(pvals.data))
    for (i in target.cells_1[-1]) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp =  kp|tpm
    }
  }else{
    kp_1 = grepl(pattern = target.cells_1[1], colnames(pvals.data))
    for (i in target.cells_1) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp_1 =  kp_1|tpm
    }

    kp_2 = grepl(pattern = target.cells_2[1], colnames(pvals.data))
    for (i in target.cells_2) {
      tpm = grepl(pattern = i, colnames(pvals.data))
      kp_2 =  kp_2|tpm
    }

    kp = kp_1 & kp_2
  }

  pos = (1:ncol(pvals.data))[kp]
  #choose_pvalues <- pvals.data[,c(c(1,5,6,8,9),pos)]
  choose_pvalues <- pvals.data[,c(c(1,2,5,6,8,9),pos)]
  #choose_means <- means.data[,c(c(1,5,6,8,9),pos)]
  choose_means <- means.data[,c(c(1,2,5,6,8,9),pos)]

  logi <- apply(choose_pvalues[,5:ncol(choose_pvalues)]<p.cutoff, 1, sum)
  # 只保留具有细胞特异性的一些相互作用对
  choose_pvalues <- choose_pvalues[logi>=1,]

  # 去掉空值
  logi1 <- choose_pvalues$gene_a != ""
  logi2 <- choose_pvalues$gene_b != ""
  logi <- logi1 & logi2
  choose_pvalues <- choose_pvalues[logi,]

  # 同样的条件保留choose_means
  choose_means <- choose_means[choose_means$id_cp_interaction %in%
                                 choose_pvalues$id_cp_interaction,]

  # 将choose_pvalues和choose_means数据宽转长
  meansdf <- choose_means %>% reshape2::melt()
  #meansdf <- data.frame(interacting_pair = paste0(meansdf$gene_a,"_",meansdf$gene_b),
  meansdf <- data.frame(interacting_pair = meansdf$interacting_pair,
                        CC = meansdf$variable,
                        means = meansdf$value)
  print(head(meansdf))
  pvalsdf <- choose_pvalues %>% reshape2::melt()
  #pvalsdf <- data.frame(interacting_pair = paste0(pvalsdf$gene_a,"_",pvalsdf$gene_b),
  pvalsdf <- data.frame(interacting_pair = meansdf$interacting_pair,
                        CC = pvalsdf$variable,
                        pvals = pvalsdf$value)
  print(head(pvalsdf))

  # 合并p值和mean文件
  pvalsdf$joinlab<- paste0(pvalsdf$interacting_pair,"_",pvalsdf$CC)
  meansdf$joinlab<- paste0(meansdf$interacting_pair,"_",meansdf$CC)
  pldf <- merge(pvalsdf,meansdf,by = "joinlab")
  colnames(pldf) =   c("joinlab","interacting_pair","Celltypes",
                       "pvals", "interacting_pair.means",
                       "CC.means", "means")

  if(!is.na(gene_a[1])){
    gene_a.df = str_split(pldf$interacting_pair,pattern = "\\_",simplify = T)[,1]
    pldf = pldf[gene_a.df %in%gene_a,]
  }

  if(!is.na(gene_b[1])){
    gene_b.df = str_split(pldf$interacting_pair,pattern = "\\_",simplify = T)[,2]
    pldf = pldf[gene_b.df %in%gene_b,]
  }
  # dotplot可视化
  summary((filter(pldf,means >0))$means)
  head(pldf)

# key <- "HSC"
 pldf$Celltypes <-factor(pldf$Celltypes,levels = unique(c(
                               # 先选择 'HSC' 在 '|' 前的部分
                               sort(pldf$Celltypes[grepl(paste0("^", key, "\\|"), pldf$Celltypes)]),
                               # 然后选择 'HSC' 在 '|' 后的部分
                               sort(pldf$Celltypes[grepl(paste0("\\|", key), pldf$Celltypes)])
                             )))
  #pcc =  pldf%>% filter(means >0) %>%

  plotdata=pldf%>%select(interacting_pair,Celltypes,pvals,means)
  write.table(plotdata,paste0("cellphoneDB_Dotplot.",key,".plotdata.xls"),row.names=F,quote=F,sep="\t")
  pcc =  pldf %>%
    ggplot(aes(Celltypes,interacting_pair) )+
    geom_point(aes(color=log2(means+1),size=-log10(pvals+0.0001))) + # add log2
    scale_size_continuous(range = c(0.3,3),name="-log10(pvals)")+
#    scale_color_gradient2(high="red",mid = "yellow",
#                          low ="darkblue",
#                          midpoint = 1,name="Means")+
#                          name="log2(Means+1)")+
    theme_bw()+
    #labs(y= ylab,x= xlab,title = title) +
    labs(y= "",x= "",title = "") +
    # scale_color_manual(values = rainbow(100))+
     scale_colour_gradientn(colors=c("#3A5978","#5689B3","#F6B31D","#DA2328","#8F342D"))+
    mytheme
  return(pcc)
}
