##hard-filter之后进行PCA
##min to xiaotong to me，服务器中有更新的R作图脚本

plink --vcf ../06_snp_filtered/GQ15/merge.51.snp.raw_vcftools_gatk_pass_sedname.vcf --recode --out merge.51.plink --const-fid --allow-extra-chr

plink --allow-extra-chr --file merge.51.plink --noweb --make-bed --out merge.51.plink

plink --allow-extra-chr --threads 20 -bfile merge.51.plink --pca 100 --out merge.51.plink

#准备作图文件(pop individual PC1 PC2 PC3 ...)
cat carp.plink.eigenvec | awk '{print $2}'  > pop.txt  #群体名称文件，具体怎么弄看文件
cat carp.plink.eigenvec | awk '{print $2"\t"$3"\t"$4"\t"$5}' > pca.txt  #提取前3个主成分
paste pop.txt pca.txt | tr " " "\t" > pca.info.txt  #合并
sed -i "1i/pop\tindividual\tPC1\tPC2\tPC3" pca.info.txt  #添加抬头信息


#R计算每个PC的解释率，修改下文的PC1、2、3的大小
eigval <- read.table("plink.eigenval", header = F)
percentage <- eigval$V1/sum(eigval$V1)*100
print(percentage)

#R出图
R

library(ggplot2)

setwd("/home/wangmin/project/01_Cyprinus_carpio/00_sanxia/09_PCA")

data <- read.table("pca.info.txt",sep = "\t",header = T)

p <- ggplot(data, aes(x = PC1, y = PC2,color = data[,1]))+
  geom_point(size=3, shape = 16) +
  theme_bw()+
  labs(x = "PC1(36.6%)", y = "PC2(22.6%)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_blank(),legend.text=element_text(size=12),title=element_text(size=12),axis.text = element_text(size = 10),text=element_text(size=16,  family="serif"))
  