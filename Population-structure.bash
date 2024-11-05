### 系统发育树
## step1: vcf2phylip
python /data01/wangyf/software/vcf2phylip-2.8/vcf2phylip.py --input /data01/wangyf/project2/CyprinusCarpio/15.pop/raw.indel5.biSNP.QUAL30.QD3.FS20.MQ55.SOR3.MQRS-5.RPRS-5.PASS.GQ10.popmiss.maxmiss0.15.AF0.05.10-3ClusterFilter.vcf.gz --output-folder . --output-prefix snp

## step 2: phylip to iqtree
# -s: input file
# -st: set the category of the input sequence
# -m: choose substitution model

conda activate iqtree
iqtree -s snp.min4.phy.varsites.phy -st DNA -m MFP+ASC -B 1000 --bnni -T 20
# -m MFP               Extended model selection followed by tree inference
# -m ...+ASC           Ascertainment bias correction
# it will report an error and output a file *."varsites.phy", then let this file as input file and run again
# 结果文件 .contree可以用Figtree或者iTol网站打开美化

############################################################################################################################################################

##  PCA
conda activate vcftools
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
conda activate r4.3
R
library(ggplot2)

setwd("/home/wangmin/project/01_Cyprinus_carpio/00_sanxia/09_PCA")
data <- read.table("pca.info.txt",sep = "\t",header = T)
p <- ggplot(data, aes(x = PC1, y = PC2,color = data[,1]))+
  geom_point(size=3, shape = 16) +
  theme_bw()+
  labs(x = "PC1(36.6%)", y = "PC2(22.6%)")+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_blank(),legend.text=element_text(size=12),title=element_text(size=12),axis.text = element_text(size = 10),text=element_text(size=16,  family="serif"))
  
######################################################################################################################

### structure
source activate plink

# .vcf.gz to plink format, get file .map .nosex .ped
plink --vcf Irtysh.vcf.gz --recode --threads 10 --out plink --const-fid --allow-extra-chr && \

# add SNP number
x=0; awk -v x="$x" '{x+=1}{print $1"\tSNP"x"\t"$3"\t"$4}' plink.map > admix.map

# plink filter and get the .bed file (already)
plink --allow-extra-chr --noweb --file plink --geno 0.05 --maf 0.05 --hwe 0.0001 --make-bed --out admixture --chr-set 1

# use awk to rewrite $1 in *.bim as 整数
awk '{$1 = substr($1, 4)} 1' admixture.bim > newadmixture.bim # delete "NC_"
awk '{$1=substr($1,1,6) "" substr($1,10)} 1' newadmixture.bim > new.admixture.bim # delete ".1"
mv new.admixture.bim admixture.bim
rm newadmixture.bim

# find a suitable k-value
source activate admixture
for i in {1..10};do admixture --cv admixture.bed $i -j10 |tee log${i}.out;done

# choose smallest k-value
grep -h CV log*.out |sort -nk4  -t ' ' > cv_error.txt
