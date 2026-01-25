## January 22, 2026

# Running easySFS to generate projections appropriate for my dataset
```
cd /scratch/a_monc/postdoc/xipho_project/GADMA
```
## Install easySFS
## https://github.com/isaacovercast/easySFS
```
conda create -n easySFS
conda activate easySFS
conda install -c conda-forge numpy pandas scipy -y
git clone https://github.com/isaacovercast/easySFS.git
cd easySFS
chmod 777 easySFS.py
./easySFS.py
```
## Create variables
```
VCF="/scratch/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_vcf_final.recode.vcf"
pop_file="/scratch/a_monc/postdoc/xipho_project/vcftools_filtering/xiph_pops.txt"
```
## estimate projections
```
./easySFS.py -i $VCF -p $pop_file -a -f --preview
```

## Results (not used)
Bel
(2, 10428)	(3, 15641)	(4, 19134)	(5, 21766)	(6, 23883)	(7, 25658)	(8, 27190)	(9, 28533)	(10, 29742)	(11, 30776)	(12, 31776)	(13, 32264)	(14, 33109)	(15, 31543)	(16, 32227)	(17, 25469)	(18, 25935)	(19, 12380)	(20, 12570)	

Tap
(2, 14537)	(3, 21806)	(4, 27141)	(5, 31509)	(6, 35280)	(7, 38636)	(8, 41687)	(9, 44470)	(10, 47087)	(11, 49231)	(12, 51539)	(13, 52259)	(14, 54292)	(15, 50889)	(16, 52575)	(17, 40050)	(18, 41185)	(19, 18653)	(20, 19107)	

Xin
(2, 11922)	(3, 17883)	(4, 22080)	(5, 25395)	(6, 28175)	(7, 30592)	(8, 32746)	(9, 34699)	(10, 36492)	(11, 38155)	(12, 39708)	(13, 41156)	(14, 42536)	(15, 43663)	(16, 44903)	(17, 45050)	(18, 46154)	(19, 43249)	(20, 44173)	(21, 35378)	(22, 36035)	(23, 21330)	(24, 21677)	(25, 7156)	(26, 7257)


For Tap, Xin, Bel, I select following projection: 14, 18, 14


# Running for the more strict VCF:
```
conda activate easySFS
VCF="/ddnA/work/a_monc/postdoc/xipho_project/vcftools_filtering/GADMA_STRICT_vcf_final.recode.vcf"
pop_file="/scratch/a_monc/postdoc/xipho_project/vcftools_filtering/xiph_pops.txt"
./easySFS.py -i $VCF -p $pop_file -a -f --preview
```
# Results for Strict VCF (VCF Dataset 5)
Bel
(2, 11678)	(3, 17516)	(4, 21273)	(5, 23990)	(6, 26095)	(7, 27804)	(8, 29237)	(9, 30469)	(10, 31548)	(11, 32508)	(12, 33373)	(13, 34159)	(14, 34881)	(15, 35549)	(16, 36170)	(17, 36751)	(18, 37297)	(19, 30262)	(20, 30650)	

Tap
(2, 15089)	(3, 22634)	(4, 27840)	(5, 31877)	(6, 35213)	(7, 38085)	(8, 40624)	(9, 42914)	(10, 45010)	(11, 46949)	(12, 48758)	(13, 50459)	(14, 52067)	(15, 53595)	(16, 55052)	(17, 56446)	(18, 57784)	(19, 45992)	(20, 46955)	

Xin
(2, 13149)	(3, 19724)	(4, 24101)	(5, 27378)	(6, 30006)	(7, 32210)	(8, 34115)	(9, 35801)	(10, 37316)	(11, 38698)	(12, 39969)	(13, 41150)	(14, 42253)	(15, 43290)	(16, 44270)	(17, 45200)	(18, 46084)	(19, 46929)	(20, 47738)	(21, 48515)	(22, 49262)	(23, 49982)	(24, 50678)	(25, 30657)	(26, 31037)


For Tap, Xin, Bel, I select following projection: 18, 24, 18