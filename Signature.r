# nature medicine stem cell signature 
gene <- c("VSNL1","AKR1B1","NOTCH4","TMEM237","SAMD5","PKD2",
"NAP1L1","PTTG1","CDK6","CDCA7","ACSL4","HELLS",
"IKBIP","PLTP","TMEM201","CACHD1","ILF3","DNMT1",
"USP31","FAM216A","SLC41A1","PFKM","KANK1","SUPT16H",
"ADCY3","FGD1","PTPN14","C20orf27","LGR6","SLC16A7",
"JAM3","FBL","NASP","RANBP1","PRNP","DSE","GPX7",
"KDELC1","FCHSD2","SLCO3A1","CCNB1IP1",
"LOC284023","NOL9","NKRF","NUP107","RCC2","ARHGAP25","DDX46","TCOF1","GMPS")
name="stem_cell"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"


# hallmark EMT 
gene <- c("ABI3BP","ACTA2","ADAM12","ANPEP",
"APLP1","AREG","BASP1","BDNF",
"BGN","BMP1","CADM1","CALD1",
"CALU","CAP2","CAPG","CD44","CD59",
"CDH11","CDH2","CDH6","COL11A1","COL12A1",
"COL16A1","COL1A1","COL1A2","COL3A1","COL4A1",
"COL4A2","COL5A1","COL5A2","COL5A3","COL6A2",
"COL6A3","COL7A1","COL8A2","COMP","COPA",
"CRLF1","CTGF","CTHRC1","CXCL1","CXCL12",
"CXCL6","CYR61","DAB2","DCN","DKK1",
"DPYSL3","DST","ECM1","ECM2","EDIL3",
"EFEMP2","ELN","EMP3","ENO2","FAP",
"FAS","FBLN1","FBLN2","FBLN5","FBN1",
"FBN2","FERMT2","FGF2","FLNA","FMOD",
"FN1","FOXC2","FSTL1","FSTL3","FUCA1",
"FZD8","GADD45A","GADD45B","GAS1","GEM",
"GJA1","GLIPR1","GLT25D1","GPC1","GPX7",
"GREM1","HTRA1","ID2","IGFBP2","IGFBP3",
"IGFBP4","IL15","IL32","IL6","IL8","INHBA",
"ITGA2","ITGA5","ITGAV","ITGB1","ITGB3",
"ITGB5","JUN","LAMA1","LAMA2","LAMA3","LAMC1",
"LAMC2","LEPRE1","LGALS1","LOX","LOXL1",
"LOXL2","LRP1","LRRC15","LUM","MAGEE1",
"MATN2","MATN3","MCM7","MEST","MFAP5","MGP",
"MMP1","MMP14","MMP2","MMP3","MSX1",
"MXRA5","MYL9","MYLK","NID2","NNMT",
"NOTCH2","NT5E","NTM","OXTR","PCOLCE",
"PCOLCE2","PDGFRB","PDLIM4","PFN2","PLAUR",
"PLOD1","PLOD2","PLOD3","PMEPA1","PMP22",
"POSTN","PPIB","PRRX1","PRSS2","PTHLH","PTX3",
"PVR","QSOX1","RGS4","RHOB","SAT1","SCG2",
"SDC1","SDC4","SERPINE1","SERPINE2",
"SERPINH1","SFRP1","SFRP4","SGCB","SGCD",
"SGCG","SLC6A8","SLIT2","SLIT3","SNAI2",
"SNTB1","SPARC","SPOCK1","SPP1","TAGLN",
"TFPI2","TGFB1","TGFBI","TGFBR3","TGM2",
"THBS1","THBS2","THY1","TIMP1","TIMP3","TNC",
"TNFAIP3","TNFRSF11B","TNFRSF12A","TPM1",
"TPM2","TPM4","VCAM1","VCAN","VEGFA","VEGFC","VIM","WIPF1","WNT5A")
name="Hallmark_EMT"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"







# Nature medicine CELLULAR_MORPHOGENESIS_DURING_DIFFERENTIATION
gene <- c("ALS2","AMIGO1","APOE","BAI1","BAIAP2",
"CDK5R1","CEP290","CNTN4","CYFIP1",
"DPYSL5","FEZ1","FEZ2","GLI2","KAL1",
"KLK8","LRRC4C","MAP1S","MAPT","NRL",
"NRP1","NRP2","NRXN1","NRXN3","NTNG1",
"NTNG2","OPHN1","OTX2","PARD3","PARD6B",
"PAX2","POU4F1","ROBO1","ROBO2","RORB",
"RTN4","RTN4RL1","RTN4RL2","S100B",
"SEMA3B","SEMA4F","SHH","SIAH1","SLIT1",
"SLIT2","SPON2","THY1","UBB","UNC5C","YWHAH")
name="MORPHOGENESIS"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"




# Nature Medicine AEP
gene <- c("SFTPC","TECTA","POPDC3","PGC","RIMS4",
"AGTR2","SFTPA1","LHFPL3","NAPSA","PTPRD",
"PDPN","LAMP3","SORCS2","AARD","RP11-22C11.2",
"PLA2G1B","MAPT","BCL2L15","MYBPHL","GGTLC1",
"RP3-467L1.6","ABCA3","SEZ6L2","WFDC12","SFTPD",
"HOXD-AS1","ODAM","NBPF3","BMP2","GPR133",
"DPP4","PARM1","C4BPA","SDR16C5","KRT121P","CLDN18",
"CEACAM5","VWA2","SLC34A2","NPR1","SYTL5",
"SMPDL3B","HOPX","RGS16","ABCC6P2","AC090616.2",
"SPNS2","CHIA","EPHA1","SHE","HSD17B6","S100A2",
"TACC2","KRT16P1","PEBP4","CRTAC1","DUOX1","LAD1",
"AC008268.1","NTN4","GKN2","SFTA2","ARSJ","ABCC13",
"FMO5","VEPH1","SEC14L6","AXIN2","C1orf210","SGPP2",
"MAOA","PDZK1IP1","ITGB6","GPRIN2","GSTA4","PLLP",
"DUOXA1","SFN","C14orf132","FLRT3","SLC22A31",
"MFSD2A","C1R","AC034193.5","CREB3L1","MMP28",
"SMAGP","CDH3","SLC22A3","KANK2","SFTA3","COL12A1",
"GPR116","SLC6A14","PLS3","RAB17","PTK7","SEPP1","YAP1",
"AMN","AK1","CCT8P1","LIMCH1","SLC39A8","SEMA3B","CHI3L1",
"SUSD2","DOK4","STXBP1","AQP4","TMPRSS2","S100A14","SCEL",
"BEX2","MUC21","LRRC36","SELENBP1","SULT1A2","LAPTM4B",
"TMEM125","CFI","HNF1B","CFTR","ADAMTS1","LMO3","MUC1",
"SLAIN1","SLC44A3","MYO1B","TTC39A","PTPN13","PEG10",
"CADM1","GRB7","SNX25","NPNT","RAI2","SLC46A2","PRSS8",
"PPP1R9A","PID1","ALS2CL","EMP2","PON3","SFTA1P","CA12",
"CTSE","C4A","LPCAT1","RNF128","VSIG2","FERMT2","NEO1",
"ASS1","EFNB2","ID1","GPC4","CAV1","VSIG10","MECOM",
"PCP4L1","LURAP1L","SMARCA1","MALL","CRB3","PRR15L",
"ELF5","NCKAP1","HIP1R","DSG2","NFIX","FGFR2","ESRP2",
"RND1","KCNS3","PARD3","EPCAM","RNF180","OCLN","SLC25A23",
"CKB","RHOBTB2","CIT","SLC25A4","ARHGAP29","SDC1","LRP5","KRT7",
"COBLL1","NGFRAP1","C8orf4","ATP13A4","MAL2","MBIP","DAPK2",
"TM7SF2","IL20RA","RNF43","PRKCZ","SH3D19","EFNA1","FAM84B",
"AP1M2","CYBRD1","C16orf89","F3","AHCYL2","TMEM97","USP54",
"MSLN","CDH1","MAGI3","LSR","FXYD3","NAB2","GTF2IRD1",
"RAB25","PARD6B","NR1D1","IGFBP4","CGN","GPRC5C","C3","ARMC9",
"OSMR","MMP7","PPARGC1A","ATP2C2","EPS8","CXCL17","PKP3",
"MYO5C","FNBP1L","CNN3","LMO7","C1orf116","MTUS1","DDAH1",
"LGALSL","TMEM30B","SH2D4A","EVPL","MET","ANKRD29","GPR125",
"TM4SF1","ESRP1","LNX2","GPRC5A","SCNN1A","ALOX15B","LIPH",
"CEP70","CYP2B7P1","FADS2","TEF","LRIG3","AMMECR1","PMM1",
"MYH14","RASEF","PCDH1","NNMT","ANKRD65","SDHAP3","RAB40C",
"PDIA5","RNF145","NR3C2","CHI3L2","TTC3P1","RAB3IP","MAP3K13",
"GPD1L","EFEMP1","PXMP4","MTND4P12","FAM20A","TANC1","SQLE",
"TBCEL","NEBL","ERBB3","PPA1","ST3GAL5","GATA6","INADL","CTTN",
"FOLR1","KRT8","BNIP3","DHCR24","FGFR1","TFPI","DLC1","LDLR",
"ABLIM1","CYR61","KLF5","PON2","AMOTL2","MYC","CDS1","ELF3",
"CLIC5","FZD6","CDC42EP1","IFT57","KRT18","PLS1","TSPAN12","CTNNBIP1")

name="AEP"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"















# PNAS stemnness signature 
# https://www.pnas.org/content/pnas/116/18/9020.full.pdf 
gene <- c("DNMT3B","PFAS","XRCC5","HAUS6","TET1","IGF2BP1","PLAA",
"TEX10","MSH6","DLGAP5","SKIV2L2","SOHLH2","RRAS2","PAICS","CPSF3",
"LIN28B","IPO5","BMPR1A","ZNF788","ASCC3","FANCB","HMGA2","TRIM24",
"ORC1","HDAC2","HESX1","INHBE","MIS18A","DCUN1D5","MRPL3","CENPH",
"MYCN","HAUS1","GDF3","TBCE","RIOK2","BCKDHB","RAD1","NREP","ADH5",
"PLRG1","ROR1","RAB3B","DIAPH3","GNL2","FGF2","NMNAT2","KIF20A",
"CENPI","DDX1","XXYLT1","GPR176","BBS9","C14orf166","BOD1","CDC123",
"SNRPD3","FAM118B","DPH3","EIF2B3","RPF2","APLP1","DACT1","PDHB",
"C14orf119","DTD1","SAMM50","CCL26","MED20","UTP6","RARS2","ARMCX2",
"RARS","MTHFD2","DHX15","HTR7","MTHFD1L","ARMC9","XPOT","IARS","HDX",
"ACTRT3","ERCC2","TBC1D16","GARS","KIF7","UBE2K","SLC25A3","ICMT",
"UGGT2","ATP11C","SLC24A1","EIF2AK4","GPX8","ALX1","OSTC","TRPC4",
"HAS2","FZD2","TRNT1","MMADHC","SNX8","CDH6","HAT1","SEC11A","DIMT1","TM2D2","FST","GBE1")

name="PNAS_stem"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,name,out=paste0(respath,name,".mod")) # make a mod file 











# MSGDB signature 
gene <- c("AKT1","BAP1","MRE11","NF2","NSMCE3","PDGFB","PIK3CA","RNF168","SMARCB1","SMARCE1","SMO","SUFU","TERT","TRAF7")

name="Radio_sens"
respath="/public/workspace/lily/MOD_file/NatureMedicine/"


source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,name,out=paste0(respath,name,".mod")) # make a mod file 





gene.f <- c("LGALS9","ARG1","CD47","IDO1","VEGFA","IL10","TGFB1","PVR")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene.f,"immunsupp",out="/public/workspace/lily/MOD_file/immunsupp.mod") # make a mod file 












# 2021-9-28
# Brain gene and Lung gene signature 
# come form https://bioinfo.uth.edu/TissGDB/download_dir/g_ctype_tissue.txt?csrt=17039100935693777611 
# Tissue-specific genes (TissGenes) with cancer type and tissue information

tmp <- read.table("~/tmp/g_ctype_tissue.txt")
unique(as.vector(tmp[which(tmp$V2%in%c("GBM","LGG")),1]))
unique(as.vector(tmp[which(tmp$V2%in%c("LUAD")),1]))







# 2021-10-14
# another radio signature 
# https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-13-348/tables/1
ACTN1	Down (0.42)
ANXA2	Down (0.36)
ANXA5	Down (0.42)
ARHGDIB	Up (−0.49)
CAPNS1	Down (0.48)
CBR1	Down (0.41)
CCND1	Down (0.54)
CD63	Down (0.51)
CORO1A	Up (−0.46)
CXCR4	Up (−0.46)
DAG1	Down (0.60)
EMP2	Down (0.41)
HCLS1	Up (−0.58)
HTRA1	Down (0.52)
ITGB5	Down (0.47)
LAPTM5	Up (−0.50)
LRMP	Up (−0.49)
MYB	Up (−0.59)
PFN2	Down (0.61)
PIR	Down (0.43)
PKM2	Down (0.44)
PTMS	Down (0.48)
PTPRC	Up (−0.55)
PTPRCAP	Up (−0.49)
PYGB	Down (0.35)
RAB13	Down (0.43)
RALB	Down (0.47)
SCRN1	Down (0.40)
SQSTM1	Down (0.48)
TWF1	Down (0.43)
WAS	Up (−0.60)

tmp <- read.table("/public/workspace/lily/tmp/tmp.txt",sep="\t")
gene <- as.numeric(gsub("−","-",(gsub("\\)","",sapply(strsplit(as.character(tmp$V2),"\\("),function(x){x[[2]]})))))
names(gene) <- as.character(tmp$V1)
gene <- ifelse(gene<0,1,(-1))

gene <- gene[which(gene>0)]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"Radio_sen_BMC",out="/public/workspace/lily/MOD_file/Radio_sen_BMC.mod") # make a mod file 











ATP4A	-2.03
COLQ	-0.42
AMDHD1	0.42
PER1	-0.3
KPNA2	0.34
TAP2	0.39
MSMO1	0.29
BTG3	0.12
TUBB2A	0.64
CMPK2	-0.04
NEK2	0.91
ZFYVE1	-0.01
ATP5H	-0.43
LOC283398	0.07
E2F8	-0.63
GGA2	-0.01
DDIT3	0.97
TXNIP	-0.33
LINS	-1.85
ZNF787	0.38
PLEK2	-0.06
PLA2G4A	0.6
SARS2	-0.4
MMD	0.58
BBS9	-0.34
FGFR1OP	-0.33
ANO3	-0.33
NUDT15	0.47
COX7B	1.2
H2AFJ	5.5
PYGB	-1.12
USP7	-0.58
PHEX	0.39
FAM114A2	-0.43

tmp <- read.table("/public/workspace/lily/tmp/tmp.txt",sep="\t")
gene <- as.numeric(tmp$V2)
names(gene) <- as.character(tmp$V1)
gene <- ifelse(gene<0,1,(-1))

# gene <- gene[which(gene>0)]
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"Radio_sen_CCR",out="/public/workspace/lily/MOD_file/Radio_sen_CCR.mod") # make a mod file 









#======================================================================================
# 2021-10-14
# make this 
# use this https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3038688/
gene <- c("AR","PRKCA","RELA","SUMO1","CDK1","HDAC1","IRF1")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"Radio_sen_test",out="/public/workspace/lily/MOD_file/Radio_sen_test.mod") # make a mod file 

sapply(gene,function(x){
    c(cor.test(mod$Brain_gene_norm,as.numeric(dat.BM[x,]),method="spearman")$estimate,
    cor.test(mod$Lung_gene_norm,as.numeric(dat.BM[x,]),method="spearman")$estimate,
    cor.test(mod$BMS_update_norm,as.numeric(dat.BM[x,]),method="spearman")$estimate)
    })

# gene <- c("KIAA0586","C14orf135","TRMT5","GLMN","MUDENG","ZMYM6","SFRS5")
# source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
# mod.generate(gene,"Radio_sen_test",out="/public/workspace/lily/MOD_file/Radio_sen_test.mod") # make a mod file 



gene <- c("RPIA","RBBP4","RGS19")
source('/public/workspace/lily/software/ssGSEA/ssgseaMOD.r')
mod.generate(gene,"Radio_sen_test",out="/public/workspace/lily/MOD_file/Radio_sen_test.mod") # make a mod file 



