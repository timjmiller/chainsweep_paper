library("RefManageR")

bib = RefManageR::ReadBib("~/work/paired_tow_studies/chainsweep_paper/paper/sweep-paper.bib")
bib = RefManageR::ReadBib("~/work/bibliographies/fish.bib")

x = RefManageR::GetBibEntryWithDOI(c(
	"10.1139/f03-029",#Fryer et al. 2003
	"10.1016/j.fishres.2017.02.012",#Kotwicki et al. 2017
	"10.1139/F09-055", #Kotwicki et al. 2009
	"10.1023/A:1008838220001" #Millar and Fryer 1999
	), delete.file=F) 

RefManageR:::modify_url("https://doi.org/", path = "10.1023/A:1008838220001")
x = ReadPDFs("~/work/pubs")