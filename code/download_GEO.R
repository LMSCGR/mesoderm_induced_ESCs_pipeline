library(RCurl)
url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198730/suppl/"
filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
ifelse(Sys.info()["sysname"]=="Windows",
       filenames <- strsplit(filenames, "\r\n"),filenames <- strsplit(filenames, "\n"))


filenames = unlist(filenames)
filenames

file_list<-c("GSE198730_aPSM_scRNA_rep1_barcodes.tsv.gz","GSE198730_aPSM_scRNA_rep1_features.tsv.gz","GSE198730_aPSM_scRNA_rep1_matrix.mtx.gz")
get_files<-intersect(filenames,file_list)

dest_dir<-paste(getwd(), "/filtered_feature_bc_matrix/",sep="")
dir.create(dest_dir, showWarnings = FALSE)

for (filename in get_files) {
  download.file(paste(url, filename, sep = ""), paste(getwd(), "/filtered_feature_bc_matrix/",filename,sep = ""), mode = "wb")
}

file.rename(paste(dest_dir,"GSE198730_aPSM_scRNA_rep1_barcodes.tsv.gz",sep=""),paste(dest_dir,"barcodes.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_aPSM_scRNA_rep1_features.tsv.gz",sep=""),paste(dest_dir,"features.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_aPSM_scRNA_rep1_matrix.mtx.gz",sep=""),paste(dest_dir,"matrix.mtx.gz",sep=""))


file_list<-c("GSE198730_aPSM_scATAC_rep1_filtered_peak_bc_matrix.h5","GSE198730_aPSM_scATAC_rep1_singlecell.csv.gz","GSE198730_aPSM_scATAC_rep1_fragments.tsv.gz")
get_files<-intersect(filenames,file_list)

dest_dir<-paste(getwd(), "/aPSM_scATAC/",sep="")
dir.create(dest_dir, showWarnings = FALSE)

for (filename in get_files) {
  download.file(paste(url, filename, sep = ""), paste(getwd(), "/aPSM_scATAC/",filename,sep = ""), mode = "wb")
}

file.rename(paste(dest_dir,"GSE198730_aPSM_scATAC_rep1_filtered_peak_bc_matrix.h5",sep=""),paste(dest_dir,"filtered_peak_bc_matrix.h5",sep=""))
file.rename(paste(dest_dir,"GSE198730_aPSM_scATAC_rep1_fragments.tsv.gz",sep=""),paste(dest_dir,"fragments.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_aPSM_scATAC_rep1_singlecell.csv.gz",sep=""),paste(dest_dir,"singlecell.csv.gz",sep=""))

url = "https://github.com/LMSCGR/mesoderm_induced_ESCs_pipeline/blob/main/code/"
filename="aPSM_fragments.tsv.gz.tbi"
download.file(paste(url, filename, sep = ""), paste(getwd(), "/aPSM_scATAC/",filename,sep = ""), mode = "wb")
file.rename(paste(dest_dir,"aPSM_fragments.tsv.gz.tbi",sep=""),paste(dest_dir,"fragments.tsv.gz.tbi",sep=""))

##integration
url = "https://raw.githubusercontent.com/LMSCGR/mesoderm_induced_ESCs_pipeline/main/code/naive_instructed_esc.csv"
filename="naive_instructed_esc.csv"
download.file(url, paste(getwd(), "./",filename,sep = ""), mode = "wb")

url = "https://raw.githubusercontent.com/LMSCGR/mesoderm_induced_ESCs_pipeline/main/code/aPSM_f.txt"
filename="aPSM_f.txt"
download.file(url, paste(getwd(), "./",filename,sep = ""), mode = "wb")

url = "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE198nnn/GSE198730/suppl/"
filenames = getURL(url, ftp.use.epsv = FALSE, dirlistonly = TRUE)
ifelse(Sys.info()["sysname"]=="Windows",
       filenames <- strsplit(filenames, "\r\n"),filenames <- strsplit(filenames, "\n"))


filenames = unlist(filenames)
filenames

file_list<-c("GSE198730_Naive_ESC_scRNA_rep1_barcodes.tsv.gz","GSE198730_Naive_ESC_scRNA_rep1_features.tsv.gz","GSE198730_Naive_ESC_scRNA_rep1_matrix.mtx.gz")
get_files<-intersect(filenames,file_list)
dest_dir<-""
dest_dir<-paste(getwd(), "/naive_scRNA",sep="")
dir.create(dest_dir, showWarnings = FALSE)
dest_dir<-paste(dest_dir, "/filtered_feature_bc_matrix/",sep="")
dir.create(dest_dir, showWarnings = FALSE)

for (filename in get_files) {
  download.file(paste(url, filename, sep = ""), paste(dest_dir,filename,sep = ""), mode = "wb")
}

file.rename(paste(dest_dir,"GSE198730_Naive_ESC_scRNA_rep1_barcodes.tsv.gz",sep=""),paste(dest_dir,"barcodes.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_Naive_ESC_scRNA_rep1_features.tsv.gz",sep=""),paste(dest_dir,"features.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_Naive_ESC_scRNA_rep1_matrix.mtx.gz",sep=""),paste(dest_dir,"matrix.mtx.gz",sep=""))

file_list<-c("GSE198730_Instructed_ESC_scRNA_rep1_barcodes.tsv.gz","GSE198730_Instructed_ESC_scRNA_rep1_features.tsv.gz","GSE198730_Instructed_ESC_scRNA_rep1_matrix.mtx.gz")
get_files<-intersect(filenames,file_list)
dest_dir<-""
dest_dir<-paste(getwd(), "/instructed_scRNA",sep="")
dir.create(dest_dir, showWarnings = FALSE)
dest_dir<-paste(dest_dir, "/filtered_feature_bc_matrix/",sep="")
dir.create(dest_dir, showWarnings = FALSE)

for (filename in get_files) {
  download.file(paste(url, filename, sep = ""), paste(dest_dir,filename,sep = ""), mode = "wb")
}

file.rename(paste(dest_dir,"GSE198730_instructed_ESC_scRNA_rep1_barcodes.tsv.gz",sep=""),paste(dest_dir,"barcodes.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_instructed_ESC_scRNA_rep1_features.tsv.gz",sep=""),paste(dest_dir,"features.tsv.gz",sep=""))
file.rename(paste(dest_dir,"GSE198730_instructed_ESC_scRNA_rep1_matrix.mtx.gz",sep=""),paste(dest_dir,"matrix.mtx.gz",sep=""))
