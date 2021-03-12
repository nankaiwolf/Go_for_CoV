# 这些患者的RNA比对中发现有SARS-CoV的物种注释。
# 为了确认这些序列的可靠性，就写了这段代码回溯找寻原始序列进行分析。
# 读入mpa格式的物种表格
tax.rna.count <- read.csv("tax.rna.mpa",
                          sep = "\t", header = F)
# 提取其中Viruses的序列数量
virus.rna.count <- tax.rna.count[grep("d__Viruses", tax.rna.count[ , 1], invert = F), ]
# 读取序列注释表格
reads.rna.table <- read.csv("tax.rna.kra",
                            sep = "\t", header = F)
# 保留已经匹配（Classified）的序列
classified.reads.rna <- reads.rna.table[grep("C", reads.rna.table[ , 1], invert = F), ]
# 在kraken2注释里SARS-CoV2的物种编号是2697049
# 这里大家先整体看看各个物种的数量情况
summary(as.factor(classified.reads.rna[ , 3]))

# 读入含有疑似序列的fastq文件，记得这一步千万要去除读取注释，因为有的质控行可能以字符“#”开头。
reads.out <- read.table("seq.rna.fastq", sep = "\n", comment.char = "")
# 把数据的四列分离出来。
reads.data <- reads.out[as.numeric(rownames(reads.out))%%4 == 1, ]
reads.data <- cbind(reads.data, reads.out[as.numeric(rownames(reads.out))%%4 == 2, ])
reads.data <- cbind(reads.data, reads.out[as.numeric(rownames(reads.out))%%4 == 3, ])
reads.data <- cbind(reads.data, reads.out[as.numeric(rownames(reads.out))%%4 == 0, ])

# 从reads数据中选出符合注释要求的序列
reads.selected <- 
  reads.data[reads.data[ , 1] %in% paste0('@', classified.reads.rna[classified.reads.rna$V3 == "2697049", 'V2']), ]
# 将这些序列写入文档
write.table(reads.selected, file = "COVID.suspected.fastq", sep = '\n',
            row.names = F, col.names = F, quote = F)
# 最后可以用Linux命令转换为fasta格式
# awk '{if(NR%4 == 1){print ">" substr($0, 2)}}{if(NR%4 == 2){print}}' COVID.suspected.fastq >COVID.suspected.fasta
