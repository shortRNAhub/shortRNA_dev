# an <- readRDS("/mnt/IM/projects/software/ssc/annotation/anno.rds")
#
# anno <- data.frame(an@elementMetadata)
#
# head(anno)
#
# ans <- split(anno, anno$transcript_type)
#
# View(ans$tRNA)
#
# # Some of the genes IDs are numeric.
# # From here, I guess the tree is to be made on the transcript IDs.
# ## Either we keep tRNA-Ala-AGC-1-1 or tRNA_Ala_AGC_1_1 (- vs _)
# ## What about this?
# ans$tRNA[grep(pattern = "tRNA-Ala-AGC-1-1|tRNA_Ala_AGC_1_1", x = ans$tRNA$transcript_id),]
#
# ## Same with miRNAs
# ans$miRNA[grep(pattern = "mmu-let-7g-5p|mmu_let_7g_5p", x = ans$miRNA$transcript_id),]
#
# ## How to deal with this while making the tree
# ans$miRNA_hairpin[grep(pattern = "let-7g", x = ans$miRNA_hairpin$transcript_id),]
#
#
# ## How to put these in tree?
# head(ans$pseudo_tRNA)
#
# ## Also these?
# head(ans$putative_tRNA)
#
# ## and these
# head(ans$Mt_tRNA)
#
# ### Similar issues with miRNA
#
# # rRNAs are available only as fasta, how should I add the sequences to the fasta file?
# # you mentioned to use a GTF, but we are not using a GTF file at all.
#
#
# ov <- readRDS("../../shortRNA_reports/uterosome/03_tse/overlap.rds")
