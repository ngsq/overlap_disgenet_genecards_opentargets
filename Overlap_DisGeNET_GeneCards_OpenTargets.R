library(rentrez)
library(dplyr)
library(readr)
library(scales)
library(ggplot2)
library(VennDiagram)
getwd()
setwd('D:/Downloads')
# Read downloaded OpenTarget files
as_data <- read_tsv("AS_OpenTarget.tsv")
pso_data <- read_tsv("PsO_OpenTarget.tsv")

# Filter genes with Relevance score >= 0.05
as_filtered <- as_data[as_data$globalScore >= 0.05, ]
pso_filtered <- pso_data[pso_data$globalScore >= 0.05, ]

# Get gene symbols
as_genes <- unique(as_filtered$symbol)
pso_genes <- unique(pso_filtered$symbol)

# Find shared genes for OpenTarget
shared_genes <- intersect(as_genes, pso_genes)
print(shared_genes)
print(length(shared_genes))
#write.csv(shared_genes,"Shared_gene_OpenTarget.csv", row.names = FALSE)


#intersect with DisGeNET

# Read downloaded DisGeNET files
as_data_1 <- read_tsv("AS_curated_0.20.tsv")
pso_data_1 <- read_tsv("PsO_curated_0.20.tsv")

# Get gene symbols
as_genes_1 <- unique(as_data_1$Gene)
pso_genes_1 <- unique(pso_data_1$Gene)

# Find shared genes for DisGeNET
shared_genes_1 <- intersect(as_genes_1, pso_genes_1)
print(shared_genes_1)
print(length(shared_genes_1))


# Intersect with GeneCards

# Read downloaded GeneCards files
as_data_2 <- read_tsv("AS_GeneCard.tsv")
pso_data_2 <- read_tsv("PsO_GeneCard.tsv")

# Filter genes with Relevance score >= 0.3
as_filtered_1 <- as_data_2[as_data_2$`Relevance score` >= 0.3, ]
pso_filtered_1 <- pso_data_2[pso_data_2$`Relevance score` >= 0.3, ]

# Get gene symbols
as_genes_2 <- unique(as_filtered_1$`Gene Symbol`)
pso_genes_2 <- unique(pso_filtered_1$`Gene Symbol`)

# Find shared genes
shared_genes_2 <- intersect(as_genes_2, pso_genes_2)
print(shared_genes_2)
print(length(shared_genes_2))


# Overlapping

# Overlap between shared genes OpenTarget and DisGeNET
shared_genes_3 <- intersect(shared_genes, shared_genes_1)
print(shared_genes_3)
print(length(shared_genes_3))

# Overlap between shared genes OpenTarget and GeneCards
shared_genes_4 <- intersect(shared_genes, shared_genes_2)
print(shared_genes_4)
print(length(shared_genes_4))

# Overlap between shared genes DisGeNET and GeneCards
shared_genes_5 <- intersect(shared_genes_1, shared_genes_2)
print(shared_genes_5)
print(length(shared_genes_5))

# Overlap between shared genes OpenTarget, DisGeNET and GeneCards
shared_genes_6 <- intersect(shared_genes_3, shared_genes_4)
print(shared_genes_6)
print(length(shared_genes_6))

# Obtaining Scores

# Filter DisGeNET data for ScoreGDA
as_disgenet <- as_data_1 %>% filter(Gene %in% shared_genes_6) %>% 
  group_by(Gene) %>% summarise(DisGeNET_ScoreGDA_AS = mean(`ScoreGDA`, na.rm = TRUE))
print(as_disgenet)

pso_disgenet <- pso_data_1 %>% filter(Gene %in% shared_genes_6) %>% 
  group_by(Gene) %>% summarise(DisGeNET_ScoreGDA_PsO = mean(`ScoreGDA`, na.rm = TRUE))

# Filter GeneCards data for Relevance score
as_genecards <- as_data_2 %>% filter(`Gene Symbol` %in% shared_genes_6) %>% 
  group_by(`Gene Symbol`) %>% summarise(GeneCard_RelevanceScore_AS = mean(`Relevance score`, na.rm = TRUE))

pso_genecards <- pso_data_2 %>% filter(`Gene Symbol` %in% shared_genes_6) %>% 
  group_by(`Gene Symbol`) %>% summarise(GeneCard_RelevanceScore_PsO = mean(`Relevance score`, na.rm = TRUE))

# Filter OpenTarget data for GlobalScore
as_opentarget <- as_data %>% filter(symbol %in% shared_genes_6) %>%
  group_by(symbol) %>% summarise(OpenTarget_GlobalScore_AS = mean(globalScore, na.rm = TRUE))

pso_opentarget <- pso_data %>% filter(symbol %in% shared_genes_6) %>%
  group_by(symbol) %>% summarise(OpenTarget_GlobalScore_PsO = mean(globalScore, na.rm = TRUE))

# Merge all datasets
final_table <- as_disgenet %>%
  inner_join(pso_disgenet, by = "Gene") %>%
  inner_join(as_genecards, by = c("Gene" = "Gene Symbol")) %>%
  inner_join(pso_genecards, by = c("Gene" = "Gene Symbol")) %>%
  inner_join(as_opentarget, by = c("Gene" = "symbol")) %>%
  inner_join(pso_opentarget, by = c("Gene" = "symbol"))

# View result
print(final_table)
print(nrow(final_table))

# Export to CSV
write_csv(final_table, "shared_genes_6_scores_table.csv")

# # Create a named list of gene sets
# venn_list <- list(
#   "Open Targets" = shared_genes,
#   "DisGeNET" = shared_genes_1,
#   "GeneCards" = shared_genes_2
# )
# 
# # Output image
# venn.plot <- venn.diagram(
#   x = venn_list,
#   category.names = c("Open Targets", "DisGeNET", "GeneCards"),
#   filename = NULL,  # Use NULL for drawing in RStudio plot viewer
#   output = TRUE,
#   imagetype = "png",
#   height = 3000,
#   width = 3000,
#   resolution = 500,
#   lwd = 2,
#   col = "black",
#   fill = c("#E69F00", "#56B4E9", "#009E73"),
#   cex = 1.4,
#   cat.cex = 1.4,
#   cat.fontface = "bold"
# )
# 
# grid::grid.draw(venn.plot)

library(VennDiagram)
library(grid)

# Define gene sets
venn_list <- list(
  "Open Targets" = shared_genes,
  "DisGeNET" = shared_genes_1,
  "GeneCards" = shared_genes_2
)

# Create Venn diagram
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Open Targets", "DisGeNET", "GeneCards"),
  filename = NULL,
  imagetype = "png",
  height = 3000,
  width = 3000,
  resolution = 600,
  
  # Circle style
  lwd = 2,
  lty = "solid",
  col = "black",
  fill = c("#E69F00", "#56B4E9", "#009E73"),
  alpha = 0.6,
  
  # Number style
  cex = 2.0,
  fontface = "bold",
  fontfamily = "sans",
  
  # Category label styling and manual adjustment
  cat.cex = 2.0,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = "black",
  cat.pos = c(-20, 20, 0),          # angles for label direction
  cat.dist = c(0.08, 0.08, 0.10),   # distances from circle edge
  cat.just = list(c(0.5, 1), c(0.5, 0), c(0.5, 0.5)),  # justification tweaks
  
  # Margin
  margin = 0.1
)

# Draw with a title
grid.newpage()
grid.draw(venn.plot)
grid.text("Overlap of Shared Genes in AS and PsO", y = unit(0.95, "npc"), gp = gpar(fontsize = 18, fontface = "bold"))


library(rentrez)

# Set your email to comply with NCBI policy
options(rentrez.email = "qsiaun@gmail.com")


# Step 1: Load your gene score CSV
df <- read_csv("shared_genes_6_scores_table.csv")

# Step 2: Normalize GeneCards scores (optional)
df <- df %>%
  mutate(
    # GeneCards_AS_norm = rescale(GeneCard_RelevanceScore_AS, to = c(0, 1)),
    # GeneCards_PsO_norm = rescale(GeneCard_RelevanceScore_PsO, to = c(0, 1)),
    DisGeNET_avg = (DisGeNET_ScoreGDA_AS + DisGeNET_ScoreGDA_PsO) / 2,
    GeneCards_avg = (GeneCards_AS_norm + GeneCards_PsO_norm) / 2,
    OT_avg = (OpenTarget_GlobalScore_AS + OpenTarget_GlobalScore_PsO) / 2
  )

# Step 3: Function to get PubMed hit count
get_pubmed_hits <- function(gene, disease) {
  # Tukar nama penyakit kepada bentuk MeSH dan sinonim jika sesuai
  disease_query <- switch(tolower(disease),
                          "ankylosing spondylitis" = "(ankylosing spondylitis[Title/Abstract] OR spondyloarthritis[Title/Abstract] OR \"spondylitis, ankylosing\"[MeSH Terms])",
                          "psoriasis" = "(psoriasis[Title/Abstract] OR \"psoriasis\"[MeSH Terms])",
                          paste0("(", disease, "[Title/Abstract])")  # default fallback
  )
  
  # Gabungkan dengan gene query
  query <- paste0("(", gene, "[Title/Abstract]) AND ", disease_query)
  
  result <- tryCatch({
    Sys.sleep(0.3)  # Untuk elak NCBI block jika banyak carian
    entrez_search(db = "pubmed", term = query, retmax = 0)$count
  }, error = function(e) 0)
  
  as.numeric(result)
}

# Step 4: Apply PubMed search to each gene (may take a few minutes)
df <- df %>%
  rowwise() %>%
  mutate(
    PubMed_AS = get_pubmed_hits(Gene, "ankylosing spondylitis"),
    PubMed_PsO = get_pubmed_hits(Gene, "psoriasis")
  ) %>%
  ungroup()


# Step 5: Assign literature support score
df <- df %>%
  mutate(
    Lit_score = case_when(
      PubMed_AS >= 5 & PubMed_PsO >= 5 ~ 2,
      PubMed_AS >= 3 | PubMed_PsO >= 3 ~ 1,
      TRUE ~ 0
    )
  )

# Step 6: Calculate Composite Score (adjust weights here)
df <- df %>%
  mutate(
    Composite = 0.3 * DisGeNET_avg +
      0.2 * GeneCards_avg +
      0.2 * OT_avg +
      0.3 * (Lit_score / 2),  # normalize Lit_score 0–2 to 0–1
    Rank = rank(-Composite),
    
    # Tier = case_when(
    #   Composite >= 0.7 ~ "Tier 1",
    #   Composite >= 0.5 ~ "Tier 2",
    #   TRUE ~ "Tier 3"
    # )
)

# Kira ambang tertile berdasarkan Composite Score
tier_cutoffs <- quantile(df$Composite, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
print(tier_cutoffs)
# Tetapkan Tier berdasarkan quantile
df <- df %>%
  mutate(
    Tier = case_when(
      Composite >= tier_cutoffs[3] ~ "Tier 1",
      Composite >= tier_cutoffs[2] ~ "Tier 2",
      Composite >= tier_cutoffs[1] ~ "Tier 3",
      TRUE ~ "Tier 4"
    )
  )

# Step 7: Save output
write_csv(df, "prioritized_shared_genes_R.csv")

