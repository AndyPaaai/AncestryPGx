################################################################################
# 0. Introduction --------------------------------------------------------------
################################################################################
# This script reproduces the original workflow. Each numbered section is 
# self-contained, so collaborators can jump directly to the parts they need.

################################################################################
# 1. Load Required Libraries ---------------------------------------------------
################################################################################
# Uncomment the following lines to install any missing libraries
# pkgs <- c("readxl", "dplyr", "nortest", "ggplot2", "reshape2", "tidytext")
# install.packages(setdiff(pkgs, rownames(installed.packages())))

library(readxl)
library(circlize)
library(dplyr)
library(nortest)
library(ggplot2)
library(reshape2)
library(tidytext)
library(rentrez)
library(httr)
library(jsonlite)
library(xml2)
library(stringr)
library(purrr)
library(countrycode)
library(readr)
library(writexl)

################################################################################
# 2. Import Data ---------------------------------------------------------------
################################################################################
# Variant annotation and study parameter data (change paths as needed)
var_drug_ann_path      <- "~/data/var_drug_ann.xlsx"
study_parameters_path  <- "~/data/study_parameters.xlsx"
SNPs_FQs_path          <- "~/data/SNPs_FQs.csv"


stopifnot(
  file.exists(var_drug_ann_path),
  file.exists(study_parameters_path),
  file.exists(SNPs_FQs_path)
)

var_drug_ann     <- readxl::read_excel(var_drug_ann_path)
study_parameters <- readxl::read_excel(study_parameters_path)
SNPs_FQs         <- read.csv(SNPs_FQs_path, stringsAsFactors = FALSE)

################################################################################
# 3. Data Cleaning & Feature Engineering ---------------------------------------
################################################################################

# 3.1 Remove duplicated Variant Annotation IDs, keeping the row with the highest (Cases + Controls)
study_parameters_filtered <- study_parameters %>%
  group_by(`Variant Annotation ID`) %>%
  filter(n() > 1) %>%
  ungroup()

study_parameters_final_repeated <- study_parameters_filtered %>%
  group_by(`Variant Annotation ID`) %>%
  slice_max(order_by = (`Study Cases` + `Study Controls`), n = 1, with_ties = FALSE) %>%
  ungroup()

study_parameters_unique <- study_parameters %>%
  group_by(`Variant Annotation ID`) %>%
  filter(n() == 1) %>%
  ungroup()

study_parameters_final <- bind_rows(study_parameters_final_repeated, study_parameters_unique)

# 3.2 Merge with annotation table by Variant Annotation ID
df_final <- merge(study_parameters_final, var_drug_ann, by = "Variant Annotation ID", all = FALSE)

# 3.3 Filter for rs variants and single-letter alleles
data <- df_final %>%
  filter(grepl("^rs", `Variant/Haplotypes`)) %>%
  filter(grepl("^[A-Z]$", Alleles))

# 3.4 Extract significant variants
data_significant <- data %>%
  filter(Significance == "yes")

# 3.5 Save SNPs IDs (optional, change paths as needed)
# SNPs_IDs <- unique(data$`Variant/Haplotypes`)
# write.csv(SNPs_IDs, "~/ancestryPGx/data/SNPs_IDs.csv", row.names = FALSE)

# 3.6 Add population columns to data
colnames(SNPs_FQs)[colnames(SNPs_FQs) == "rs"] <- "Variant/Haplotypes"
populations <- c("PLQ", "ATQCES", "CLM", "CHG", "ATQPGC")
for (population in populations) {
  data[[population]] <- NA
}
for (i in 1:nrow(data)) {
  variant <- data$`Variant/Haplotypes`[i]
  for (population in populations) {
    subset_SNPs <- SNPs_FQs[SNPs_FQs$`Variant/Haplotypes` == variant & SNPs_FQs$Population.Group == population, ]
    if (nrow(subset_SNPs) > 0) {
      allele <- data$Alleles[i]
      frequency <- subset_SNPs[[allele]][1]
      data[i, population] <- frequency
    }
  }
}

# 3.7 Add population columns to significant data
for (population in populations) {
  data_significant[[population]] <- NA
}
for (i in 1:nrow(data_significant)) {
  variant <- data_significant$`Variant/Haplotypes`[i]
  for (population in populations) {
    subset_SNPs <- SNPs_FQs[SNPs_FQs$`Variant/Haplotypes` == variant & SNPs_FQs$Population.Group == population, ]
    if (nrow(subset_SNPs) > 0) {
      allele <- data_significant$Alleles[i]
      frequency <- subset_SNPs[[allele]][1]
      data_significant[i, population] <- frequency
    }
  }
}

# 3.8 Calculate mean frequency by ancestry for significant data
african_populations <- c("PLQ", "CHG")
european_populations <- c("ATQCES", "CLM", "ATQPGC")
data_significant$AV_AFR <- rowMeans(data_significant[african_populations], na.rm = TRUE)
data_significant$AV_EUR <- rowMeans(data_significant[european_populations], na.rm = TRUE)

# 3.9 Remove rows with all NA values in population columns
cols_to_check <- populations
data_significant <- data_significant[rowSums(is.na(data_significant[, cols_to_check])) != length(cols_to_check), ]

################################################################################
# 4. Correlation Analysis: Significant Data ------------------------------------
################################################################################

# 4.1 Count non-NA values per population
population_counts <- sapply(populations, function(p) sum(!is.na(data_significant[[p]])))
print(population_counts)

# 4.2 Normality Test (Shapiro-Wilk)
shapiro_results <- lapply(populations, function(p) shapiro.test(na.omit(data_significant[[p]])))
names(shapiro_results) <- populations
print(shapiro_results)

# 4.3 Spearman Correlation Matrix
data_populations_sig <- data_significant[, populations]
correlation_matrix <- cor(data_populations_sig, method = "spearman", use = "pairwise.complete.obs")
print(correlation_matrix)

################################################################################
# 5. Visualization: Correlation Plot -------------------------------------------
################################################################################
# Create heatmap for correlation matrix
correlation_melted <- melt(correlation_matrix)
order <- c("ATQCES", "ATQPGC", "CLM", "PLQ", "CHG")
correlation_melted$Var1 <- factor(correlation_melted$Var1, levels = order)
correlation_melted$Var2 <- factor(correlation_melted$Var2, levels = order)
correlation_melted <- correlation_melted[as.numeric(correlation_melted$Var1) <= as.numeric(correlation_melted$Var2), ]

ggplot(data = correlation_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab",
                       name = "Correlation") +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 50, vjust = 1, hjust = 1),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank()) +
  coord_fixed() +
  geom_text(aes(label = round(value, 2)), color = "black", size = 3.7) +
  guides(fill = guide_colorbar(title = "Correlation")) +
  scale_x_discrete(limits = order) +
  scale_y_discrete(limits = rev(order))

################################################################################
# 6. Bar Plots: Drugs ----------------------------------------------------------
################################################################################
# Filter and aggregate data for plotting drug-related SNPs
filtered_data <- data_significant %>%
  filter(!is.na(`Drug(s)`), !is.na(`Phenotype Category`), !is.na(`Direction of effect`))

top_10_drugs <- filtered_data %>%
  count(`Drug(s)`, sort = TRUE) %>%
  top_n(10, n) %>%
  pull(`Drug(s)`)

filtered_data <- filtered_data %>%
  filter(`Drug(s)` %in% top_10_drugs) %>%
  mutate(`Phenotype_Direction` = paste(`Phenotype Category`, `Direction of effect`, sep = " - "))

agg_data <- filtered_data %>%
  group_by(`Drug(s)`, `Phenotype_Direction`) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  group_by(`Drug(s)`) %>%
  mutate(total_count = sum(count)) %>%
  ungroup() %>%
  arrange(desc(total_count), `Drug(s)`)

# Define colors for phenotype categories
keep_values <- c("Dosage - decreased", "Dosage - increased",
                 "Efficacy - decreased", "Efficacy - increased",
                 "Metabolism/PK - increased", "Metabolism/PK - decreased")

agg_data <- agg_data %>%
  mutate(Phenotype_Direction = ifelse(Phenotype_Direction %in% keep_values,
                                      Phenotype_Direction, "Other"))

color_values <- c(
  "Dosage - decreased" = "#E63946",
  "Dosage - increased" = "#F4A261",
  "Efficacy - decreased" = "#2A9D8F",
  "Efficacy - increased" = "#264653",
  "Metabolism/PK - decreased" = "#6A0572",
  "Metabolism/PK - increased" = "#457B9D",
  "Other" = "#B0B0B0"
)

# Generate bar plot for top 10 drugs
ggplot(agg_data, aes(x = count, y = reorder(`Drug(s)`, total_count), fill = `Phenotype_Direction`)) +
  geom_vline(xintercept = c(20, 40, 60, 120, 140),
             linetype = "dashed", color = "grey50", size = 0.5) +
  geom_bar(stat = "identity", position = 'stack', width = 0.8) +
  scale_fill_manual(values = color_values) +
  scale_x_continuous(breaks = c(20, 40, 60, 120, 140)) +
  labs(x = "SNPs", y = "Drugs", fill = "Phenotypes") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 7. Bar Plots: Genes ----------------------------------------------------------
################################################################################
# Filter and aggregate data for plotting gene-related SNPs
colnames(data_significant) <- gsub("\\s+", "_", colnames(data_significant))

filtered_data <- data_significant %>%
  filter(!is.na(Gene)) %>%
  mutate(Phenotype_Group = case_when(
    Phenotype_Category %in% c("Dosage", "Efficacy", "Metabolism/PK") &
      Direction_of_effect %in% c("decreased", "increased") ~ paste(Phenotype_Category, "-", Direction_of_effect),
    TRUE ~ "Other"
  ))

top_genes <- filtered_data %>%
  count(Gene, name = "Total") %>%
  arrange(desc(Total)) %>%
  slice_head(n = 15)

plot_data <- filtered_data %>%
  filter(Gene %in% top_genes$Gene) %>%
  count(Gene, Phenotype_Group, name = "Frequency")

gene_order <- plot_data %>%
  group_by(Gene) %>%
  summarise(Total_Frequency = sum(Frequency)) %>%
  arrange(Total_Frequency) %>%
  pull(Gene)

# Define colors for phenotype groups
color_values <- c(
  "Dosage - decreased" = "#E63946",
  "Dosage - increased" = "#F4A261",
  "Efficacy - decreased" = "#2A9D8F",
  "Efficacy - increased" = "#264653",
  "Metabolism/PK - decreased" = "#6A0572",
  "Metabolism/PK - increased" = "#457B9D",
  "Other" = "#B0B0B0"
)

# Generate bar plot for top 15 genes
ggplot(plot_data, aes(x = Frequency, y = factor(Gene, levels = gene_order), fill = Phenotype_Group)) +
  geom_vline(xintercept = seq(0, 70, by = 10), linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = color_values) +
  scale_x_continuous(breaks = seq(0, 70, by = 10), limits = c(0, 80)) +
  labs(x = "SNPs", y = "Genes", fill = "Phenotypes") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank())

################################################################################
# 8. Jitter Plot: Ancestry Frequency -------------------------------------------
################################################################################
# Categorize significant SNPs by ancestry predominance
data_significant$Group <- "No predominance"
for (i in 1:nrow(data_significant)) { 
  if (!is.na(data_significant$PLQ[i]) && !is.na(data_significant$AV_EUR[i])) {
    if (data_significant$PLQ[i] > 0.75 && data_significant$AV_EUR[i] < 0.5) {
      data_significant$Group[i] <- "PLQ predominance"
    } else if (data_significant$AV_EUR[i] > 0.75 && data_significant$PLQ[i] < 0.5) {
      data_significant$Group[i] <- "European predominance"
    }
  }
}

# Create the concatenation of rs + allele
data_significant <- data_significant %>%
  mutate(Variant_Allele = paste(`Variant/Haplotypes`, Alleles, sep = "-"))

# Filter significant SNPs by group
PLQ_predominance <- data_significant %>%
  filter(Group == "PLQ predominance")

European_predominance <- data_significant %>%
  filter(Group == "European predominance")

# Optional: Export group data to Excel
# write_xlsx(PLQ_predominance, "PLQ_predominance.xlsx")
# write_xlsx(European_predominance, "European_predominance.xlsx")

data_significant_unique <- data_significant %>%
  distinct(Variant_Allele, .keep_all = TRUE)

# Generate jitter plot for ancestry frequencies
ggplot(data_significant_unique, aes(x = AV_EUR, y = PLQ, color = Group)) +
  geom_jitter(width = 0.03, height = 0.05) +
  scale_color_manual(values = c(
    "PLQ predominance" = "red",
    "European predominance" = "blue",
    "No predominance" = "gray"
  )) +
  geom_hline(yintercept = c(0.5, 0.75), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(0.5, 0.75), linetype = "dashed", color = "black") +
  labs(x = "European ancestry",
       y = "African ancestry (PLQ)") +
  theme_classic() +
  theme(legend.position = "none")

################################################################################
# 9. Top SNPs by Population ----------------------------------------------------
################################################################################
# Extract and order top SNPs for African and European ancestry

# Top 5 SNPs with highest AV_AFR in African ancestry
top_afr <- data_significant[data_significant$Group == "PLQ predominance", ]
top_afr <- top_afr[order(-top_afr$AV_AFR), ]

# Top 5 SNPs with highest AV_EUR in European ancestry
top_eur <- data_significant[data_significant$Group == "European predominance", ]
top_eur <- top_eur[order(-top_eur$AV_EUR), ]

# Count unique SNPs for African and European populations
unique_snps_afr <- length(unique(top_afr$`Variant/Haplotypes`))
unique_snps_eur <- length(unique(top_eur$`Variant/Haplotypes`))

# Calculate absolute difference between AV_AFR and AV_EUR
top_afr$dif_AV <- abs(top_afr$AV_AFR - top_afr$AV_EUR)
top_eur$dif_AV <- abs(top_eur$AV_AFR - top_eur$AV_EUR)

# Order by largest difference
top_afr_ordered <- top_afr[order(-top_afr$dif_AV), ]
top_eur_ordered <- top_eur[order(-top_eur$dif_AV), ]

################################################################################
# 10. Ancestry Analysis --------------------------------------------------------
################################################################################
# Analyze ancestry distribution and retain rows with highest sample size

# Retain rows with the highest total sample size per PMID
data_unique <- data %>%
  mutate(total_n = `Study Cases` + `Study Controls`) %>%
  group_by(PMID) %>%
  slice_max(total_n, n = 1, with_ties = FALSE) %>%
  ungroup()

################################################################################
# 11. Ancestry Analysis --------------------------------------------------------
################################################################################
# This section analyzes ancestry categories in the dataset, assigning each row 
# to a major population group and summarizing sample sizes for each ancestry.

################################################################################
# 11.1 Load Required Libraries -------------------------------------------------
################################################################################
# Install and load all required libraries for ancestry analysis
packs <- c("dplyr", "tidyr", "stringr", "purrr", "tibble", "circlize")
install.packages(setdiff(packs, rownames(installed.packages())), quiet = TRUE)
invisible(lapply(packs, library, character.only = TRUE))

################################################################################
# 11.2 Retain Largest Sample Per PMID ------------------------------------------
################################################################################
# For each PMID, keep the row with the largest total sample size
data_unique <- data %>%
  mutate(total_n = `Study Cases` + `Study Controls`) %>%
  group_by(PMID) %>%
  slice_max(total_n, n = 1, with_ties = FALSE) %>%
  ungroup()

################################################################################
# 11.3 Define Synonym Dictionary for Ancestry Categories -----------------------
################################################################################
# List of key terms and synonyms for each ancestry category
dict <- list(
  European = c(
    "caucasian", "white(?!\\s*and)", "white\\s*european", "european",
    "europeans", "british", "french", "german", "italian", "dutch",
    "spanish", "portuguese", "swiss", "nordic", "finnish", "swedish",
    "norwegian", "greek", "turkish", "jewish", "israeli", "arab", "iberian"
  ),
  Asian = c(
    "asian", "east\\s*asian", "central/south\\s*asian",
    "chinese", "han", "korean", "japanese", "malay", "singapore",
    "indian", "hindu", "pakistani", "punjabi", "pathan", "kashmiri",
    "thai", "taiwanese", "filipino", "hmong", "ughur",
    "maori", "polynesian", "pacific(?!\\s*islander)", "pacific\\s*islander"
  ),
  American = c(
    "\\bamerican\\b", "latino", "hispanic", "mexican", "brazilian",
    "chilean", "colombian", "peruvian", "argentinian", "puerto\\s*rican",
    "native\\s*american", "amerindian", "indigenous", "torres\\s*strait",
    "brown\\b", "pardo"
  ),
  African = c(
    "african", "afro", "black", "sub\\s*-?saharan",
    "ethiopian", "eritrean", "somali", "kenyan", "tanzanian", "ugandan",
    "rwandan", "burund(i|ian)", "ivory\\s*coast", "liberian", "ghanaian",
    "nigerian", "senegalese", "afro[-\\s]?caribbean"
  ),
  Other = c(
    "other(?!\\s*wise)", "mixed", "multiracial", "biracial",
    "admixed", "various", "miscellaneous"
  ),
  Unknown = c("unknown", "^na$", "not\\s*reported", "unreported", "unspecified")
)
patterns <- purrr::map(dict, ~ paste0("(?i)\\b(", paste(.x, collapse = "|"), ")\\b"))

################################################################################
# 11.4 Helper Function: Assign Category ----------------------------------------
################################################################################
# Assign a text segment to an ancestry category based on the dictionary
get_category <- function(text_segment) {
  for (cat in names(patterns)) {
    if (str_detect(text_segment, patterns[[cat]])) return(cat)
  }
  return("Other")
}

################################################################################
# 11.5 Helper Function: Parse Row ----------------------------------------------
################################################################################
# Parse a row to estimate ancestry breakdown based on description and sample size
parse_row <- function(desc, cases, controls) {
  if (is.na(desc) || str_detect(desc, patterns$Unknown)) {
    return(
      tibble(Ancestry = "Unknown",
             Cases     = cases,
             Controls  = controls)
    )
  }
  total_sample  <- cases + controls
  w_cases       <- ifelse(total_sample == 0, 0, cases    / total_sample)
  w_controls    <- ifelse(total_sample == 0, 0, controls / total_sample)
  desc_l        <- str_replace_all(tolower(desc), "[\n\r]", " ")
  token_regex   <- "(\\d+\\.?\\d*)\\s*%?[^;,.\\n]*?"
  tokens        <- str_extract_all(desc_l, token_regex)[[1]]
  out <- tibble()
  used_any <- FALSE
  for (tok in tokens) {
    cat_here <- get_category(tok)
    if (cat_here %in% c("Other", "Unknown")) next
    number_val <- as.numeric(str_extract(tok, "\\d+\\.?\\d*"))
    is_pct     <- str_detect(tok, "%")
    if (is_pct) {
      case_n    <- cases    * number_val / 100
      ctrl_n    <- controls * number_val / 100
    } else {
      case_n    <- number_val * w_cases
      ctrl_n    <- number_val * w_controls
    }
    out <- bind_rows(out,
                     tibble(Ancestry = cat_here,
                            Cases     = case_n,
                            Controls  = ctrl_n))
    used_any <- TRUE
  }
  if (!used_any) {
    cat_row <- get_category(desc_l)
    out     <- tibble(Ancestry = cat_row,
                      Cases     = cases,
                      Controls  = controls)
  }
  out
}

################################################################################
# 11.6 Apply Ancestry Parsing to All Rows --------------------------------------
################################################################################
# Apply the parsing function to each row in the dataset
parsed <- pmap_dfr(
  list(
    desc     = data_unique$`Biogeographical Groups`,
    cases    = data_unique$`Study Cases`,
    controls = data_unique$`Study Controls`
  ),
  parse_row
)

################################################################################
# 11.7 Summarize and Complete Category Counts ----------------------------------
################################################################################
# Aggregate case/control totals for each ancestry, ensuring all categories exist
totals <- parsed %>%
  group_by(Ancestry) %>%
  summarise(
    Total_Cases    = sum(Cases,    na.rm = TRUE),
    Total_Controls = sum(Controls, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Total_Sample = Total_Cases + Total_Controls) %>%
  complete(
    Ancestry = c("European", "Asian", "American", "African", "Other", "Unknown"),
    fill     = list(Total_Cases = 0, Total_Controls = 0, Total_Sample = 0)
  ) %>%
  arrange(match(Ancestry, c("European","Asian","American","African","Other","Unknown")))

################################################################################
# 11.8 Print Final Ancestry Summary --------------------------------------------
################################################################################
print(totals)
# Optionally, remove the 'Other' category if desired:
# totals <- totals %>% filter(Ancestry != "Other")

################################################################################
# 11.9 Chord Diagram: Ancestry x Case/Control ----------------------------------
################################################################################
# Visualize ancestry categories and sample sizes using a chord diagram

# Prepare data for chord diagram
links <- totals %>%
  rename(from = Ancestry) %>%
  pivot_longer(cols = c("Total_Cases", "Total_Controls"),
               names_to = "to", values_to = "value") %>%
  mutate(to = recode(to,
                     "Total_Cases"    = "Cases",
                     "Total_Controls" = "Controls")) %>%
  filter(value > 0) %>%
  select(from, to, value)

# Determine sector order for chord diagram
ancestry_order <- totals %>%
  arrange(desc(Total_Sample)) %>%
  pull(Ancestry)
sector.order <- c("Cases", "Controls", ancestry_order)

# Define colors and transparency for sectors
grid.col <- c(
  "Cases"     = "#999999",
  "Controls"  = "#222222",
  "European"  = "#5fa8d3",
  "Asian"     = "#1b4965",
  "American"  = "#62b6cb",
  "African"   = "#e63946",
  "Other"     = "#bee9e8",
  "Unknown"   = "#cae9ff"
)[sector.order]
transparency.map <- c(
  European = 0,
  Asian    = 1,
  American = 1,
  African  = 1,
  Other    = 1,
  Unknown  = 1
)
trans_vec <- ifelse(links$from %in% names(transparency.map),
                    transparency.map[as.character(links$from)],
                    0.1)

# Convert columns to factors with correct levels
links$from <- factor(links$from, levels = sector.order)
links$to   <- factor(links$to,   levels = sector.order)

# Draw the chord diagram
circos.clear()
chordDiagram(
  links,
  order        = sector.order,
  grid.col     = grid.col,
  transparency = trans_vec,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.01)
)

# Count NAs in the 'Biogeographical Groups' column for reference
na_count <- sum(is.na(data_unique$`Biogeographical Groups`))

################################################################################
# 12. Bibliometric Analysis ----------------------------------------------------
################################################################################
# This section gathers bibliometric data for the publications in the dataset.
# It retrieves details such as DOI, title, publication date, journal, country, 
# number of authors, open access status, and citation count using APIs.

################################################################################
# 12.1 Load Required Libraries -------------------------------------------------
################################################################################
library(httr)         # HTTP requests
library(xml2)         # XML parsing
library(jsonlite)     # JSON parsing
library(stringr)      # String utilities
library(dplyr)        # Data manipulation
library(countrycode)  # Country normalization
library(purrr)        # Functional programming

################################################################################
# 12.2 Set up API Configuration ------------------------------------------------
################################################################################
my_email <- "your_email@example.com"                # <--- Place your email here
my_tool  <- "biblio_script"
my_key   <- "your_ncbi_api_key"                     # <--- Place your NCBI API key here

# Null-coalescing operator for convenience
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

# Empty template
empty_info <- list(
  doi       = NA_character_,
  title     = NA_character_,
  pub_date  = NA_character_,
  journal   = NA_character_,
  country   = NA_character_,
  n_authors = NA_integer_
)

################################################################################
# 12.3 Safe GET Helper for API Calls -------------------------------------------
################################################################################
safe_get <- function(url, q, tries = 3) {
  RETRY(
    "GET", url,
    query  = q,
    times  = tries,
    pause_base = 0.4,
    pause_cap  = 2,
    timeout(20)
  )
}

################################################################################
# 12.4 PubMed Article Parsing --------------------------------------------------
################################################################################
parse_article <- function(node) {
  doi <- xml_text(xml_find_first(node, ".//ArticleId[@IdType='doi']"))
  if (is.na(doi) || doi == "") {
    doi_alt <- xml_text(xml_find_first(node, ".//ELocationID[@EIdType='doi']"))
    if (!is.na(doi_alt) && doi_alt != "") doi <- doi_alt
  }
  doi <- ifelse(doi == "", NA_character_, doi)
  
  title   <- xml_text(xml_find_first(node, ".//ArticleTitle"))
  journal <- xml_text(xml_find_first(node, ".//Journal//Title"))
  
  date_node <- xml_find_first(node, ".//PubDate")
  pub_date  <- if (!inherits(date_node, "xml_missing")) {
    vals <- c(
      xml_text(xml_find_first(date_node, ".//Year")),
      xml_text(xml_find_first(date_node, ".//Month")),
      xml_text(xml_find_first(date_node, ".//Day")))
    vals <- vals[vals != ""]
    if (length(vals) == 0) NA_character_ else paste(vals, collapse = "-")
  } else NA_character_
  
  n_authors <- length(xml_find_all(node, ".//Author"))
  
  aff <- xml_text(xml_find_first(
    node, ".//AuthorList/Author[1]/AffiliationInfo/Affiliation"))
  country <- NA_character_
  if (!is.na(aff) && aff != "") {
    last_chunk <- str_trim(tail(str_split(aff, "[,;]")[[1]], 1))
    cc <- suppressWarnings(countrycode(last_chunk,
                                       origin = "country.name",
                                       destination = "country.name"))
    country <- ifelse(is.na(cc), last_chunk, cc)
  }
  
  list(
    doi       = doi,
    title     = title   %||% NA_character_,
    pub_date  = pub_date,
    journal   = journal %||% NA_character_,
    country   = country,
    n_authors = n_authors
  )
}

################################################################################
# 12.5 Retrieve Metadata for Blocks of PMIDs -----------------------------------
################################################################################
get_pubmed_block <- function(pmids) {
  ids  <- paste(pmids, collapse = ",")
  url  <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  res  <- safe_get(
    url,
    q = list(
      db      = "pubmed",
      id      = ids,
      retmode = "xml",
      rettype = "abstract",
      tool    = my_tool,
      email   = my_email,
      api_key = my_key
    ))
  if (status_code(res) != 200)
    stop("HTTP error ", status_code(res), " while getting PMIDs")
  
  xml_doc  <- read_xml(content(res, "raw"))
  articles <- xml_find_all(xml_doc, ".//PubmedArticle")
  
  purrr::map_dfr(articles, parse_article,
                 .id = "PMID") %>%
    mutate(PMID = pmids[as.integer(PMID)])
}

################################################################################
# 12.6 Open Access and Citation Count Helpers ----------------------------------
################################################################################
get_open_access_status <- function(doi) {
  if (is.na(doi) || doi == "") return(NA_character_)
  url <- paste0("https://api.unpaywall.org/v2/",
                URLencode(doi, reserved = TRUE))
  res <- safe_get(url, q = list(email = my_email))
  if (status_code(res) != 200) return(NA_character_)
  data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  if (is.null(data$is_oa)) NA_character_ else ifelse(data$is_oa, "Open", "Closed")
}

get_citation_count <- function(doi) {
  if (is.na(doi) || doi == "") return(NA_integer_)
  url <- paste0("https://api.semanticscholar.org/graph/v1/paper/",
                URLencode(doi, reserved = TRUE))
  res <- safe_get(url, q = list(fields = "citationCount"))
  if (status_code(res) != 200) return(NA_integer_)
  data <- fromJSON(content(res, "text", encoding = "UTF-8"))
  data$citationCount %||% NA_integer_
}

################################################################################
# 12.7 Build the Bibliometric Data Frame ---------------------------------------
################################################################################
pmid_vec <- data %>%
  mutate(PMID = as.character(PMID)) %>%
  distinct(PMID) %>%
  pull(PMID)

chunks <- split(pmid_vec, ceiling(seq_along(pmid_vec) / 200))

# Collect core metadata
biblio_core <- map_dfr(chunks, get_pubmed_block)

# Add open access and citations (can be parallelized)
biblio_final <- biblio_core %>%
  rowwise() %>%
  mutate(
    Access    = get_open_access_status(doi),
    Citations = get_citation_count(doi)
  ) %>%
  ungroup() %>%
  rename_with(~ stringr::str_to_title(.x), everything()) %>%
  select(Pmid, Doi, Title, Pub_date, Journal,
         Country, N_authors, Access, Citations)

# Convert publication date to year (numeric)
biblio_final$Pub_date <- as.numeric(sub("-.*", "", biblio_final$Pub_date))

# Optional: Save as Excel
# write_xlsx(biblio_final, "biblio_final.xlsx")

################################################################################
# 13. Journal Metrics: Quartile and H-index Matching ---------------------------
################################################################################
# This section matches the bibliometric data with SJR (Scimago Journal & Country Rank)
# metrics to annotate journals by quartile and H-index.

################################################################################
# 13.1 Load Required Libraries for Journal Metrics -----------------------------
################################################################################
library(tidyverse)
library(janitor)
library(fuzzyjoin)
library(future.apply)

################################################################################
# 13.2 Load SJR Data and Prepare for Join --------------------------------------
################################################################################
sjr_dir   <- "/path/to/SJR"
sjr_files <- list.files(sjr_dir, "^scimagojr \\d{4}\\.csv$", full.names = TRUE)

sjr_all <- map_dfr(sjr_files, function(f) {
  yr <- as.integer(str_extract(basename(f), "\\d{4}"))
  read_delim(f, delim = ";", show_col_types = FALSE, locale = locale(encoding = "UTF-8")) %>% 
    clean_names() %>% 
    select(title, h_index, sjr_best_quartile) %>% 
    mutate(
      year        = yr,
      title_clean = title %>% 
        str_to_lower() %>% 
        str_replace_all("[^a-z0-9]", "") %>% 
        str_squish()
    )
})

################################################################################
# 13.3 Prepare Bibliometric Data for Join --------------------------------------
################################################################################
biblio_prepped <- biblio_final %>% 
  mutate(
    year          = as.integer(Pub_date),
    journal_clean = Journal %>% 
      str_to_lower() %>% 
      str_replace_all("[^a-z0-9]", "") %>% 
      str_squish()
  )

################################################################################
# 13.4 Exact Join by Year and Journal Name -------------------------------------
################################################################################
joined_exact <- biblio_prepped %>% 
  left_join(
    sjr_all %>% select(year, title_clean, h_index, sjr_best_quartile),
    by = c("year", "journal_clean" = "title_clean")
  )

################################################################################
# 13.5 Fuzzy Join for Unmatched Journals ---------------------------------------
################################################################################
need_fuzzy <- joined_exact %>% filter(is.na(h_index))

if (nrow(need_fuzzy) > 0) {
  key_unmatched <- need_fuzzy %>% 
    distinct(year, journal_clean)
  future::plan(multisession, workers = parallel::detectCores() - 1)
  fuzzy_matches <- future_lapply(
    split(key_unmatched, key_unmatched$year),
    function(df_year) {
      yr        <- df_year$year[1]
      sjr_year  <- sjr_all %>% 
        filter(year == yr) %>%    
        select(title_clean, h_index, sjr_best_quartile)
      stringdist_left_join(
        df_year, sjr_year,
        by           = c("journal_clean" = "title_clean"),
        method       = "jw",
        max_dist     = 0.10,
        distance_col = "dist"
      ) %>% 
        group_by(year, journal_clean) %>% 
        slice_min(dist, n = 1, with_ties = FALSE) %>% 
        ungroup() %>% 
        select(-dist)
    }
  ) %>% bind_rows()
  future::plan(sequential)
  joined_all <- joined_exact %>% 
    left_join(
      fuzzy_matches,
      by = c("year", "journal_clean")
    ) %>% 
    mutate(
      h_index            = coalesce(h_index.x,  h_index.y),
      sjr_best_quartile  = coalesce(sjr_best_quartile.x,
                                    sjr_best_quartile.y)
    ) %>% 
    select(-ends_with(".x"), -ends_with(".y"))
} else {
  joined_all <- joined_exact
}

################################################################################
# 13.6 Finalize Bibliometric Data with Metrics ---------------------------------
################################################################################
biblio_final <- joined_all %>% 
  mutate(
    H_index_SJR       = h_index,
    Best_Quartile_SJR = sjr_best_quartile
  ) %>% 
  select(-journal_clean)

# Print unmatched count
cat("Unmatched records after fuzzy join:", biblio_final %>% filter(is.na(H_index_SJR)) %>% nrow(), "\n")

# Delete repeated columns
biblio_final$H_index_SJR <- NULL
biblio_final$Best_Quartile_SJR <- NULL
biblio_final$year <- NULL
biblio_final$title_clean <- NULL


# Optional: Save as Excel
# write_xlsx(biblio_final, "biblio_final_results.xlsx")

################################################################################
# 14. Country Cleaning and Region/Income Annotation ----------------------------
################################################################################
# This section standardizes country names and annotates each publication with
# WHO region and World Bank income group using manual classifications.

################################################################################
# 14.1 Load Required Libraries -------------------------------------------------
################################################################################
library(countrycode)
library(dplyr)

################################################################################
# 14.2 Country Cleaning Function -----------------------------------------------
################################################################################
clean_country <- function(country){
  country_clean <- countrycode(country, "country.name", "country.name")
  country_clean <- case_when(
    grepl("UK|Scotland|England|King's College", country, ignore.case = TRUE) ~ "United Kingdom",
    grepl("Brasil|Brazil|Instituto Estadual do Cérebro|Paulo Niemeyer|rumotecnologia", country, ignore.case = TRUE) ~ "Brazil",
    grepl("México|Mexico", country, ignore.case = TRUE) ~ "Mexico",
    grepl("Shanghai|Beijing|Xinjiang|Sichuan", country, ignore.case = TRUE) ~ "China",
    grepl("Hong Kong", country, ignore.case = TRUE) ~ "Hong Kong",
    grepl("Taichung|Taiwan|Yilan", country, ignore.case = TRUE) ~ "Taiwan",
    grepl("Sorbonne|Paris", country, ignore.case = TRUE) ~ "France",
    grepl("Pisa|benedetti.francesco", country, ignore.case = TRUE) ~ "Italy",
    grepl("Seattle|Calif|CA\\.|California|Davis", country, ignore.case = TRUE) ~ "United States",
    grepl("OH\\.|Illinois|IL\\.|Maryland|NY\\.|Minnesota|Florida|Ohio|North Carolina|Delaware|Missouri|Arkansas|Philadelphia|TN\\.|Texas|Colorado|Michigan|PA\\.|MA\\.|Iowa|Bethesda|St\\. Louis|Univ\\.|Johns Hopkins", country, ignore.case = TRUE) ~ "United States",
    grepl("Hanoi", country, ignore.case = TRUE) ~ "Vietnam",
    grepl("Spgrain", country, ignore.case = TRUE) ~ "Spain",
    grepl("Lund|Stockholm", country, ignore.case = TRUE) ~ "Sweden",
    grepl("Amedeo di Savoia", country, ignore.case = TRUE) ~ "Italy",
    grepl("Patras", country, ignore.case = TRUE) ~ "Greece",
    grepl("Montreal", country, ignore.case = TRUE) ~ "Canada",
    grepl("Haiti", country, ignore.case = TRUE) ~ "Haiti",
    grepl("Pharmacology|Pharmacy|Pharmaceutical|Biochemistry|Genetics|Department|Hospital|Medicine|Veterans Affairs|Health|Data Solutions|Genotyping|Affiliations|Electronic|address|gmail|hotmail|yahoo|edu|com|uk|us|net", country, ignore.case = TRUE) ~ NA_character_,
    is.na(country_clean) ~ NA_character_,
    TRUE ~ country_clean
  )
  return(country_clean)
}

################################################################################
# 14.3 Apply Country Cleaning --------------------------------------------------
################################################################################
biblio_final <- biblio_final %>%
  mutate(Country = clean_country(Country))

################################################################################
# 14.4 Annotate WHO Region and Income Group ------------------------------------
################################################################################
# Manual WHO Region classification
region_manual <- c(
  "Japan" = "Western Pacific", "Canada" = "Americas", "Israel" = "European", 
  "France" = "European", "United States" = "Americas", "Germany" = "European", 
  "China" = "Western Pacific", "United Kingdom" = "European", "Taiwan" = "Western Pacific",
  "Italy" = "European", "Finland" = "European", "Portugal" = "European", 
  "Slovenia" = "European", "Sweden" = "European", "South Korea" = "Western Pacific", 
  "India" = "South-East Asia", "Netherlands" = "European", "Spain" = "European", 
  "Austria" = "European", "Poland" = "European", "Australia" = "Western Pacific", 
  "Benin" = "African", "Greece" = "European", "Brazil" = "Americas", 
  "Singapore" = "Western Pacific", "Malaysia" = "Western Pacific", "Belgium" = "European", 
  "Hungary" = "European", "Lebanon" = "Eastern Mediterranean", "Russia" = "European", 
  "Egypt" = "Eastern Mediterranean", "Norway" = "European", "Mexico" = "Americas", 
  "Turkey" = "European", "Croatia" = "European", "Argentina" = "Americas", 
  "New Zealand" = "Western Pacific", "Hong Kong" = "Western Pacific", "Switzerland" = "European", 
  "Haiti" = "Americas", "Chile" = "Americas", "Estonia" = "European", "Romania" = "European", 
  "Thailand" = "South-East Asia", "South Africa" = "African", "Nigeria" = "African", 
  "Denmark" = "European", "Ireland" = "European", "Slovakia" = "European", 
  "Iran" = "Eastern Mediterranean", "Colombia" = "Americas", "Morocco" = "Eastern Mediterranean", 
  "Bosnia & Herzegovina" = "European", "Tanzania" = "African", "Tunisia" = "Eastern Mediterranean", 
  "Saudi Arabia" = "Eastern Mediterranean", "Serbia" = "European", 
  "Palestinian Territories" = "Eastern Mediterranean", "Czechia" = "European", 
  "Pakistan" = "Eastern Mediterranean", "United Arab Emirates" = "Eastern Mediterranean", 
  "Kenya" = "African", "Kazakhstan" = "European", "Peru" = "Americas", "Vietnam" = "Western Pacific", 
  "Latvia" = "European", "Bangladesh" = "South-East Asia", "Indonesia" = "South-East Asia", 
  "Algeria" = "African"
)
# Manual World Bank Income Group classification
income_manual <- c(
  "Japan" = "High income", "Canada" = "High income", "Israel" = "High income",
  "France" = "High income", "United States" = "High income", "Germany" = "High income",
  "China" = "Upper middle income", "United Kingdom" = "High income", "Taiwan" = "High income",
  "Italy" = "High income", "Finland" = "High income", "Portugal" = "High income",
  "Slovenia" = "High income", "Sweden" = "High income", "South Korea" = "High income",
  "India" = "Lower middle income", "Netherlands" = "High income", "Spain" = "High income",
  "Austria" = "High income", "Poland" = "High income", "Australia" = "High income",
  "Benin" = "Low income", "Greece" = "High income", "Brazil" = "Upper middle income",
  "Singapore" = "High income", "Malaysia" = "Upper middle income", "Belgium" = "High income",
  "Hungary" = "High income", "Lebanon" = "Upper middle income", "Russia" = "Upper middle income",
  "Egypt" = "Lower middle income", "Norway" = "High income", "Mexico" = "Upper middle income",
  "Turkey" = "Upper middle income", "Croatia" = "High income", "Argentina" = "Upper middle income",
  "New Zealand" = "High income", "Hong Kong" = "High income", "Switzerland" = "High income",
  "Haiti" = "Low income", "Chile" = "High income", "Estonia" = "High income",
  "Romania" = "Upper middle income", "Thailand" = "Upper middle income", "South Africa" = "Upper middle income",
  "Nigeria" = "Lower middle income", "Denmark" = "High income", "Ireland" = "High income",
  "Slovakia" = "High income", "Iran" = "Upper middle income", "Colombia" = "Upper middle income",
  "Morocco" = "Lower middle income", "Bosnia & Herzegovina" = "Upper middle income", 
  "Tanzania" = "Low income", "Tunisia" = "Lower middle income", "Saudi Arabia" = "High income",
  "Serbia" = "Upper middle income", "Palestinian Territories" = "Lower middle income", 
  "Czechia" = "High income", "Pakistan" = "Lower middle income", "United Arab Emirates" = "High income",
  "Kenya" = "Lower middle income", "Kazakhstan" = "Upper middle income", "Peru" = "Upper middle income",
  "Vietnam" = "Lower middle income", "Latvia" = "High income", "Bangladesh" = "Lower middle income",
  "Indonesia" = "Lower middle income", "Algeria" = "Lower middle income"
)

biblio_final <- biblio_final %>%
  mutate(
    WHO_region = region_manual[Country],
    Income_Group = income_manual[Country]
  )

################################################################################
# 14. Import dataframe of the result of bibliometric analysis ------------------
################################################################################
# This is done in case you hace the excel file of bibliometric analysis and want
# to import it instead of running the code above

# biblio_final <- read_excel("~/ancestryPGx/data/biblio_final.xlsx")

################################################################################
# 15. Bibliometric Barplots and Visualization of bibliometric data -------------
################################################################################
# This section generates visualizations (barplots) of bibliometric data 
# by WHO region, World Bank income group, journal quartile, and open access status.

################################################################################
# 15.1 Load Required Libraries -------------------------------------------------
################################################################################
library(tidyverse)
library(lubridate)
library(ggplot2)
library(forcats)
library(scales)

################################################################################
# 15.2 Prepare Factors and Intervals -------------------------------------------
################################################################################
# Create a 5-year interval column, ordered from most recent to oldest
biblio_final <- biblio_final %>%
  mutate(
    Interval = case_when(
      Pub_date >= 2021 ~ "2021-2025",
      Pub_date >= 2016 ~ "2016-2020",
      Pub_date >= 2011 ~ "2011-2015",
      Pub_date >= 2006 ~ "2006-2010",
      Pub_date >= 2001 ~ "2001-2005",
      Pub_date >= 1996 ~ "1996-2000",
      Pub_date < 1996  ~ "Before 1996",
      TRUE ~ NA_character_
    ),
    Interval = factor(Interval, levels = c(
      "Before 1996", "1996-2000", "2001-2005", "2006-2010",
      "2011-2015", "2016-2020", "2021-2025"
    ))
  )

# Standardize and factorize Access column
biblio_final <- biblio_final %>%
  mutate(
    Access = tolower(Access),
    Access = case_when(
      Access == "open" ~ "Yes",
      Access == "closed" ~ "No",
      is.na(Access) ~ "NA",
      TRUE ~ "NA"
    ),
    Access = factor(Access, levels = c("Yes", "No", "NA"))
  )

# Factorize for plotting
biblio_final$WHO_region        <- fct_infreq(biblio_final$WHO_region)
biblio_final$Income_Group      <- fct_infreq(biblio_final$Income_Group)
biblio_final$sjr_best_quartile <- fct_infreq(biblio_final$sjr_best_quartile)
biblio_final$Access            <- fct_infreq(biblio_final$Access)

################################################################################
# 15.3 Create and Display Barplots ---------------------------------------------
################################################################################

## --- Plot A: Publications per 5-Year Interval by WHO Region ---
p1 <- ggplot(biblio_final, aes(x = Interval, fill = WHO_region)) +
  geom_bar() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  labs(x = "5-Year Interval", y = "Publications", fill = "WHO Region") +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 350, 50)) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linetype = "dashed", linewidth = 0.7),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "European" = "#F25C75",
      "Americas" = "#96C33D",
      "Western Pacific" = "#F9D24A",
      "Eastern Mediterranean" = "#6B4E9B",
      "South-East Asia" = "#2282C5",
      "African" = "#49B583",
      "NA" = "#C8C8C8"
    ),
    na.value = "#C8C8C8"
  )

## --- Plot B: Publications by Income Group ---
p2 <- ggplot(biblio_final, aes(x = Interval, fill = Income_Group)) +
  geom_bar() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  labs(x = "5-Year Interval", y = "Publications", fill = "Income Group") +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 350, 50)) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linetype = "dashed", linewidth = 0.7),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "High income" = "#e3f1ff",
      "Low income" = "#C6DBEF",
      "Lower middle income" = "#6BAED6",
      "Upper middle income" = "#2171B5",
      "NA" = "#C8C8C8"
    ),
    na.value = "#C8C8C8"
  )

## --- Plot C: Publications by Journal Quartile ---
p3 <- ggplot(biblio_final, aes(x = Interval, fill = sjr_best_quartile)) +
  geom_bar() +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  labs(x = "5-Year Interval", y = "Publications", fill = "Journal Quartile") +
  scale_y_continuous(limits = c(0, 400), breaks = seq(0, 350, 50)) +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linetype = "dashed", linewidth = 0.7),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    legend.position = "none"
  ) +
  scale_fill_manual(
    values = c(
      "Q1" = "#FFD700",
      "Q2" = "#4682B4",
      "Q3" = "#FF69B4",
      "Q4" = "#9B59B6",
      "NA" = "#C8C8C8"
    ),
    na.value = "#C8C8C8"
  )

## --- Plot D: Proportion Open Access Publications (stacked bar) ---
p4 <- ggplot(biblio_final, aes(x = Interval, fill = Access)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  labs(x = "5-Year Interval", y = "Proportion", fill = "Open Access") +
  theme_minimal() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey70", linetype = "dashed", linewidth = 0.7),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 18, hjust = -0.08),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("Yes" = "#FF69B4", "No" = "#003366", "NA" = "#C8C8C8")) +
  scale_y_continuous(labels = percent_format())

# Display all plots
print(p1)
print(p2)
print(p3)
print(p4)

################################################################################
# 16. Export, Session Info, and User Guidance ----------------------------------
################################################################################
# This section provides convenient export commands and recommendations for future users.

################################################################################
# 16.1 Export Main Results -----------------------------------------------------
################################################################################
# Save main analysis tables to Excel or CSV if needed

# Uncomment and set your desired paths to export relevant results:
# write_xlsx(data_significant, "results_data_significant.xlsx")
# write_xlsx(data_significant_unique, "results_data_significant_unique.xlsx")
# write_xlsx(top_afr_ordered, "results_top_afr_ordered.xlsx")
# write_xlsx(top_eur_ordered, "results_top_eur_ordered.xlsx")
# write_xlsx(biblio_final, "results_bibliometric_final.xlsx")

################################################################################
# 16.2 Save Plots --------------------------------------------------------------
################################################################################
# Example: Save the bibliometric plots as images
# ggsave("plot_who_region.png", plot = p1, width = 8, height = 6)
# ggsave("plot_income_group.png", plot = p2, width = 8, height = 6)
# ggsave("plot_journal_quartile.png", plot = p3, width = 8, height = 6)
# ggsave("plot_open_access.png", plot = p4, width = 8, height = 6)

################################################################################
# 16.3 Session Info ------------------------------------------------------------
################################################################################
# Print R session info for reproducibility
sessionInfo()

################################################################################
# 16.4 Guidance for Future Users -----------------------------------------------
################################################################################
# - Ensure all input file paths are correct before running the script.
# - All libraries required are loaded at the start of each section; install as needed.
# - API keys and emails for bibliometric analysis must be set for your own credentials.
# - For large datasets or API requests, consider parallel execution or breaking into chunks.
# - All sections are modular: you can run only the sections you need.
# - If you modify the workflow, update section numbers and comments for clarity.