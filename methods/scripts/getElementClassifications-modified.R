# Load necessary libraries
library(data.table)  # Efficient data manipulation
library(dplyr)       # Data manipulation
library(tidyr)       # Data tidying
library(stringr)     # String manipulation

# Define custom operator
`%notin%` <- Negate(`%in%`)

# Read the CSV file
mobileOGs <- fread("D:/桌面/文件/博一下/文章2/beatrix-1-6_v1_all/mobileOG-db-beatrix-1.6-All.csv")

# Select relevant columns
mobileOGs_classes <- mobileOGs %>%
  select(`mobileOG Entry Name`, `Major mobileOG Category`, `Minor mobileOG Categories`, ACLAME, immedb, GPD, pVOG, ICE, AICE, CIME, IME, ISFinder, PlasmidRefSeq, COMPASS)

# Classify entries
TE <- mobileOGs_classes %>%
  filter(ISFinder != 0) %>%
  mutate(MGE_Class = "Insertion Sequence (IS)")

IGE <- mobileOGs_classes %>%
  filter((AICE != 0 | ICE != 0 | CIME != 0 | immedb != 0 | IME != 0) & ISFinder == 0) %>%
  mutate(MGE_Class = "Integrative Element (IGE)")

Ph <- mobileOGs_classes %>%
  filter(GPD != 0 | pVOG != 0) %>%
  mutate(MGE_Class = "Phage")

plasmidic <- mobileOGs_classes %>%
  filter((COMPASS != 0 | PlasmidRefSeq != 0) & `mobileOG Entry Name` %notin% IGE$`mobileOG Entry Name` & ISFinder == 0) %>%
  mutate(MGE_Class = "Plasmid")

CE <- mobileOGs_classes %>%
  filter(`Major mobileOG Category` == "transfer" & str_detect(`Minor mobileOG Categories`, "conjugati")) %>%
  mutate(MGE_Class = "Conjugative Element (CE)")

# Combine classified data
merged <- bind_rows(TE, IGE, Ph, plasmidic, CE)

# Find entries only detected in ACLAME
ACLAME <- mobileOGs_classes %>%
  filter(`mobileOG Entry Name` %notin% merged$`mobileOG Entry Name`)

# Further classification within ACLAME
Ph1 <- ACLAME %>%
  filter(`Major mobileOG Category` == "phage") %>%
  mutate(MGE_Class = "Phage")

IGEE1 <- ACLAME %>%
  filter(`Major mobileOG Category` == "integration/excision") %>%
  mutate(MGE_Class = "Integrative Element (IGE)")

# Add further classified ACLAME entries to ACLAME
ACLAME_classified <- bind_rows(Ph1, IGEE1)

# Entries not classified by other categories
ACLAME_only <- ACLAME %>%
  filter(`mobileOG Entry Name` %notin% ACLAME_classified$`mobileOG Entry Name`) %>%
  mutate(MGE_Class = "Others")

# Combine all classified entries
merged <- bind_rows(merged, ACLAME_classified, ACLAME_only)

# Aggregate MGE classes for each entry
out <- merged %>%
  group_by(`mobileOG Entry Name`) %>%
  summarize(MGE_Class = paste(unique(MGE_Class), collapse = ','), .groups = 'drop')

# Write the result to a CSV file
write.csv(out, "D:/桌面/文件/博一下/文章2/beatrix-1-6_v1_all/mobileOG_class_information-rev.csv", row.names = FALSE)

