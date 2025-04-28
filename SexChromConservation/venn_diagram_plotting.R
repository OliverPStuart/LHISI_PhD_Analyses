### Plotting venn diagram of sex chromosome gene content

# Environment
HOME_DIR="/Volumes/Alter/LHISI"
WORKING_DIR=paste0(HOME_DIR,"/Analyses/SexChromConservation")
setwd(WORKING_DIR)

library(ggvenn)

# Function make venn diagram
# Taken from https://github.com/DarrenJParker/Timema_sex_chr_evol_code/blob/main/Orthologs_on_the_X.R#L48

five_sp_DE_venn <- function(Samer,Sgreg,Spice,Chookeri,LHISI){
    venny.plot <- venn.diagram(
      list("Samer" = Samer, "Sgreg" = Sgreg, "Spice" = Spice, "Chookeri" = Chookeri, "LHISI" = LHISI ), filename = NULL,
      fill = c("chocolate", "gold", "goldenrod1", "olivedrab3", "palegreen3"),
      cat.col = c("chocolate", "gold", "goldenrod1", "olivedrab3", "palegreen3"),
      margin = 0.2, cat.dist = 0.3, cat.cex = 2)
    return(venny.plot)
  }
  
# Sharing of X chromosome hits

sex_sharing <- five_sp_DE_venn(Samer=readLines("Samer_sex_hits"),
                               Sgreg=readLines("Sgreg_sex_hits"),
                               Spice=readLines("Spice_sex_hits"),
                               Chookeri=readLines("Chookeri_sex_hits"),
                               LHISI=readLines("LHISI_sex_hits")
                               )

ggsave(file="sex_chromosome_hit_sharing.pdf", sex_sharing, width=9,height=9,units='in')

# Sharing of Phasmid X chromosome hits with the other big hit chromosome in schistocerca

sex_sharing <- five_sp_DE_venn(Samer=readLines("Samer_supp_hits"),
                               Sgreg=readLines("Sgreg_supp_hits"),
                               Spice=readLines("Spice_supp_hits"),
                               Chookeri=readLines("Chookeri_sex_hits"),
                               LHISI=readLines("LHISI_sex_hits")
)

ggsave(file="sex_chromosome_supp_hit_sharing.pdf", sex_sharing, width=9,height=9,units='in')

# Sharing of phasmid X chromosomes with schistocerca X + alternative chromosome

sex_sharing <- five_sp_DE_venn(Samer=readLines("Samer_both_hits"),
                               Sgreg=readLines("Sgreg_both_hits"),
                               Spice=readLines("Spice_both_hits"),
                               Chookeri=readLines("Chookeri_sex_hits"),
                               LHISI=readLines("LHISI_sex_hits")
)

ggsave(file="sex_chromosome_both_hit_sharing.pdf", sex_sharing, width=9,height=9,units='in')
