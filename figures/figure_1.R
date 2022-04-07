#read in magma results
SCZ_geneset <- read.table("SCZ_complete_genesets_transpose.gmt.gsa.out", header = T)

SCZ_geneset$FULL_NAME <- as.character(SCZ_geneset$FULL_NAME)
str(SCZ_geneset)

SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "exPFC1"] <- "exPFC1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "exPFC2"] <- "exPFC2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "exCA1"] <- "exCA1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "exCA3"] <- "exCA3"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GABA1"] <- "GABA1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GABA2"] <- "GABA2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "exDG"] <- "exDG"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "ASC1"] <- "ASC1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "ASC2"] <- "ASC2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "ODC1"] <- "Oligodendrocytes"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "OPC"] <- "Oligodendrocytes precursor cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "MG"] <- "Microglia"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "NSC"] <- "Neuronal stem cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "END"] <- "Endothelial cells"

SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PI_genes"] <- "PI Genes"

SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxexPFC1"] <- "PI x exPFC1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxexPFC2"] <- "PI x exPFC2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxexCA1"] <- "PI x exCA1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxexCA3"] <- "PI x exCA3"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxGABA1"] <- "PI x GABA1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxGABA2"] <- "PI x GABA2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxexDG"] <- "PI x exDG"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxASC1"] <- "PI x ASC1"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxASC2"] <- "PI x ASC2"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxODC1"] <- "PI x Oligodendrocytes"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxOPC"] <- "PI x Oligodendrocytes precursor cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxMG"] <- "PI x Microglia"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxNSC"] <- "PI x Neuronal stem cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxEND"] <- "PI x Endothelial cells"

SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "astrocytes_ependymal"] <- "Astrocytes/ependymal"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Dopaminergic_Adult"] <- "Dopaminergic adult"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Dopaminergic_Neuroblast"] <- "Dopaminergic neuroblast"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Embryonic_Dopaminergic_Neuron"] <- "Embryonic dopaminergic neuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Embryonic_GABAergic_Neuron"] <- "Embryonic GABAergic neuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Embryonic_midbrain_nucleus_n..."] <- "Embryonic midbrain nucleus neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "endothelial-mural"] <- "Endothelial-mural"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Hypothalamic_Dopaminergic_Ne..."] <- "Hypothalamic dopaminergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Hypothalamic_GABAergic_Neuro..."] <- "Hypothalamic GABAergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Hypothalamic_Glutamatergic_N..."] <- "Hypothalamic glutamatergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "interneurons"] <- "Interneurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Medium_Spiny_Neuron"] <- "Medium spiny neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "microglia"] <- "Microglia"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Neural_Progenitors"] <- "Neural progenitors"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Neuroblasts"] <- "Neuroblasts"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Oligodendrocyte_Precursor"] <- "Oligodendrocyte precursor"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Oligodendrocytes"] <- "Oligodendrocytes"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Oxytocin_and_Vasopressin_Exp..."] <- "Oxytocin/vasopressin neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "pyramidal_CA1"] <- "Pyramidal (CA1)"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "pyramidal_SS"] <- "Pyramidal (S1)"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Radial_glia_like_cells"] <- "Radial glia-like cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Serotonergic_Neuron"] <- "Serotonergic Neuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Striatal_Interneuron"] <- "Striatal interneuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "Vascular_Leptomeningeal_Cells"] <- "Vascular leptomeningeal"


SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxastrocytes_ependymal"] <- "PI x Astrocytes/ependymal"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxDopaminergic_Adult"] <- "PI x Dopaminergic adult"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxDopaminergic_Neuroblast"] <- "PI x Dopaminergic neuroblast"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxEmbryonic_Dopaminergic_Ne..."] <- "PI x Embryonic dopaminergic neuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxEmbryonic_GABAergic_Neuron"] <- "PI x Embryonic GABAergic neuron" 
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxEmbryonic_midbrain_nucleus"] <- "PI x Embryonic midbrain nucleus neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxendothelial_mural"] <- "PI x Endothelial-mural"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxHypothalamic_Dopaminergic"] <- "PI x Hypothalamic dopaminergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxHypothalamic_GABAergic"] <- "PI x Hypothalamic GABAergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxHypothalamic_Glutamatergic"] <- "PI x Hypothalamic glutamatergic neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxinterneurons"] <- "PI x Interneurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxMedium_Spiny_Neuron"] <- "PI x Medium spiny neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxmicroglia"] <- "PI x Microglia"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxNeural_Progenitors"] <- "PI x Neural progenitors"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxNeuroblasts"] <- "PI x Neuroblasts"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxOligodendrocyte_Precursor"] <- "PI x Oligodendrocyte precursor"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxOxytocin_and_Vasopressin"] <- "PI x Oxytocin/vasopressin neurons"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxpyramidal_CA1"] <- "PI x Pyramidal (CA1)"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxpyramidal_SS"] <- "PI x Pyramidal (S1)"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxRadial_glia_like_cells"] <- "PI x Radial glia-like cells"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxSerotonergic_Neuron"] <- "PI x Serotonergic Neuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxStriatal_Interneuron"] <- "PI x Striatal interneuron"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "PIxVascular_Leptomeningeal_C..."] <- "PI x Vascular leptomeningeal"

SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0007416"] <- "Synapse assembly"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0050808"] <- "Synapse organization"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:synprocess"] <- "Synaptic process"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099536"] <- "Synaptic signaling"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099537"] <- "Trans-synaptic signaling"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0050804"] <- "Modulation of chemical synaptic transmission"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0007268"] <- "Chemical synaptic transmission"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:metabolism"] <- "Metabolism"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:presynprocess"] <- "Presynaptic process"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099504"] <- "Synaptic vesicle cycle"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:postsynprocess"] <- "Postsynaptic process"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099072"] <- "Regulation of postsynaptic membrane neurotransmitter"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099173"] <- "Postsynapse organization"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0060078"] <- "Regulation of postsynaptic membrane potential"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0016079"] <- "Synaptic vesicle exocytosis"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:translation_postsynapse"] <- "Translation postsynapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:translation_synapse"] <- "Translation synapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:translation_presynapse"] <- "Translation presynapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099572"] <- "Postsynaptic specialization"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0098794"] <- "Postsynapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0045202"] <- "Synapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0014069"] <- "Postsynaptic density"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0098793"] <- "Presynapse"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0098839"] <- "Postsynaptic density membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099061"] <- "Integral component of postsynaptic density membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0042734"] <- "Presynaptic membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099056"] <- "Integral component of presynaptic membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0099055"] <- "Integral component of postsynaptic membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0045211"] <- "Postsynaptic membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0048786"] <- "Presynaptic active zone"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0048787"] <- "Presynaptic active zone membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0008021"] <- "Synaptic vesicle"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "GO:0030672"] <- "Synaptic vesicle membrane"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:postsyn_ribosome"] <- "Postsynaptic ribosome"
SCZ_geneset$FULL_NAME[SCZ_geneset$VARIABLE == "SYNGO:presyn_ribosome"] <- "Presynaptic ribosome"
#SCZ_geneset[113,] <- SCZ_geneset[1,]

####group gene sets####
human_brain_cells <- c("exPFC1",
                       "exPFC2",
                       "exCA1",
                       "exCA3",
                       "GABA1",
                       "GABA2",
                       "exDG",
                       "ASC1",
                       "ASC2",
                       "ODC1",
                       "OPC",
                       "MG",
                       "NSC",
                       "END")
pi_genes <- c("PI_genes")
pix_human <- c("PIxexPFC1",
               "PIxexPFC2",
               "PIxexCA1",
               "PIxexCA3",
               "PIxGABA1",
               "PIxGABA2",
               "PIxexDG",
               "PIxASC1",
               "PIxASC2",
               "PIxODC1",
               "PIxOPC",
               "PIxMG",
               "PIxNSC",
               "PIxEND")
mouse_brain_cells <- c(
"astrocytes_ependymal",
"Dopaminergic_Adult",
"Dopaminergic_Neuroblast",
"Embryonic_Dopaminergic_Neuron",
"Embryonic_GABAergic_Neuron",
"Embryonic_midbrain_nucleus_n...",
"endothelial-mural",
"Hypothalamic_Dopaminergic_Ne...",
"Hypothalamic_GABAergic_Neuro...",
"Hypothalamic_Glutamatergic_N...",
"interneurons",
"Medium_Spiny_Neuron",
"microglia",
"Neural_Progenitors",
"Neuroblasts",
"Oligodendrocyte_Precursor",
"Oligodendrocytes",
"Oxytocin_and_Vasopressin_Exp...",
"pyramidal_CA1",
"pyramidal_SS",
"Radial_glia_like_cells",
"Serotonergic_Neuron",
"Striatal_Interneuron",
"Vascular_Leptomeningeal_Cells")
pix_mouse <- c("PIxastrocytes_ependymal",
               "PIxDopaminergic_Adult",
               "PIxDopaminergic_Neuroblast",
               "PIxEmbryonic_Dopaminergic_Ne...",
               "PIxEmbryonic_GABAergic_Neuron",
               "PIxEmbryonic_midbrain_nucleus",
               "PIxendothelial_mural",
               "PIxHypothalamic_Dopaminergic",
               "PIxHypothalamic_GABAergic",
               "PIxHypothalamic_Glutamatergic",
               "PIxinterneurons",
               "PIxMedium_Spiny_Neuron",
               "PIxmicroglia",
               "PIxNeural_Progenitors",
               "PIxNeuroblasts",
               "PIxOligodendrocyte_Precursor",
               "PIxOligodendrocytes",
               "PIxOxytocin_and_Vasopressin",
               "PIxpyramidal_CA1",
               "PIxpyramidal_SS",
               "PIxRadial_glia_like_cells",
               "PIxSerotonergic_Neuron",
               "PIxStriatal_Interneuron",
               "PIxVascular_Leptomeningeal_C...")
syngo_cc <- c("GO:0045202",	"GO:0098793",	"GO:0098794",	"GO:0048786",	"GO:0008021",	"GO:0042734",	"SYNGO:presyn_ribosome",	"GO:0099572",	"GO:0045211",	
              "SYNGO:postsyn_ribosome",	
              "GO:0048787",	"GO:0030672",	"GO:0099056",	"GO:0014069",	"GO:0099055",	"GO:0098839",	"GO:0099061")
syngo_bp <- c("SYNGO:synprocess",	"SYNGO:presynprocess",	"SYNGO:postsynprocess",	"GO:0099536",	"GO:0050808",	"SYNGO:metabolism",	"GO:0099504",	"GO:0060078",	
              "GO:0099072",	"GO:0099537",	"GO:0099173",	"GO:0007416",	"SYNGO:translation_synapse",	"GO:0016079",	"GO:0007268",	"GO:0099175",	
              "SYNGO:translation_presynapse",	
              "SYNGO:translation_postsynapse",	"GO:0050804")


SCZ_geneset$category[SCZ_geneset$VARIABLE %in% human_brain_cells] <- "Human brain cells"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% pi_genes] <- "PI"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% pix_human] <- "PI x human brain cells"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% mouse_brain_cells] <- "Mouse brain cells"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% pix_mouse] <- "PI x mouse brain cells"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% syngo_cc] <- "SynGo cellular components"
SCZ_geneset$category[SCZ_geneset$VARIABLE %in% syngo_bp] <- "SynGo biological processes"
SCZ_geneset$category <- factor(SCZ_geneset$category, levels = c("PI", "Human brain cells", "PI x human brain cells", "Mouse brain cells", "PI x mouse brain cells", "SynGo cellular components","SynGo biological processes"))
#SCZ_geneset$category <- factor(SCZ_geneset$category, levels = c("PI", "Human brain cells", "PI x human genes"))

SCZ_geneset <- SCZ_geneset[complete.cases(SCZ_geneset), ]

#add fdr_p values  
SCZ_geneset$fdr_p_common <- p.adjust(SCZ_geneset$P, method = "fdr")
SCZ_geneset$significant_common <- ifelse(SCZ_geneset$fdr_p_common < 0.05, "yes", "no")

#to exclude syn_go results 
SCZ_geneset_syngo <- subset(SCZ_geneset, SCZ_geneset$category == "SynGo cellular components" | SCZ_geneset$category == "SynGo biological processes")
SCZ_geneset <- subset(SCZ_geneset, SCZ_geneset$category != "SynGo cellular components" & SCZ_geneset$category != "SynGo biological processes")

####plot####
library(lattice)
library(ggplot2)
theme_set(theme_light())

#plot to match (U)RV plot - ngenes in brackets next to gene set names 
tiff('figure_1.tiff',width=4000,height=6500,res=600)
ggplot(SCZ_geneset, aes(y=BETA, x=paste0(FULL_NAME, " (",paste0(NGENES), ")"))) + 
  geom_point(color="#00AFBB", size=2.5) +
  geom_errorbar(aes(ymin=BETA-1.96*SE, ymax=BETA+1.96*SE), width=.2,color="#00AFBB") +
  geom_point(data = SCZ_geneset[SCZ_geneset$fdr_p_common < .05, ], color="black", shape = 8, size=1.3, show.legend = F) +
  
  theme(strip.text.y = element_text(size = 12, face = "bold")) +
  facet_grid(category ~ ., scales = "free", space = "free_y") +
  
  ylab("Beta") +
  
  geom_hline(yintercept = 0) +
  
  #theme_light() + #must be turned off to change facet label xteristics
  theme(axis.title.y=element_blank(),
        legend.title=element_blank(), 
        axis.title.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y=element_text(size=9.5, face = "bold")) +
  #theme(panel.background = element_rect(fill = NA)) +
  #theme(panel.grid.major = element_line(colour = "light grey"), panel.grid.minor = element_line(colour = "light grey")) +
  #theme(panel.border = element_rect(colour = "light grey", fill = NA)) +
  #theme(axis.ticks = element_line(colour = "light grey")) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid 
  coord_flip() 

dev.off()

#syngo
tiff('figure_1_syngo.tiff',width=4000,height=6500,res=600)
ggplot(SCZ_geneset_syngo, aes(y=BETA, x=paste0(FULL_NAME, " (",paste0(NGENES), ")"))) + 
  geom_point(color="#00AFBB", size=2.5) +
  geom_errorbar(aes(ymin=BETA-1.96*SE, ymax=BETA+1.96*SE), width=.2,color="#00AFBB") +
  geom_point(data = SCZ_geneset_syngo[SCZ_geneset_syngo$fdr_p_common < .05, ], color="black", shape = 8, size=1.3, show.legend = F) +
  
  theme(strip.text.y = element_text(size = 12, face = "bold")) +
  facet_grid(category ~ ., scales = "free", space = "free_y") +
  
  ylab("Beta") +
  
  geom_hline(yintercept = 0) +
  
  #theme_light() + #must be turned off to change facet label xteristics
  theme(axis.title.y=element_blank(),
        legend.title=element_blank(), 
        axis.title.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y=element_text(size=9.5, face = "bold")) +
  #theme(panel.background = element_rect(fill = NA)) +
  #theme(panel.grid.major = element_line(colour = "light grey"), panel.grid.minor = element_line(colour = "light grey")) +
  #theme(panel.border = element_rect(colour = "light grey", fill = NA)) +
  #theme(axis.ticks = element_line(colour = "light grey")) +
  #theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # remove grid 
  coord_flip() 

dev.off()
