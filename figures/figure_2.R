#script to plot rare and ultra-rare variants in the same figure 

library(data.table)
library(lattice)
library(ggplot2)
library(readxl)

####new "base" regression corrected for total exome wide PTVs####
#plot of URVs and updated RVs (AF < 0.1% + URVs excluded) regression analyses, with total exome-wide PTVs, PCs and sex as covariates

#read in all sets of data 
scz_urv_geneset <- fread("scz_urv_ptv_geneset_estimates_ewptvcov.tsv", sep='\t', header=TRUE, data.table=FALSE)
scz_rv_geneset <- fread("scz_rv_ptv_geneset_estimates_noURV_ewptvcov.tsv", sep='\t', header=TRUE, data.table=FALSE)

#RV####
scz_rv_geneset <- scz_rv_geneset[scz_rv_geneset$geneset %like% "n_RV_PTV", ]
scz_rv_geneset$geneset <- gsub("n_RV_PTV__RVs_|_intervals_bed_tsv", "", scz_rv_geneset$geneset)

scz_rv_geneset$FULL_NAME <- scz_rv_geneset$geneset
scz_rv_geneset$FULL_NAME <- as.character(scz_rv_geneset$FULL_NAME)
str(scz_rv_geneset)

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "exPFC1"] <- "exPFC1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "exPFC2"] <- "exPFC2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "exCA1"] <- "exCA1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "exCA3"] <- "exCA3"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GABA1"] <- "GABA1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GABA2"] <- "GABA2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "exDG"] <- "exDG"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "ASC1"] <- "ASC1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "ASC2"] <- "ASC2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "ODC1"] <- "Oligodendrocytes"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "OPC"] <- "Oligodendrocytes precursor cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "MG"] <- "Microglia"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "NSC"] <- "Neuronal stem cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "END"] <- "Endothelial cells"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PI_genes"] <- "PI Genes"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxexPFC1"] <- "PI x exPFC1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxexPFC2"] <- "PI x exPFC2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxexCA1"] <- "PI x exCA1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxexCA3"] <- "PI x exCA3"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxGABA1"] <- "PI x GABA1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxGABA2"] <- "PI x GABA2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxexDG"] <- "PI x exDG"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxASC1"] <- "PI x ASC1"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxASC2"] <- "PI x ASC2"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxODC1"] <- "PI x Oligodendrocytes"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxOPC"] <- "PI x Oligodendrocytes precursor cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxMG"] <- "PI x Microglia"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxNSC"] <- "PI x Neuronal stem cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxEND"] <- "PI x Endothelial cells"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "astrocytes_ependymal"] <- "Astrocytes/ependymal"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Dopaminergic_Adult"] <- "Dopaminergic adult"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Dopaminergic_Neuroblast"] <- "Dopaminergic neuroblast"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Embryonic_Dopaminergic_Neuron"] <- "Embryonic dopaminergic neuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Embryonic_GABAergic_Neuron"] <- "Embryonic GABAergic neuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Embryonic_midbrain_nucleus_neurons"] <- "Embryonic midbrain nucleus neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "endothelial-mural"] <- "Endothelial mural"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Hypothalamic_Dopaminergic_Neurons"] <- "Hypothalamic dopaminergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Hypothalamic_GABAergic_Neurons"] <- "Hypothalamic GABAergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Hypothalamic_Glutamatergic_Neurons"] <- "Hypothalamic glutamatergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "interneurons"] <- "Interneurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Medium_Spiny_Neuron"] <- "Medium spiny neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "microglia"] <- "Microglia"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Neural_Progenitors"] <- "Neural progenitors"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Neuroblasts"] <- "Neuroblasts"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Oligodendrocyte_Precursor"] <- "Oligodendrocyte precursor"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Oligodendrocytes"] <- "Oligodendrocytes"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Oxytocin_and_Vasopressin_Expressing_Neurons"] <- "Oxytocin/vasopressin neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "pyramidal_CA1"] <- "Pyramidal (CA1)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "pyramidal_SS"] <- "Pyramidal (S1)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Radial_glia_like_cells"] <- "Radial glia-like cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Serotonergic_Neuron"] <- "Serotonergic Neuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Striatal_Interneuron"] <- "Striatal interneuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Vascular_Leptomeningeal_Cells"] <- "Vascular leptomeningeal"


scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxastrocytes_ependymal"] <- "PI x Astrocytes/ependymal"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxDopaminergic_Adult"] <- "PI x Dopaminergic adult"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxDopaminergic_Neuroblast"] <- "PI x Dopaminergic neuroblast"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxEmbryonic_Dopaminergic_Neuron"] <- "PI x Embryonic dopaminergic neuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxEmbryonic_GABAergic_Neuron"] <- "PI x Embryonic GABAergic neuron" 
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxEmbryonic_midbrain_nucleus"] <- "PI x Embryonic midbrain nucleus neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxendothelial_mural"] <- "PI x Endothelial mural"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxHypothalamic_Dopaminergic"] <- "PI x Hypothalamic dopaminergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxHypothalamic_GABAergic"] <- "PI x Hypothalamic GABAergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxHypothalamic_Glutamatergic"] <- "PI x Hypothalamic glutamatergic neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxinterneurons"] <- "PI x Interneurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxMedium_Spiny_Neuron"] <- "PI x Medium spiny neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxmicroglia"] <- "PI x Microglia"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxNeural_Progenitors"] <- "PI x Neural progenitors"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxNeuroblasts"] <- "PI x Neuroblasts"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxOligodendrocyte_Precursor"] <- "PI x Oligodendrocyte precursor"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxOxytocin_and_Vasopressin"] <- "PI x Oxytocin/vasopressin neurons"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxpyramidal_CA1"] <- "PI x Pyramidal (CA1)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxpyramidal_SS"] <- "PI x Pyramidal (S1)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxRadial_glia_like_cells"] <- "PI x Radial glia like cells"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxSerotonergic_Neuron"] <- "PI x Serotonergic Neuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxStriatal_Interneuron"] <- "PI x Striatal interneuron"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxVascular_Leptomeningeal_Cell"] <- "PI x Vascular leptomeningeal"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0007416"] <- "Synapse assembly"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0050808"] <- "Synapse organization"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_synprocess"] <- "Synaptic process"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099536"] <- "Synaptic signaling"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099537"] <- "Trans-synaptic signaling"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0050804"] <- "Modulation of chemical synaptic transmission"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0007268"] <- "Chemical synaptic transmission"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_metabolism"] <- "Metabolism"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_presynprocess"] <- "Presynaptic process"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099504"] <- "Synaptic vesicle cycle"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_postsynprocess"] <- "Postsynaptic process"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099072"] <- "Regulation of postsynaptic membrane neurotransmitter"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099173"] <- "Postsynapse organization"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0060078"] <- "Regulation of postsynaptic membrane potential"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0016079"] <- "Synaptic vesicle exocytosis"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_translation_postsynapse"] <- "Translation postsynapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_translation_synapse"] <- "Translation synapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_translation_presynapse"] <- "Translation presynapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099572"] <- "Postsynaptic specialization"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0098794"] <- "Postsynapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0045202"] <- "Synapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0014069"] <- "Postsynaptic density"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0098793"] <- "Presynapse"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0098839"] <- "Postsynaptic density membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099061"] <- "Integral component of postsynaptic density membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0042734"] <- "Presynaptic membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099056"] <- "Integral component of presynaptic membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0099055"] <- "Integral component of postsynaptic membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0045211"] <- "Postsynaptic membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0048786"] <- "Presynaptic active zone"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0048787"] <- "Presynaptic active zone membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0008021"] <- "Synaptic vesicle"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "GO_0030672"] <- "Synaptic vesicle membrane"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_postsyn_ribosome"] <- "Postsynaptic ribosome"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "SYNGO_presyn_ribosome"] <- "Presynaptic ribosome"

#rename oligodendrocytes and microglia for mouse and humans in rv data
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "ODC1"] <- "Oligodendrocytes (human)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "MG"] <- "Microglia (human)"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxODC1"] <- "PI x Oligodendrocytes (human)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxMG"] <- "PI x Microglia (human)"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "microglia"] <- "Microglia (mouse)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "Oligodendrocytes"] <- "Oligodendrocytes (mouse)"

scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxmicroglia"] <- "PI x Microglia (mouse)"
scz_rv_geneset$FULL_NAME[scz_rv_geneset$geneset == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes (mouse)"

#change colnames to differentiate estimates for variant types
names(scz_rv_geneset)[3] <- "beta_rv"
names(scz_rv_geneset)[4] <- "se_rv"
names(scz_rv_geneset)[5] <- "z_rv"
names(scz_rv_geneset)[6] <- "OR_rv"
names(scz_rv_geneset)[7] <- "p_rv"
names(scz_rv_geneset)[8] <- "fdr_p_rv"
names(scz_rv_geneset)[9] <- "significant_rv"

#URV####
scz_urv_geneset <- scz_urv_geneset[scz_urv_geneset$geneset %like% "n_URV_PTV", ]
scz_urv_geneset$geneset <- gsub("n_URV_PTV__URVs_|_intervals_bed_tsv", "", scz_urv_geneset$geneset)

scz_urv_geneset$FULL_NAME <- scz_urv_geneset$geneset
scz_urv_geneset$FULL_NAME <- as.character(scz_urv_geneset$FULL_NAME)
str(scz_urv_geneset)

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "exPFC1"] <- "exPFC1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "exPFC2"] <- "exPFC2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "exCA1"] <- "exCA1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "exCA3"] <- "exCA3"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GABA1"] <- "GABA1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GABA2"] <- "GABA2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "exDG"] <- "exDG"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "ASC1"] <- "ASC1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "ASC2"] <- "ASC2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "ODC1"] <- "Oligodendrocytes"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "OPC"] <- "Oligodendrocytes precursor cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "MG"] <- "Microglia"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "NSC"] <- "Neuronal stem cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "END"] <- "Endothelial cells"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PI_genes"] <- "PI Genes"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxexPFC1"] <- "PI x exPFC1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxexPFC2"] <- "PI x exPFC2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxexCA1"] <- "PI x exCA1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxexCA3"] <- "PI x exCA3"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxGABA1"] <- "PI x GABA1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxGABA2"] <- "PI x GABA2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxexDG"] <- "PI x exDG"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxASC1"] <- "PI x ASC1"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxASC2"] <- "PI x ASC2"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxODC1"] <- "PI x Oligodendrocytes"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOPC"] <- "PI x Oligodendrocytes precursor cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxMG"] <- "PI x Microglia"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxNSC"] <- "PI x Neuronal stem cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxEND"] <- "PI x Endothelial cells"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "astrocytes_ependymal"] <- "Astrocytes/ependymal"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Dopaminergic_Adult"] <- "Dopaminergic adult"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Dopaminergic_Neuroblast"] <- "Dopaminergic neuroblast"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Embryonic_Dopaminergic_Neuron"] <- "Embryonic dopaminergic neuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Embryonic_GABAergic_Neuron"] <- "Embryonic GABAergic neuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Embryonic_midbrain_nucleus_neurons"] <- "Embryonic midbrain nucleus neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "endothelial-mural"] <- "Endothelial mural"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Hypothalamic_Dopaminergic_Neurons"] <- "Hypothalamic dopaminergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Hypothalamic_GABAergic_Neurons"] <- "Hypothalamic GABAergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Hypothalamic_Glutamatergic_Neurons"] <- "Hypothalamic glutamatergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "interneurons"] <- "Interneurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Medium_Spiny_Neuron"] <- "Medium spiny neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "microglia"] <- "Microglia"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Neural_Progenitors"] <- "Neural progenitors"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Neuroblasts"] <- "Neuroblasts"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Oligodendrocyte_Precursor"] <- "Oligodendrocyte precursor"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Oligodendrocytes"] <- "Oligodendrocytes"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Oxytocin_and_Vasopressin_Expressing_Neurons"] <- "Oxytocin/vasopressin neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "pyramidal_CA1"] <- "Pyramidal (CA1)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "pyramidal_SS"] <- "Pyramidal (S1)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Radial_glia_like_cells"] <- "Radial glia-like cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Serotonergic_Neuron"] <- "Serotonergic Neuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Striatal_Interneuron"] <- "Striatal interneuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Vascular_Leptomeningeal_Cells"] <- "Vascular leptomeningeal"


scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxastrocytes_ependymal"] <- "PI x Astrocytes/ependymal"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxDopaminergic_Adult"] <- "PI x Dopaminergic adult"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxDopaminergic_Neuroblast"] <- "PI x Dopaminergic neuroblast"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxEmbryonic_Dopaminergic_Neuron"] <- "PI x Embryonic dopaminergic neuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxEmbryonic_GABAergic_Neuron"] <- "PI x Embryonic GABAergic neuron" 
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxEmbryonic_midbrain_nucleus"] <- "PI x Embryonic midbrain nucleus neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxendothelial_mural"] <- "PI x Endothelial mural"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxHypothalamic_Dopaminergic"] <- "PI x Hypothalamic dopaminergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxHypothalamic_GABAergic"] <- "PI x Hypothalamic GABAergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxHypothalamic_Glutamatergic"] <- "PI x Hypothalamic glutamatergic neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxinterneurons"] <- "PI x Interneurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxMedium_Spiny_Neuron"] <- "PI x Medium spiny neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxmicroglia"] <- "PI x Microglia"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxNeural_Progenitors"] <- "PI x Neural progenitors"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxNeuroblasts"] <- "PI x Neuroblasts"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOligodendrocyte_Precursor"] <- "PI x Oligodendrocyte precursor"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOxytocin_and_Vasopressin"] <- "PI x Oxytocin/vasopressin neurons"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxpyramidal_CA1"] <- "PI x Pyramidal (CA1)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxpyramidal_SS"] <- "PI x Pyramidal (S1)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxRadial_glia_like_cells"] <- "PI x Radial glia like cells"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxSerotonergic_Neuron"] <- "PI x Serotonergic Neuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxStriatal_Interneuron"] <- "PI x Striatal interneuron"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxVascular_Leptomeningeal_Cell"] <- "PI x Vascular leptomeningeal"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0007416"] <- "Synapse assembly"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0050808"] <- "Synapse organization"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_synprocess"] <- "Synaptic process"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099536"] <- "Synaptic signaling"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099537"] <- "Trans-synaptic signaling"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0050804"] <- "Modulation of chemical synaptic transmission"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0007268"] <- "Chemical synaptic transmission"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_metabolism"] <- "Metabolism"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_presynprocess"] <- "Presynaptic process"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099504"] <- "Synaptic vesicle cycle"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_postsynprocess"] <- "Postsynaptic process"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099072"] <- "Regulation of postsynaptic membrane neurotransmitter"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099173"] <- "Postsynapse organization"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0060078"] <- "Regulation of postsynaptic membrane potential"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0016079"] <- "Synaptic vesicle exocytosis"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_translation_postsynapse"] <- "Translation postsynapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_translation_synapse"] <- "Translation synapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_translation_presynapse"] <- "Translation presynapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099572"] <- "Postsynaptic specialization"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0098794"] <- "Postsynapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0045202"] <- "Synapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0014069"] <- "Postsynaptic density"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0098793"] <- "Presynapse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0098839"] <- "Postsynaptic density membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099061"] <- "Integral component of postsynaptic density membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0042734"] <- "Presynaptic membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099056"] <- "Integral component of presynaptic membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0099055"] <- "Integral component of postsynaptic membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0045211"] <- "Postsynaptic membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0048786"] <- "Presynaptic active zone"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0048787"] <- "Presynaptic active zone membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0008021"] <- "Synaptic vesicle"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "GO_0030672"] <- "Synaptic vesicle membrane"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_postsyn_ribosome"] <- "Postsynaptic ribosome"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "SYNGO_presyn_ribosome"] <- "Presynaptic ribosome"

#rename oligodendrocytes and microglia for mouse and humans in urv data
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "ODC1"] <- "Oligodendrocytes (human)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "MG"] <- "Microglia (human)"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxODC1"] <- "PI x Oligodendrocytes (human)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxMG"] <- "PI x Microglia (human)"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "microglia"] <- "Microglia (mouse)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Oligodendrocytes"] <- "Oligodendrocytes (mouse)"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxmicroglia"] <- "PI x Microglia (mouse)"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes (mouse)"

#change colnames to differentiate estimates for variant types
names(scz_urv_geneset)[3] <- "beta_urv"
names(scz_urv_geneset)[4] <- "se_urv"
names(scz_urv_geneset)[5] <- "z_urv"
names(scz_urv_geneset)[6] <- "OR_urv"
names(scz_urv_geneset)[7] <- "p_urv"
names(scz_urv_geneset)[8] <- "fdr_p_urv"
names(scz_urv_geneset)[9] <- "significant_urv"

####group gene sets in urv data####
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
  "Embryonic_midbrain_nucleus_neurons",
  "endothelial-mural",
  "Hypothalamic_Dopaminergic_Neurons",
  "Hypothalamic_GABAergic_Neurons",
  "Hypothalamic_Glutamatergic_Neurons",
  "interneurons",
  "Medium_Spiny_Neuron",
  "microglia",
  "Neural_Progenitors",
  "Neuroblasts",
  "Oligodendrocyte_Precursor",
  "Oligodendrocytes",
  "Oxytocin_and_Vasopressin_Expressing_Neurons",
  "pyramidal_CA1",
  "pyramidal_SS",
  "Radial_glia_like_cells",
  "Serotonergic_Neuron",
  "Striatal_Interneuron",
  "Vascular_Leptomeningeal_Cells")
pix_mouse <- c("PIxastrocytes_ependymal",
               "PIxDopaminergic_Adult",
               "PIxDopaminergic_Neuroblast",
               "PIxEmbryonic_Dopaminergic_Neuron",
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
               "PIxVascular_Leptomeningeal_Cell")
syngo_cc <- c("GO_0045202",	"GO_0098793",	"GO_0098794",	"GO_0048786",	"GO_0008021",	"GO_0042734",	"SYNGO_presyn_ribosome",	"GO_0099572",	"GO_0045211",	
              "SYNGO_postsyn_ribosome",	
              "GO_0048787",	"GO_0030672",	"GO_0099056",	"GO_0014069",	"GO_0099055",	"GO_0098839",	"GO_0099061")
syngo_bp <- c("SYNGO_synprocess",	"SYNGO_presynprocess",	"SYNGO_postsynprocess",	"GO_0099536",	"GO_0050808",	"SYNGO_metabolism",	"GO_0099504",	"GO_0060078",	
              "GO_0099072",	"GO_0099537",	"GO_0099173",	"GO_0007416",	"SYNGO_translation_synapse",	"GO_0016079",	"GO_0007268",	"GO_0099175",	
              "SYNGO_translation_presynapse",	
              "SYNGO_translation_postsynapse",	"GO_0050804")


scz_urv_geneset$category[scz_urv_geneset$geneset %in% human_brain_cells] <- "Human brain cells"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% pi_genes] <- "PI"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% pix_human] <- "PI x human brain cells"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% mouse_brain_cells] <- "Mouse brain cells"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% pix_mouse] <- "PI x mouse brain cells"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% syngo_cc] <- "SynGo cellular components"
scz_urv_geneset$category[scz_urv_geneset$geneset %in% syngo_bp] <- "SynGo biological processes"

#merge all three dataframes on full names
merged_genesets <- merge(scz_urv_geneset,scz_rv_geneset, by="FULL_NAME")

#add ngenes column
ngenes <- read_excel("H:/first_projects/Exome sequencing/03_data_analyses/genesets_intervals/geneset_size_urvs.xlsx")
merged_genesets <- merge(merged_genesets, ngenes, by="geneset.x")
merged_genesets$NEW_NAME <- paste0(merged_genesets$FULL_NAME, " (",paste0(merged_genesets$NGENES), ")")

merged_genesets$category <- factor(merged_genesets$category, levels = c("PI", "Human brain cells", "PI x human brain cells", "Mouse brain cells", "PI x mouse brain cells", "SynGo cellular components","SynGo biological processes"))

#to exclude syn_go stuff 
merged_genesets_syngo <- subset(merged_genesets, merged_genesets$category == "SynGo cellular components" | merged_genesets$category == "SynGo biological processes")
merged_genesets <- subset(merged_genesets, merged_genesets$category != "SynGo cellular components" & merged_genesets$category != "SynGo biological processes")

#plot results
theme_set(theme_light())

nudge1 <- position_nudge(x = .1, y = 0)
nudge2 <- position_nudge(x = -.1, y = 0)

#PI and brain gene sets 
tiff('figure_2.tiff',width=5500,height=8000,res=600)
ggplot(merged_genesets, aes(y=beta_urv, x=NEW_NAME)) + 
  scale_colour_manual( name="Data", values=c("1" = "#E7B800","2" = "#CC79A7"), labels = c("Rare variants", "Ultra-rare variants")) + 
  
  geom_point(aes(y=beta_rv, x=NEW_NAME,color="1"), size=2.5, position = nudge1) +
  geom_errorbar(aes(ymin=beta_rv-1.96*se_rv, ymax=beta_rv+1.96*se_rv,color="1"), width=.3, size = 0.8, show.legend=FALSE, position = nudge1) +
  geom_point(aes(y=beta_rv, x=NEW_NAME), data = merged_genesets[merged_genesets$fdr_p_rv < .05, ], color="black", shape = 8, size=1.3, show.legend = F, position = nudge1) +
  
  geom_point(aes(y=beta_urv, x=NEW_NAME,color="2"), size=2.5, position = nudge2) +
  geom_errorbar(aes(ymin=beta_urv-1.96*se_urv, ymax=beta_urv+1.96*se_urv,color="2"), width=.3, size = 0.8, show.legend=FALSE, position = nudge2) +
  geom_point(aes(y=beta_urv, x=NEW_NAME), data = merged_genesets[merged_genesets$fdr_p_urv < .05, ], color="black", shape = 8, size=1.3, show.legend = F, position = nudge2) +
  
  theme(strip.text.y = element_text(size = 12, face = "bold")) +
  facet_grid(category ~ ., scales = "free", space = "free_y") +
  
  ylab("Beta") +
  
  geom_hline(yintercept = 0) +
  
  #theme_light() + #must be turned off to change facet label xteristics
  theme(axis.title.y=element_blank(),
        legend.title=element_blank(), 
        axis.title.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y=element_text(size=9.5, face = "bold")) +
  theme(legend.position="right") +
  theme(legend.text = element_text(size = 11)) +
  coord_flip() 

dev.off()

#syngo gene sets 
tiff('figure_2_syngo.tiff',width=5500,height=8000,res=600)
ggplot(merged_genesets_syngo, aes(y=beta_urv, x=NEW_NAME)) + 
  scale_colour_manual( name="Data", values=c("1" = "#E7B800","2" = "#CC79A7"), labels = c("Rare variants", "Ultra-rare variants")) + 
  
  geom_point(aes(y=beta_rv, x=NEW_NAME,color="1"), size=2.5, position = nudge1) +
  geom_errorbar(aes(ymin=beta_rv-1.96*se_rv, ymax=beta_rv+1.96*se_rv,color="1"), width=.3, size = 0.8, show.legend=FALSE, position = nudge1) +
  geom_point(aes(y=beta_rv, x=NEW_NAME), data = merged_genesets_syngo[merged_genesets_syngo$fdr_p_rv < .05, ], color="black", shape = 8, size=1.3, show.legend = F, position = nudge1) +
  
  geom_point(aes(y=beta_urv, x=NEW_NAME,color="2"), size=2.5, position = nudge2) +
  geom_errorbar(aes(ymin=beta_urv-1.96*se_urv, ymax=beta_urv+1.96*se_urv,color="2"), width=.3, size = 0.8, show.legend=FALSE, position = nudge2) +
  geom_point(aes(y=beta_urv, x=NEW_NAME), data = merged_genesets_syngo[merged_genesets_syngo$fdr_p_urv < .05, ], color="black", shape = 8, size=1.3, show.legend = F, position = nudge2) +
  
  theme(strip.text.y = element_text(size = 12, face = "bold")) +
  facet_grid(category ~ ., scales = "free", space = "free_y") +
  
  ylab("Beta") +
  
  geom_hline(yintercept = 0) +
  
  #theme_light() + #must be turned off to change facet label xteristics
  theme(axis.title.y=element_blank(),
        legend.title=element_blank(), 
        axis.title.x = element_text(size = 11, face = "bold")) +
  theme(axis.text.y=element_text(size=9.5, face = "bold")) +
  theme(legend.position="right") +
  theme(legend.text = element_text(size = 11)) +
  coord_flip() 

dev.off()

