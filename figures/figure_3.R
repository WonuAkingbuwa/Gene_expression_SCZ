library(wCorr)
library(data.table)

#script to plot beta's of common vs rare variants 
scz_urv_geneset <- fread("scz_urv_ptv_geneset_estimates_ewptvcov.tsv", sep='\t', header=TRUE, data.table=FALSE)
scz_common_geneset <- read.table("SCZ_complete_genesets_transpose.gmt.gsa.out", header = T)

#rename common and rare variant genesets so they match 
scz_common_geneset$FULL_NAME <- as.character(scz_common_geneset$FULL_NAME)
str(scz_common_geneset)

scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "exPFC1"] <- "exPFC1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "exPFC2"] <- "exPFC2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "exCA1"] <- "exCA1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "exCA3"] <- "exCA3"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GABA1"] <- "GABA1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GABA2"] <- "GABA2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "exDG"] <- "exDG"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "ASC1"] <- "ASC1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "ASC2"] <- "ASC2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "ODC1"] <- "Oligodendrocytes_human"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "OPC"] <- "Oligodendrocytes precursor cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "MG"] <- "Microglia_human"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "NSC"] <- "Neuronal stem cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "END"] <- "Endothelial cells"

scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PI_genes"] <- "PI Genes"

scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxexPFC1"] <- "PI x exPFC1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxexPFC2"] <- "PI x exPFC2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxexCA1"] <- "PI x exCA1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxexCA3"] <- "PI x exCA3"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxGABA1"] <- "PI x GABA1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxGABA2"] <- "PI x GABA2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxexDG"] <- "PI x exDG"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxASC1"] <- "PI x ASC1"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxASC2"] <- "PI x ASC2"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxODC1"] <- "PI x Oligodendrocytes_human"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxOPC"] <- "PI x Oligodendrocytes precursor cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxMG"] <- "PI x Microglia_human"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxNSC"] <- "PI x Neuronal stem cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxEND"] <- "PI x Endothelial cells"

scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "astrocytes_ependymal"] <- "Astrocytes/ependymal"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Dopaminergic_Adult"] <- "Dopaminergic adult"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Dopaminergic_Neuroblast"] <- "Dopaminergic neuroblast"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Embryonic_Dopaminergic_Neuron"] <- "Embryonic dopaminergic neuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Embryonic_GABAergic_Neuron"] <- "Embryonic GABAergic neuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Embryonic_midbrain_nucleus_n..."] <- "Embryonic midbrain nucleus neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "endothelial-mural"] <- "Endothelial-mural"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Hypothalamic_Dopaminergic_Ne..."] <- "Hypothalamic dopaminergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Hypothalamic_GABAergic_Neuro..."] <- "Hypothalamic GABAergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Hypothalamic_Glutamatergic_N..."] <- "Hypothalamic glutamatergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "interneurons"] <- "Interneurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Medium_Spiny_Neuron"] <- "Medium spiny neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "microglia"] <- "Microglia_mouse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Neural_Progenitors"] <- "Neural progenitors"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Neuroblasts"] <- "Neuroblasts"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Oligodendrocyte_Precursor"] <- "Oligodendrocyte precursor"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Oligodendrocytes"] <- "Oligodendrocytes_mouse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Oxytocin_and_Vasopressin_Exp..."] <- "Oxytocin/vasopressin neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "pyramidal_CA1"] <- "Pyramidal (CA1)"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "pyramidal_SS"] <- "Pyramidal (S1)"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Radial_glia_like_cells"] <- "Radial glia-like cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Serotonergic_Neuron"] <- "Serotonergic Neuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Striatal_Interneuron"] <- "Striatal interneuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "Vascular_Leptomeningeal_Cells"] <- "Vascular leptomeningeal"


scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxastrocytes_ependymal"] <- "PI x Astrocytes/ependymal"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxDopaminergic_Adult"] <- "PI x Dopaminergic adult"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxDopaminergic_Neuroblast"] <- "PI x Dopaminergic neuroblast"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxEmbryonic_Dopaminergic_Ne..."] <- "PI x Embryonic dopaminergic neuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxEmbryonic_GABAergic_Neuron"] <- "PI x Embryonic GABAergic neuron" 
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxEmbryonic_midbrain_nucleus"] <- "PI x Embryonic midbrain nucleus neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxendothelial_mural"] <- "PI x Endothelial-mural"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxHypothalamic_Dopaminergic"] <- "PI x Hypothalamic dopaminergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxHypothalamic_GABAergic"] <- "PI x Hypothalamic GABAergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxHypothalamic_Glutamatergic"] <- "PI x Hypothalamic glutamatergic neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxinterneurons"] <- "PI x Interneurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxMedium_Spiny_Neuron"] <- "PI x Medium spiny neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxmicroglia"] <- "PI x Microglia_mouse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxNeural_Progenitors"] <- "PI x Neural progenitors"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxNeuroblasts"] <- "PI x Neuroblasts"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxOligodendrocyte_Precursor"] <- "PI x Oligodendrocyte precursor"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes_mouse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxOxytocin_and_Vasopressin"] <- "PI x Oxytocin/vasopressin neurons"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxpyramidal_CA1"] <- "PI x Pyramidal (CA1)"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxpyramidal_SS"] <- "PI x Pyramidal (S1)"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxRadial_glia_like_cells"] <- "PI x Radial glia like cells"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxSerotonergic_Neuron"] <- "PI x Serotonergic Neuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxStriatal_Interneuron"] <- "PI x Striatal interneuron"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "PIxVascular_Leptomeningeal_C..."] <- "PI x Vascular leptomeningeal"

scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0007416"] <- "Synapse assembly"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0050808"] <- "Synapse organization"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:synprocess"] <- "Synaptic process"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099536"] <- "Synaptic signaling"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099537"] <- "Trans-synaptic signaling"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0050804"] <- "Modulation of chemical synaptic transmission"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0007268"] <- "Chemical synaptic transmission"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:metabolism"] <- "Metabolism"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:presynprocess"] <- "Presynaptic process"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099504"] <- "Synaptic vesicle cycle"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:postsynprocess"] <- "Postsynaptic process"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099072"] <- "Regulation of postsynaptic membrane neurotransmitter"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099173"] <- "Postsynapse organization"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0060078"] <- "Regulation of postsynaptic membrane potential"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0016079"] <- "Synaptic vesicle exocytosis"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:translation_postsynapse"] <- "Translation postsynapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:translation_synapse"] <- "Translation synapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:translation_presynapse"] <- "Translation presynapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099572"] <- "Postsynaptic specialization"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0098794"] <- "Postsynapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0045202"] <- "Synapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0014069"] <- "Postsynaptic density"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0098793"] <- "Presynapse"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0098839"] <- "Postsynaptic density membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099061"] <- "Integral component of postsynaptic density membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0042734"] <- "Presynaptic membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099056"] <- "Integral component of presynaptic membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0099055"] <- "Integral component of postsynaptic membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0045211"] <- "Postsynaptic membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0048786"] <- "Presynaptic active zone"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0048787"] <- "Presynaptic active zone membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0008021"] <- "Synaptic vesicle"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "GO:0030672"] <- "Synaptic vesicle membrane"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:postsyn_ribosome"] <- "Postsynaptic ribosome"
scz_common_geneset$FULL_NAME[scz_common_geneset$VARIABLE == "SYNGO:presyn_ribosome"] <- "Presynaptic ribosome"

#URV
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
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "endothelial-mural"] <- "Endothelial-mural"
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
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxendothelial_mural"] <- "PI x Endothelial-mural"
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
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "ODC1"] <- "Oligodendrocytes_human"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "MG"] <- "Microglia_human"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxODC1"] <- "PI x Oligodendrocytes_human"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxMG"] <- "PI x Microglia_human"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "microglia"] <- "Microglia_mouse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "Oligodendrocytes"] <- "Oligodendrocytes_mouse"

scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxmicroglia"] <- "PI x Microglia_mouse"
scz_urv_geneset$FULL_NAME[scz_urv_geneset$geneset == "PIxOligodendrocytes"] <- "PI x Oligodendrocytes_mouse"

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
#change beta colnames to differentiate common and rare variants
names(scz_common_geneset)[4] <- "beta_common"
names(scz_common_geneset)[6] <- "se_common"
names(scz_common_geneset)[7] <- "p_common"

#merge common and rare variants on full names
merged_genesets <- merge(scz_urv_geneset,scz_common_geneset, by="FULL_NAME")

#add columns for updated fdr adjusted p
merged_genesets$fdr_p_urv <- p.adjust(merged_genesets$p, method = "fdr")
merged_genesets$fdr_p_common <- p.adjust(merged_genesets$p_common, method = "fdr")

#add column for weights and other bits
merged_genesets$weights <- 1/(merged_genesets$se*merged_genesets$se_common)
merged_genesets$significant_urv <- ifelse(merged_genesets$fdr_p_urv < 0.05, "yes", "no")
merged_genesets$significant_common <- ifelse(merged_genesets$fdr_p_common < 0.05, "yes", "no")
merged_genesets$label_plot <- ifelse(merged_genesets$significant_common == "yes" | merged_genesets$significant_urv == "yes", merged_genesets$FULL_NAME, "")
merged_genesets$bigger_label <- ifelse(merged_genesets$significant_common == "yes" & merged_genesets$significant_urv == "yes", merged_genesets$FULL_NAME, "")
table(merged_genesets$bigger_label) #no. of overlapping sig enrichment = 21
`%notin%` <- Negate(`%in%`)
merged_genesets$smaller_label <- ifelse(merged_genesets$label_plot %notin% merged_genesets$bigger_label, merged_genesets$FULL_NAME, "")

#function for computing CIs for spearman correlations
spearman_CI <- function(x, y, alpha = 0.05){
  rs <- cor(x, y, method = "spearman", use = "complete.obs")
  n <- sum(complete.cases(x, y))
  sort(tanh(atanh(rs) + c(-1,1)*sqrt((1+rs^2/2)/(n-3))*qnorm(p = alpha/2)))
}

#correlate betas
cor(merged_genesets$beta,merged_genesets$beta_common, method = "spearman") 
weightedCorr(merged_genesets$beta, merged_genesets$beta_common, weights=merged_genesets$weights, method = "Spearman") 

cor.test(merged_genesets$beta,merged_genesets$beta_common, method = "spearman", exact=FALSE)
spearman_CI(merged_genesets$beta,merged_genesets$beta_common, alpha = 0.05)

#scatter plot
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(ggpubr)

theme_set(theme_light())

#tiff('scz_urv_vs_common_plot_both.tiff',width=7500,height=5500,res=600)
jpeg('scz_urv_vs_common_plot_both.jpeg',width=10000,height=9000,res=950)
ggplot(merged_genesets, aes(x = beta, y = beta_common, size = weights)) + #make point size proportional to weights
 geom_point(aes(color = category), alpha = 0.5) +#colour points by geneset type (category)
  xlab("Beta (ultra-rare variants)") +
  theme(axis.title = element_text(size=10, face = "plain"))+
  #  theme(axis.text.x = element_text(face="plain",size=9.5, angle = 45)) +
  # scale_x_continuous(breaks=merged_genesets$beta,
  #                 labels=merged_genesets$FULL_NAME)+  #put axis ticks at each plotted value 
  guides(size="none") + #remove weights legend
  ylab("Beta (common variants)") +
  theme(axis.title = element_text(size=10, face = "plain")) +
  # theme(axis.text.y = element_text(face="plain",size=9.5, angle = 45))+
  #scale_y_continuous(breaks=merged_genesets$beta_common,
  #                labels=merged_genesets$FULL_NAME)
  # geom_text_repel(data = merged_genesets, aes(label=smaller_label,color = category), show.legend = F, size = 3) + #label points on plot and colour by geneset type
  geom_text_repel(data = merged_genesets, aes(label=bigger_label,color = category, size = 300), show.legend = F, max.overlaps = 15) + #i added size on 18/01/22. if you remove it, the labels are proportional to the size of the cirles
  labs(color = "Geneset") + #rename genesets legend
  scale_color_manual(values = c("#00AFBB", "#E7B800", "wheat4", "#009E73","#0072B2", "#D55E00", "#CC79A7"))
dev.off()