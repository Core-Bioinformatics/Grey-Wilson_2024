# Grey-Wilson_2024

## Cellular heterogeneity allows self-organogenesis of human liver bud in vitro
C. Grey-Wilson, M. Larraz, C. Morell, Á. Blanes-Rodríguez, J. Martínez-García, I. Talon, S. Ghimire, A. Osnato, A. Munteanu, A. Shahsavari, M. Vila Gonzalez, C. Gribben, F. Roos, P. Baptista, I.
Mohorianu & L. Vallier

Uncovering the mechanisms directing the first step of organ formation is essential to understand developmental disorders or regenerative processes. However, early organogenesis in humans remains challenging to study due to the difficulty to access primary tissues and the lack of in vitro models, especially in the context of the hepatobiliary system.  Here, we address this limitation by identifying culture conditions allowing the production of a heterogeneous population of pancreatic, hepatic and biliary progenitors combined with cardiac and endothelial cells from human pluripotent stem cells, which mimic the complex cellular composition of the developing hepato-biliary system. When grown in 3D, these different cells self-organised into ‘HepBud’ organoids recapitulating functional structures of the hepatobiliary region of the primitive gut tube. Single-nuclei RNA sequencing confirms that the differentiation of these organoids follows a natural path of development while uncovering novel cellular interactions involved in this process. Critically, we show that the tissues generated by our approach are functional since they contain organ-specific foetal stem cells, such as hepatoblasts, which can differentiate into hepatocytes and cholangiocytes following in vitro differentiation or intra-splenic injection into mice. To conclude, the HepBud system demonstrates that early progenitors display the capacity to self-organise into complex organ buds with high developmental potential, thereby providing an accessible model to study mechanisms currently impossible to apprehend in humans.

- `Alignment` - CellRanger script to align fastq files

- `QC_Preprocess` - read count matrices into Seurat, and preprocess gene expression data. Add cell type annotations, and apply Clustassess . Expected run order is:
  * 1_QC_norm.R
  * 2_subsetPC_ClustAssess.R
  * DoubletDetection.py
  * 3_Annotate2D.R
  * 4_Annotate3D.R

- `Pseudotime` - Apply custom scripts to generate pseudotime results. Expected run order is:
  * 1_MakeDatasets.R
  * 2_MakeApps.R

- `inVivoComparisons` - Run comparison against in vivo datasets. Expected run order is:
  * CompareToInVivo.R
  * MultiOmicsIntegration-FinalMethod.R

- `GO` - Find and plot GO terms:
  * GO_terms.R

