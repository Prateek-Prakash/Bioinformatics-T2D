# Bioinformatics-T2D
CSC 7300 - Final Project

## Introduction

### Project Motivation
According to the Center for Disease Control, as of 2015, just about 30.3 million people or 9.4% of the United States population has diabetes [1]. This includes children, women, and men of all ages. According to the World Health organization, the number of people with diabetes world wide has risen from 108 million in 1980 to 422 million in 2014 [4]. That’s a staggering increase and spread within the span of 34 years. Diabetes is a major cause of blindness, kidney flare, heart attacks, stroke, and lower limb amputation. It is vital that we gain a better understanding of it as well as possible medicine to help control it. Aside from the widespread importance of diabetes however, I also have a personal motivation for choosing this topic. My dad was diagnosed with type 2 diabetes (T2D) back in 2000, thus this topic is just as much personal. T2D also comprises the majority of the people with diabetes around the world. Both of these reasons are why I have chosen to specifically focus specifically on T2D for this project.  

### Paper Summary
[Paper Link](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2903480/)  

**What is the goal of the study?**  
The goal of the study was to provide insight into abnormal insulin secretion and turnover associated with pancreatic beta-cells from T2D. In a general sense, the study was conducted to enhance the understanding of the pathogenesis of T2D.

**What type of data the authors analyzed?**  
The authors used microarray data for the study. The platform they used was Affymetrix Human X3P Array.

**What analysis did they performed?**  
The authors performed microarray data analysis on the dataset. This involved assessing differentially expressed genes, performing hierarchical clustering, and conducting principal component analysis (PCA) to assess the variation in the expression of genes among the samples. Changes in the expression of functionally related genes were then detected by using Gene Set Enrichment Analysis (GSEA).

**What were the conclusions of the study?**  
At the end of the study, the authors were able to identify many novel changes in gene expression that enhanced their understanding of the pathogenesis of T2D. They were also able to identify a number of gene loci associated with T2D.

## Methods (Overview)
1. Load Dataset
2. Log2 Transform
3. Differential Expression Analysis (Empirical Bayes Moderated T-Statistics Test)
4. Hierarchical Clustering (RAW Data)
5. Hierarchical Clustering (Differentially Expressed Data With CI.L Cutoff)
6. Hierarchical Clustering (Differentially Expressed Data With P.Value Cutoff)
7. Volcano Plot
8. Pathway Analysis

## Methods (Detailed)
To begin, I first ran the selected dataset through GEO2R. I grouped the samples into diabetic and control groups. After this, I started my code with the basic R script that was generated through GEO2R. The script by itself goes through and does differential expression analysis on the dataset after normalizing the expression values. The test that was used for differential expression analysis was an empirical Bayes moderated t-statistics test. After this, the normalized expression data was plotted using a boxplot.

After this, I annotated the resulting data with gene symbols. During this step, I also viewed the resulting data in two ways. The first way was by ordering the data by the left confidence interval (CI.L). This way I could compare the genes of interest based off the cuttoff giving in the paper. The other way was by ordering the data by the P value (P.Value) for similar reasons.

Next I wanted to try doing some clustering to see if I could get a definite split between the two groups as the authors of the paper did. First, I tried the clustering three different times. The first time I tried with the RAW data. For the other two times, I used a subset of the differentially expressed data. One subset was based off the 1.2 CI.L cutoff and another subset was based off the 0.05 P.Value cutoff.

Finally, I created a volcano plot to display possible genes of interested. Then using GeneMania, I took the genes of interest from the plot and ran it through their application to get significant pathways [3]. I took this same list and also used it as input in Enrichr, to see if I could get any significant pathways [9].

## Data Description
[Dataset Link](https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS3782)  

* Title: Type 2 Diabetes: Pancreatic Beta-Cells
* Organism: Homo Sapiens
* Platform: Affymetrix Human X3P Array
* Count: 20 Samples (10 T2D, 10 Non-Diabetic Control)

## Results
### Figure 1: Normalized Expression Plot
![Normalized Expression](/Results/Exp-Norm.png?raw=true "Normalized Expression")  
### Figure 2: Hierarchical Clustering (RAW Data)
![Hierarchical Clustering (RAW Data)](/Results/HC-RAW.png?raw=true "Hierarchical Clustering (RAW Data)")  
### Figure 3: Hierarchical Clustering (Diff. Exp. Data :: CI.L)
![Hierarchical Clustering (Diff. Exp. Data :: CI.L)](/Results/HC-Diff-Exp-CI.png?raw=true "Hierarchical Clustering (Diff. Exp. Data :: CI.L)")  
### Figure 4: Hierarchical Clustering (Diff. Exp. Data :: P.Value)
![Hierarchical Clustering (Diff. Exp. Data :: P.Value)](/Results/HC-Diff-Exp-P.png?raw=true "Hierarchical Clustering (Diff. Exp. Data :: P.Value)")  
### Figure 5: Volcano Plot
![Volcano Plot](/Results/Volcano-Plot.png?raw=true "Volcano Plot")  
### Figure 6: GeneMania Network
![GeneMania Network](/Results/GeneMania-Network.jpg?raw=true "GeneMania Network")  

## Disucssion
Initially, the results looked good from the differential expression analysis. However, upon comparing my results to that of the authors of the study, they did not match at all. As mentioned before, gene expressions which were considered important were selected based off CI.L and P.Value. Using the same 1.2 cutoff for CI.L and 0.05 cutoff for P.Value as the study, I onyl got three genes that actually met those requirements. Table 2 of the study contains the genes of interest and none of the three genes I found were at the top of the list. Also, the values they have in that table don't match up with my CI.L nor my P.Value.

With that in mind, I decided to contain forward to see if I could get some good results with clustering even though my results looked different than the study. Unfortunately, it seems as though that also wouldn't turn out right. As I mentioned previously, I used three different sets in the clustering to see what results I would get. The RAW data resulted in possibly the worst clustering of the three. With using a subset of the processed and differentially expressed data, I got some better grouping however I still couldn't get two distinct groups of T2D and control.

Now, I'd like to discuss why my results might have been completely off compared to the study. Once again this is all speculation as I don't know if these are exactly the reasons. One possiblity is that the normalization methods varied. In the study, the data was normalized using invariant set normalization. This was done with dChip, however I don't know if the data from GEO was post-normalization [2]. Along with this, I don't know if they transformed their data either but that was the normalization I had to resort to. Another possibility is that their may have been differences in the t-statistics test as well. There was no mention of the moderated t-statistics test in the study however, using GEO2R results in an empirical Bayes moderated t-statistics test. Finally, the differential expression analysis in the study was done using DNA-Chip Analyzer software. This could be a huge reason why we may have gotten different results. The methods the DNA-Chip software utilizes could be different than the methods I used. Overall, any or all of these reasons could be why my results were completely different.

Finally, there was the results of the pathway analysis methods that I chose to do. Unfortunately, the GeneMania search didn't result in any good results. SERPING1 and F12 were the only genes that go highlighted, however neither of them seem to be relevant when it comes to T2D. This left me with Enrichr. Using Enrichr, I did get one pathway which seemed to be significant and relevant to T2D. HSA04920 Adipocytokine Signaling Pathway seems to be directly related to the deleveopment of T2D amongst other insulin related issues [8]. The pathway seems to be related to an insulin signaling pathway as well and can be seen on [KEGG](http://www.genome.jp/kegg-bin/show_pathway?hsa04920) [5, 6, 7].

## References
1. Center for Disease Control. (2017, July 17). National Diabetes Statistics Report, 2017 [Online]. Available: https://www.cdc.gov/diabetes/pdfs/data/statistics/national-diabetes-statistics-report.pdf
2. Edgar R, Domrachev M, Lash AE.  
Gene Expression Omnibus: NCBI Gene Expression And Hybridization Array Data Repository.  
Nucleic Acids Res. 2002 Jan 1;30(1):207-10
3. Warde-Farley D, Donaldson SL, Comes O, Zuberi K, Badrawi R, Chao P, Franz M, Grouios C, Kazi F, Lopes CT, Maitland A, Mostafavi S, Montojo J, Shao Q, Wright G, Bader GD, Morris Q.  
The GeneMANIA Prediction Server: Biological Network Integration For Gene Prioritization And Predicting Gene Function.  
Nucleic Acids Res. 2010 Jul 1;38 Suppl:W214-20
4. World Health Organization. (2016, Nov.). Diabetes. [Online]. Available: http://www.who.int/mediacentre/factsheets/fs312/en/
5. Kanehisa, Furumichi, M., Tanabe, M., Sato, Y., Morishima, K.  
KEGG: New Perspectives On Genomes, Pathways, Diseases And Drugs.  
Nucleic Acids Res. 45, D353-D361 (2017).
6. Kanehisa, M., Sato, Y., Kawashima, M., Furumichi, M., and Tanabe, M.  
KEGG As A Reference Resource For Gene And Protein Annotation.  
Nucleic Acids Res. 44, D457-D462 (2016).
7. Kanehisa, M. and Goto, S.  
KEGG: Kyoto Encyclopedia Of Genes And Genomes.  
Nucleic Acids Res. 28, 27-30 (2000).
8. Eduardo Esteve, Wifredo Ricart, José Manuel Fernández-Real.  
Adipocytokines And Insulin Resistance: The Possible Role Of Lipocalin-2, Retinol Binding Protein-4, And Adiponectin.  
Diabetes Care. 2009 Nov; 32(Suppl 2): S362–S367.
9. Edward Y Chen, Christopher M Tan, Yan Kou, Qiaonan Duan, Zichen Wang, Gabriela Vaz Meirelles, Neil R Clark And Avi Ma’ayan.  
Enrichr: Interactive And Collaborative HTML5 Gene List Enrichment Analysis Tool.  
