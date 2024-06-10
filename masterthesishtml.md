```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r, echo=FALSE}
knitr::include_graphics("COVER TESI PNG.png")
```
\newpage

\tableofcontents

\newpage

\thispagestyle{empty}

\vspace*{\fill}
\begin{center}
    \textit{To peaceful intentions only.}
\end{center}
\vspace*{\fill}

\newpage

# Abstract

After its long reputation as an extremely polluted water river and anthropogenic mining waste dump, studies support Rio Tinto (southern Spain, Iberian Pyrite Belt) to be an extremely acidic environment where life relentlessly thrived long before human history. It has become clear that in this extreme environment, there are strong relationships between living and nonliving components at the microscale, resulting in the formation of micro-niche-based ecosystems. The site has therefore become a terrestrial analog of first interest for astrobiological and planetary science studies in terms of the search for life on Mars. The acidic environment is the product of the chemolithotrophic activity of microorganisms aggressively targeting sulfides (pyrite, chalcopyrite), here abundant, causing the leaching of iron and sulfur. This contributes to the formation of a variety of minerals, mainly gypsum, jarosite, goethite, and hematite, all of which have been detected on the Red Planet.

Identifying and discretizing sulfides and iron-bearing sulfates from orbit and landed missions has been a relevant method for searching for life on Mars, notably distinguished by its iron-sulfur-rich composition. Similar mapping of sulfide and sulfate distributions on easy-to-access terrestrial analog is critical to improving our ability to interpret data from other worlds and contextualize astrobiolgical observations.

In this work, we present the spectroscopic analysis of remote sensing data over Rio Tinto, focusing on mapping the distribution of sulfides and sulfates as a proxy for the presence of biosignatures. We have studied multi- and hyper-spectral data from orbital and airborne spectrometers, cross-checking evidence from different datasets. The results have been cartographically formatted to support a geologic fieldwork campaign held at the Rio Tinto in November 2023.

\newpage

# 1 INTRODUCTION

Rio Tinto (figure 1) is a 100 km acidic river that streams in the southern Spanish region of Andalusia digging its way through the Iberian Pyrite Belt (IPB), eventually heading to the mediterranean sea. The northern part of the basin is characterized by a large mining area, whose ores have been exploited since the Roman colonization. The intense anthropic activity, contamination and acid mine drainage had long kept the biological interest of science away. Recent studies, since early 2000s, demonstrated it to be a biologically active extreme environment, whose chemolithotrophic ecosystem is now an established astrobiological analog for Mars \cite{roach2006}.  These extreme conditions found in the river are the direct consequence of the active metabolism of chemolithotrophic microorganisms thriving in the polymetallic sulphides present in the IPB in high concentration \cite{gomez2011}. Bacteria from all major taxonomic groups have been found in metal-rich habitats such as mineral ores, acidic soils, or polluted environments. Chemolithotrophic Fe and S oxidizers such as Thiobacillus ferrooxidans (now named Acidithiobacillus ferrooxidans) and Leptospirillum ferrooxidans solubilize Fe and S contained in metal sulphides, with the concomitant release of associated metals \cite{glasauer2013}. The bio-oxidation of sulphides such as pyrite (FeS$_2$), produces a sulphate-enriched acidic solution (pH between 0.8 and 3) that prevents ferric iron precipitation, contrary to what occurs in neutral conditions \cite{fernandez2004}. This is also the reason why the waters of the river appear of a red colour, which gives the name to the whole stream. The consequently high Fe$^{3+}$ and sulphates concentration is responsible for the formation of a variety of minerals, mainly gypsum, jarosite, goethite, and hematite, all of which have been detected on different regions of Mars, and some of these, such as the case of the sulphate jarosite, even facilitate life establishment, as they provide micro-niches against UV radiation \cite{gomez2011}.

```{r, label = "fig:1", echo=FALSE, fig.cap="Rio Tinto streams in the northern area of Ravine of Fools. Photo by Giacomo Panza."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/RAVINE OF FOOLS.jpg")
```

In the studies of \cite{fernandez2004}, the river and its environment provide one of the best-known terrestrial analogues to Noachian deposits on Meridiani plain on Mars. However, it is to be considered, as written by \cite{amils2014}, that the high sulphates and ferric-sulphates concentration in the rocks on Meridiani outcrops indicates diagenesis in environments with extremely limited amounts of water. The words from \cite{roach2006} consider Rio Tinto less a specific analog to Mars and more a model system in which minerals and their distribution can be studied on a variety of scales (fig. 2). Underlaying these interpretations, the common point is that similar lithological units were observed between our planet and Mars. Deepening the knowledge about their formation and distribution is undoubtedly interesting and necessary to shed the light about their origin, whether partly biogenetic or not. At Rio Tinto, the thread along the sulphide ores long mineralized in the IPB, their bio-oxidation thanks to prokaryotic species, and the iron-sulphates subsequently produced feeds the astrobiological potential of this scenario and the approach of our study. Indeed, our interest is devoted to the identification of the two extreme mineral components of this bio-geo-chain mechanism. This will help us better describe and further understand how these patterns have developed and can be identified on other terrestrial areas as well as on other rocky planets. Different sources of remote sensing data were used to cross check results at different scales. The instruments with lower resolution will be exploited to draw an overview map using band ratios, whereas hyperspectral data will give us detailed spectra of the mineral composition for a discrete classification.

```{r, label = "fig:2", echo=FALSE, fig.cap="A view over Pequeño Grand Canyon. Here it is possible to distinguish different mineralogic compositions. Reddish layers are oxidized rocks, mostly chunks of oxidized slates. In the lowest part, clays accumulate as grey weathered sands. Photo by Giacomo Panza."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PGZ LANDSCAPE.jpg")
```

## 1.1 Geologic Setting
The Iberian Pyrite Belt (IPB) extends from Portugal to Spain (fig. 3B), covering an area of around 250 km by up to 60 km \cite{martin2015}. It lays in the South Portuguese Zone (SPZ), which it is part of the Iberian Massive (IM) (fig. 3A). The IM was formed in Devonian and Carboniferous periods due to collision of Laurentia with Gondwana. The IPB contains volcanic massive sulphide (VMS) and stockwork deposits, which are strata-bound accumulations of sulphide minerals that precipitated at or near the seafloor in spatial, temporal, and genetic association with contemporaneous volcanism \cite{franklin2005}. According to \cite{ross2014}, this area holds the 23\% of the global VMS tonnage and can be associated to a volcanism of felsic domes and lavas. The succession here presents, from ancient to recent, the Phyllite-Quartzite Group (PQ), the Volcano Sedimentary Complex (VSC) and the Flysch Group (F), as mentioned by \cite{demello2022}. The first succession is dominated by shales and quartzite of late Devonian, meant to be deposited in a shallow continental platform. The overlaying VSC formed between Famennian and Visean, and it presents volcanic rocks interbedded with shales, sedimentary rocks, and massive sulphides. The Flysch Group is a turbidite-like sequence of siliceous and carbonaceous shales and conglomerates of late Visean age.

```{r, label = "fig:3", echo=FALSE, fig.cap="A) Simplified geological map of the Iberian Massif showing the location of the South-Portuguese Zone. B) Simplified geological map of the South-Portuguese Zone showing the location of the Rio Tinto–Nerva Unit (black box) \\cite{valenzuela2011}."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/IPB A. Valenzuela.png")
```

\newpage

The Rio Tinto district, in the Nerva Unit, is part of the IPB and, according to the models proposed by \cite{martin2015}, it originally formed as a rollover anticline that behaved like a graben, slipped down between a Southern Fault and a Northern Fault during a first transtensional period and was then gradually squeezed upwards during a later traspressional period (Fig. 4). Along the two major faults hydrothermal flow was facilitated, and there it is where we can find the massive sulphides mineralization and the Copper and Pyrite stockworks.

```{r, label = "fig:4", echo=FALSE, fig.cap="The tectonic momvements that contributed to the actual configuration of Rio Tinto district, as reported by \\cite{martin2015}"}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/Martin Izard 3D.png")
```

\newpage

## 1.2 Remote Sensing
Recalling what \cite{goetz1981} conceived since the dawn of this science, remote sensing refers to all the arts and techniques of measurement and interpretation of phaenomena from afar. A powerful subset of this scientific approach is imaging spectrometry, explained by the same author in 2009 as the acquisition of images in hundreds of contiguous, registered spectral bands such that for each pixel a radiance spectrum can be derived. This contiguous image sampling makes it possible to compare pixel spectra with spectral databases. Materials can thus be described and identified by their spectral signature, which will present diagnostic reflectance responses at certain wavelengths. The reflectance spectrum in the region from 0.4 to 2.5 µm can be used to identify a large range of surface cover materials, as well as minerals \cite{goetz1985}.  For this research, three instruments were used, as presented in the materials and methods section. The core difference is their spectral resolution, or rather the ability to sample multiple bands. Briefly (fig. 5), hyperspectral technologies allow a much larger number of samples for far more detailed spectra compared to multispectral systems. The advantages of experimenting with these instruments is the amount of data that they can produce for Earth and Planetary studies, and the more affordability they have when compared to landing missions, rovers and transport spacecrafts. The limits are that without a proper ground truthing or tangible cross-check, only modelling is possible.

```{r, label = "fig:5", echo=FALSE, fig.cap="A brief walkthrough that shows how multi and hyperspectral tools work. The only difference is the spectral resolving power, much higher for hyperspectral tools than multi."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/REMOTE SENSING SCHEME.png")
```

\newpage

# 2 MATERIALS AND METHODS

## 2.1 Data Source
As mentioned, the raw datasets come from three different sensors mounted on satellites that orbit around Earth or built for airborne. Everything was downloaded with the standard procedures indicated by the official provider. A brief explanation of the machinery is written below.

### 2.1.1 Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER)
The Advanced Spaceborne Thermal Emission and Reflection Radiometer (ASTER, fig. 6) is an advanced multispectral imager that was launched on board NASA’s Terra spacecraft in December 1999. ASTER covers a wide spectral region with 14 bands from the visible to the thermal infrared with high spatial, spectral, and radiometric resolution \cite{abrams2002}. Spatial resolution is 15 m for the VNIR, 30 m for the SWIR and 90 m for the TIR wavelengths. It counts 14 bands, 3 for the VNIR, 6 for the SWIR and 5 for the TIR). The AST\textunderscore L1T granule is a multi-file product, which includes an HDF-EOS2 science data product file, full-resolution images, and a metadata file. The AST\textunderscore L1T product is created by performing the geometric and radiometric corrections on the original AST\textunderscore L1A image data. The result is projected onto a rotated map (rotated from “path oriented” coordinate to UTM grid north-up) at full-instrument resolutions \cite{asterproductspecification}. Our dataset consists of one L1T product acquired in December 2007.

```{r, label = "fig:6", echo=FALSE, fig.cap="Artistic impression of Terra Spacecraft carrying ASTER machinery. Image courtesy of NASA."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/ASTER SPACECRAFT.png")
```

### 2.1.2 Airborne Visible InfraRed Imaging Spectrometer - Next Generation (AVIRIS-NG)
AVIRIS-NG, a new airborne hyperspectral sensor based on the previous AVIRIS model, has 425 spectral bands with high spatial (ca.~4 m) and spectral resolution (5 nm) in 0.4\textmu m–2.5\textmu m spectral range \cite{tripathi2019}. Signal to noise ratio of this instrument is about 300 times higher. It can deliver L2 data products, which are atmospherically corrected surface reflectance data products. The atmospheric correction is based on the Atmosphere Removal Algorithm (ATREM) program developed by \cite{gao2009}. The dataset used for the purpose of this research consists of a georeferenced mosaic based on two different images dated July 2021.

### 2.1.3 Hyperspectral Precursor of the Application Mission (PRISMA)
PRISMA (fig. 7) is a medium-resolution hyperspectral imaging satellite, developed, owned, and operated by ASI (Agenzia Spaziale Italiana). The Hyperspectral Camera (HYC) sensor is a prism spectrometer for two bands, VIS/NIR (Visible/Near Infrared) and NIR/SWIR (Near Infrared/Shortwave Infrared), with a total of 237 channels across both bands. The HYC module has a spatial resolution of 30 m and a swath width of 30 km. Our dataset consists of a L2D product, which are atmospherically corrected and geocoded, dated June 2021.

```{r, label = "fig:7", echo=FALSE, fig.cap="Artistic impression of PRISMA. Image courtesy of ASI (Agenzia Spaziale Italiana)."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMA.jpg")
```

\newpage

## 2.2 Data Preparation
Different methods were used for multispectral and hyperspectral data. The software we have used to perform calculations on imagery was RStudio (based on R open-source programming language) for some basic initial analysis and IDL ENVI® 5.3 for the effective processing and data visualization. At the bottom of this paragraph (Glossary section), descriptions of the different computing techniques and functions (marked with superscripts) can be found.

### 2.2.1 Multispectral data
From ASTER dataset VNIR and SWIR wavelengths only were used. Data was pre-processed in ENVI® \cite{envitutorials2000} by performing a radiometric calibration, followed by a layer stacking and an interleave conversion to BIL data format to finally proceed with FLAASH atmospheric correction (Fast Line-of-sight Atmospheric Analysis of Spectral Hypercubes) and retrieve reflectance from radiance. Lastly, the following IDL equation (eq. 1) was applied to scale the range of reflectance values between 0 and 1:

$$
(B1 le 0)*0+(B1 ge 10000)*1+(B1 gt 0 and B1 lt 10000)*float (b1)/10000
$$
New data was finally produced thanks to the RGB projection of significant band-ratios (from \cite{asterindexmanual2004}) that locate coarse mineral distribution. A further classification was executed by defining regions of interest (ROI), extracting endmembers as the mean spectrum for each region and mapping with SAM$^a$ tool. In figure 8 it is possible to visualize ASTER data processing scheme. A video guide is available at \href{https://drive.google.com/drive/folders/1EXbIKkHOI9dLKn1H_c6MsLGVfMZDPBQD?usp=drive_link}{ASTER\textunderscore data\textunderscore preparation}.

```{r, label = "fig:8", echo=FALSE, fig.cap="A workflow scheme for ASTER dataset."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/ASTER SCHEME.png")
```

\newpage

### 2.2.2 Hyperspectral data
The PRISMA dataset is stored in he5 format (an extension belonging to hierarchical data format 5, HDF5). This encryption makes it difficult to extract information with ENVI version 5.3, and a transformation via R programming language was executed, thanks to the \textit{prismaread} package \cite{prismaread2020}. With this code (shown in figure 9) files were converted to ENVI format and further processed.

```{r, label = "fig:9", echo=FALSE, fig.cap="RStudio view of the pr_convert function from \\cite{prismaread2020}."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMAREAD.png")
```

Fields compiled as shown in the figure above allowed us to obtain an adequate output for hyperspectral processing. PAN image was generated to have a high-resolution image of the area. VNIR and SWIR layers were generated separately because otherwise there would have been errors on the file. These two channels were later joined in ENVI as a single layer-stack. Beware, VNIR and SWIR numeration is nonlinear: 1 - 66 for VNIR and 1 - 173 for SWIR, for a total of 239 bands (whose 5 with null value as they are used for overlapping the spectra). With that in mind, in figure 10 the removed bands are reported, whereas in the following figure 11 a spectrum extracted from a random pixel is reported, showing bad bands removal gaps. PRISMA dataset was used to obtain a map via SAM algorithm but giving as input endmembers a list of library spectra extracted from USGS spectral library and reported in figure 12. SAM threshold angle was 0.200 for all endmembers. In figure 13 it is possible to visualize PRISMA data processing scheme.

```{r, label = "fig:10", echo=FALSE, fig.cap="A list of removed bands for the two PRISMA channels, VNIR and SWIR."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMA BR.png")
```







```{r, label = "fig:11", echo=FALSE, fig.cap="Single pixel spectrum for PRISMA dataset. Here it is possible to see gaps belonging to removed bands of table 1."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMA BR GRAPH.png")
```






```{r, label = "fig:12", echo=FALSE, fig.cap="Name of different library spectra. These were chosen and grouped by mineral class.}"}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/ENDMEMBERS PRISMA.png")
```






```{r, label = "fig:13", echo=FALSE, fig.cap="A workflow scheme for PRISMA dataset."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMA SCHEME.png")
```

\newpage

For AVIRIS-NG datasets, no pre-process was needed. The files already contained reflectance data, which means it was already atmospherically corrected, so, after mosaicking and sub setting, bad bands were discarded after a first overlook of the file. In figure 14 and figure 15 it is possible to see which bands were removed and how the spectrum looked like. Data production, instead, aimed to the identification of a pixel based mineral distribution. The steps required by the Spectral Hourglass procedure were respectively a MNF$^b$ followed by a PPI$^c$, which helped highlight significant data. Then, NDV$^d$ tool helps visualize how this data is distributed, so the user can identify clusters (classes) and extract a mean spectrum for each, thus endmembers can be defined. 


* MNF takes all the n good bands as input and returns n MNF bands as output. Of the first 20 data-rich output MNF bands, bands number 12 and 14 were discarded as they resulted corrupted. A total of 18 bands were accepted.
* PPI takes the 18 MNF bands as input. 10000 iterations were planned, 250 per cell with a 2.50 threshold. A band threshold to ROI filtered only pixels with value between 10 and 180.
* NDV was executed on the 18 MNF bands for the selected PPI output pixels and classes were found. For each class the mean true spectrum and MNF spectrum was extracted.

\newpage

The endmembers are named after their matching spectrum according to a standard spectral library comparison (USGS spectral libraries and JPL spectral libraries were resampled and compared). The last step is to map a distribution of the endmembers on the original image using the SAM, MTMF$^e$ and LSU$^f$ techniques. Outputs were implemented into QGIS to present the data. In the figure 16 it is possible to visualize AVIRIS-NG data processing scheme.

```{r, label = "fig:14", echo=FALSE, fig.cap="A list of removed bands for AVIRIS-NG dataset."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/AVIRIS BR.png")
```


```{r, label = "fig:15", echo=FALSE, fig.cap="Single pixel spectrum for AVIRIS-NG dataset. Here it is possible to see gaps belonging to removed bands of table 3."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/AVIRIS BR GRAPH.png")
```



```{r, label = "fig:16", echo=FALSE, fig.cap="A workflow scheme for AVIRIS-NG dataset."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/AVIRIS SCHEME.png")
```

\newpage


## 2.3 Glossary


### *2.3.1 Spectral Angle Mapper (SAM)$^a$*
SAM is an automated method for comparing image spectra to individual spectra that determines the similarity between two spectra by calculating the spectral angle between them, treating them as vectors in a space with dimensionality equal to the number of bands \cite{envitutorials2000}.

### *2.3.2 Minimum Noise Fraction (MNF)$^b$*
MNF Transform is a method like Principal Components used to segregate noise in the data, determine inherent data dimensionality, and reduce computational requirements for subsequent processing \cite{green1988,Boardmankruse1994}. For hyperspectral data, MNF divides data space into two parts: one with large eigenvalues and coherent eigen images and the second with near-unity eigenvalues and noise-dominated images. It puts most of the interesting information into just a few spectral bands \cite{envitutorials2000}.

### *2.3.3 Pixel Purity Index (PPI)$^c$*
The PPI function finds the most spectrally pure or “extreme” pixels in multispectral and hyperspectral data \cite{Boardmankruse1994}. These correspond to the materials with spectra that combine linearly to produce all the spectra in the image. The output is an image (the PPI image) in which the digital number (DN) of each pixel in the image corresponds to the number of times that pixel was recorded as extreme. Thus, bright pixels in the image show the spatial location of spectral endmembers \cite{envitutorials2000}.

### *2.3.4 N-Dimensional Visualizer (NDV)$^d$*
ENVI’s N-Dimensional Visualizer is an interactive n-dimensional scatter plotting paradigm that allows real-time rotation of scatterplots from n-d dimensions \cite{Boardman1995}, where n is the number of bands. The coordinates of the points in n-space consist of “n” values that are simply the spectral radiance or reflectance values in each band for a given pixel \cite{envitutorials2000}.

### *2.3.5 Mixture Tuned Matched Filtering (MTMF)$^e$*
Mixture Tuned Matched Filtering (MTMF) The MTMF method \cite{Boardman1998} improves upon the MF method by providing better selectivity of targets. It is useful for detecting and discriminating among multiple rare targets whose spectral signatures are similar to the background. The MTMF method recognizes that targets actually replace some of the background signature in a pixel, not add to them. It uses the same statistical method as MF but incorporates elements of a linear mixing model. It can provide accurate mapping of very small sub-pixel targets with a low number of false positives \cite{boardmankruse2011}. The MTMF method requires an MNF-transformed image for input. It calculates an MF Score along with an infeasibility measure, which describes the likelihood of each pixel being a mixture of the known target and background materials. The infeasibility measure allows analysts to identify and reject false positives \cite{NV5}.

### *2.3.6 Linear Spectral Unmixing (LSU)$^f$*
LSU is used to determine the relative abundance of materials that are depicted in multispectral or hyperspectral imagery based on the materials’ spectral characteristics. The reflectance at each pixel of the image is assumed to be a linear combination of the reflectance of each material (or endmember) present within the pixel. For example, if 25\% of a pixel contains material A, 25\% of the pixel contains material B, and 50\% of the pixel contains material C, the spectrum for that pixel is a weighted average of 0.25 times the spectrum of material A, plus 0.25 times the spectrum of material B, plus 0.5 times the spectrum of material C. So, given the resulting spectrum (the input data) and the endmember spectra, Linear Spectral Unmixing solves for the abundance values of each endmember for every pixel. The number of endmembers must be less than the number of spectral bands, and all the endmembers in the image must be used \cite{envitutorials2000}.

### *2.3.7 Normalized Difference Vegetation Index (NDVI)$^g$*
It is an index used to measure vegetation density and is based on the principle that green plants reflect more light in the near infra-red than they do in the visible red spectrum. Thus, it is calculated by the following (eq. 2) band ratio:

$$
	(VNIR - RED)/(VNIR + RED)= NDVI 
	\label{eq:2}
$$

\newpage

# 3 RESULTS AND DISCUSSION 
The data we obtained reveal some interesting features but must be interpreted carefully. The three instruments were used in order of spectral resolution, so to deepen the analysis. The outputs are raster files showing maps of the designated areas, each with its key for data reading and interpretation. 

## 3.1 ASTER 
The band ratio applied to ASTER bands is shown in figure 17 and it was a starting point to discriminate the most interesting spots on our area of interest. The subset shows the mining zone and some of the surrounding environment. It was chosen to combine NDVI$^g$ band ratio with band ratios that highlight alteration and ferric iron (as reported by \cite{asterindexmanual2004} in ASTER Mineral Index Processing Manual), so to spot major silicate distribution and sulfates and oxides distribution. Highlighting vegetated areas was important to understand where any investigation is worthless. It must be mentioned that the red colour (related to vegetation) also covers areas where shadow is strongly present. This is because the low visible light reflected by a shady zone causes similar feedback to that of dense vegetation for the NDVI$^g$. This can be encountered for example in the waters of Corta Atalaya pit, where the red NDVI-linked colour is generated by higher reflectance values in the band 3 (VNIR) than in the band 2 (RED) just for a matter of amount of reflected light. The water bodies are indicated by the dark and dense bluish colour, which is given by high absorbance on all band ratios combined. Furthermore, the other band ratios, mapped on green and blue, emit significantly less signal in shady areas. With this in mind, we can clearly identify an interesting, exposed area, which is supposed to host silicates and ferric iron bearing minerals. 

```{r, label = "fig:17", echo=FALSE, fig.cap="The colour composite image showing NDVI band-ratio, Alteration band-ratio and Ferric Iron band-ratio."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/BAND RATIO ASTER.png")
```

The figure 19 shows a useful comparison map obtained with SAM technique by training data with ROIs mean spectra. SAM thresholds are reported in figure 18.


```{r, label = "fig:18", echo=FALSE, fig.cap="SAM threshold values for each endmember."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/ENDMEMBERS ASTER.png")
```

\newpage

Some blue areas in figure X reporting ferric iron are now indicated as anthropized, and considering the building materials used in a town, we can now better understand where the limit of band ratioing is. Flooded areas and water basins are also better highlighted, and we can see some anthropogenic material spread in the northwestern ponds as these are artificially excavated and delimited by roads, whose spectrum contributes to defining the anthropized endmember, as its ROI is picked up in the towns, where some are built unpaved or by using gossan material \cite{bedini2022}. Some pixels belonging to Cerro Colorado open pit in figure X were identified as water. The area is almost completely shaded and not accessible to public, so it is not possible to recognise water basins, but it is likely that ponds were present, as it happens, even nowadays (figure 20), in the other open pits (Corta Atalaya and Peña de Hierro).


```{r, label = "fig:19", echo=FALSE, fig.cap="The SAM map reporting endmembers distribution. Data was trained to find endmembers using ROIs."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/ASTER ROIS.png")
```


```{r, label = "fig:20", echo=FALSE, fig.cap="Open pit of Peña de Hierro (norhtern area). Water is filling much of the hole. Photo by Giacomo Panza."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/Pena de Hierro.jpg")
```

\newpage

## 3.2 PRISMA

PRISMA hyperspectral data was used with focus on the bare ground spots previously highlighted by ASTER. Because now we deal with a hyperspectral dataset, SAM was used to map input minerals previously mentioned in figure 14. As shown in figure 21, the clays class here confirms the distribution of what is seen as alteration for ASTER data, as both classes indicate silicates. In fact, clays are widespread and sulfates (most of which Fe$^{3+}$ bearing) lay in smaller and thinner areas. In the north-west, where mine tailings and contaminated ponds are \cite{bedini2022}, it is shown the presence of stripe-accumulated sulfates. In these areas water evaporates, and crystallization can occur. The intense blue colour in the RGB image from ASTER agrees with the presence of ferric iron bearing sulfates such as Jarosite, Coquimbite, Copiapite, but distribution is much different, probably because PRISMA dataset is dated 2022. Though on the northern slope of Corta Atalaya sulfates were mapped (agreeing with the distribution reported by \cite{bedini2022}), no oxides were found in this spot. Oxides, however, were mapped in blue above Cerro Colorado and show relationship with pinkish-to-brown rocks found in this area. Sulfides distribution is often related to the presence of contaminated waters. The two spectra are similar for showing a reflectance peak in the visible light and then a constant spectrum in the remaining wavelengths, so it is still hard to determine the presence of this mineral class.


```{r, label = "fig:21", echo=FALSE, fig.cap="SAM map of the library spectra for the minerals from table 2 chosen as endmembers."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/PRISMA ROIS.png")
```


\newpage

## 3.3 AVIRIS-NG
Our PRISMA hyperspectral classification, combined with the previously analysed ASTER multispectral dataset, has led us to consider the Pequeño Grand Canyon zone (PGZ) as target for this last survey. This was an interesting area as it is free of large anthropized structures in the bare ground dominated surface, which hosts exposed altered material including clays and sulfates. So, according to the processing in figure 16, endmembers were found. In figure 22, main spectra are reported and compared with their USGS library reference spectra. Clays class showed a spectrum like that of montmorillonite for the absorption feature near 1.45 µm region (also reported by \cite{clark1990}), and near 2.2 µm region caused by the Al-OH bond. Also, montmorillonite is reported to be present in the Rio Tinto district, as well as illite, and kaolinite \cite{bedini2022}. It must be said that the lower the resolution, these three minerals present a very similar spectrum. With our instrument it impossible to determine the exact mineral species, so we grouped these possibilities into the so-called Clays class. The same can be said for the Sulfates class spectrum, which is mainly compared to copiapite library spectrum, but also Jarosite and Coquimbite. The absorption feature in the 0.8 to 1 µm region is caused by the presence of $Fe^{3+}$ for Jarosite, as well as the tight absorption around 1.45 µm region caused by OH \cite{clark1990}. For the Oxides class, hematite, and goethite spectra from USGS library were the best matches, with hematite being the first. Our spectrum has much more noise accumulated in the 1.5 to 2.5 µm region, and we assumed it to be flat as the reference spectrum. The $Fe^{3+}$ absorbance feature around 0.9 µm is relevant, as well as the similarity of the whole region from 0.5 to 1.0 µm for both spectra. The Sulfates and oxides class is interpreted to represent both mineral classes, as the spectrum was showing features for both classes. This is shown by the absorption feature around 1.45 µm matching with the jarosite spectrum, and the 0.5 to 1.0 µm region agreeing with both oxides and sulfates, with peaks matching slightly better the ones of the goethite spectrum. Goethite was indeed the first match according to the weight methods used (described in materials and methods section), but jarosite was anyway considered.


```{r, label = "fig:22", echo=FALSE, fig.cap="Comparison of endmembers' spectra and library spectra. In red the mean spectrum for the given pixel classes. In black or yellow the different library spectra. See how well some features overlay, or how shifted they are sometimes."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/AVIRIS SPECTRA.png")
```

Firstly (figure 23), a ROI based classification was created to identify a general landcover type. As a result of this, we want to focus our mineralogic interest on the non-classified pixels. Next to this, the SAM1 algorithm describes the mineralogy of the area. The oxides class has highlighted an area with pinkish to rocky brown cover like that reported by PRISMA in the previous analysis. Sulfates and oxides are found as main cover in the high Ravine of Fools zone (top left corner). This area is covered by reddish oxidated rocks and can host occasional water flows. Clays cover most of the bare ground in the PGZ. These slopes appear greyish to whitish, made by loose sediment and rocks in the higher spots and by black and white laminated hills of sandy to silty sediment in the lower grounds. Sulfates are found in the lowest part of the PGZ, where water was reported flowing in November 2023 by our team. The harsh climate from spring to summer brings in a lot of evaporation events, both at large and small scale. Biologically speaking, sulphury crystallized popcorn like structures of bacteria \cite{gomez2011} were reported on wet rocky surface both in open air locations and niche-like environments such as small caves or concavities in the host rock, where evaporation clearly occurs, and the air is humid. This bacterial efflorescence shown in figure 24 was in a spot where sulfates are identified by the SAM algorithm, which is an interesting point for remote sensing applied to astrobiology. Furthermore, efflorescence is described by \cite{ferrari2023}  as iron bearing sulfates (such as jarosite, copiapite, coquimbite). The presence of copiapite, for example, has been reported by \cite{roach2006} as a participant of this path of desiccation and later oxidation. MTMF result gives a good confirm to the SAM algorithm when mapping sulfates. In fact, it is possible to see the concentration of this mineral type along the riverbed mostly in the southern part, whereas clays lay on the slopes heading to the riverbed. The oxides mound found in the SAM image is coloured as sulfates in the MTMF result (fig. 25), but this is because oxides were assigned no colour during the RGB projection, and it was therefore assumed by the algorithm to be sulfate-like. When mapping LSU results (fig. 26), vegetated areas are well defined as hybrid colour (pink-violet), and the bare ground area of the PGZ shows sulfates according to SAM and MTMF results. 


```{r, label = "fig:23", echo=FALSE, fig.cap="Map 1 is a SAM classification after defining endmembers based on ROIs. Map 2 is a SAM classification of the mineral classes found with the methods previously described."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/AVIRIS MAP.png")
```

Clays surround the riverbed and are concentrated in the western side of the southern PGZ. Here our survey, according to \cite{ferrari2023} as well, found many exposed laminated sands and silts alternating with black and greyish layers. The exposure is important, as remote sensing nowadays relies on what the instrument can detect in the first micrometres of surface. It is in fact to be mentioned, that none of these are subsurface data. Sulfides mapping was considered unreliable, as this mineral class has an ambiguous spectrum that can give false positives, as it doesn’t have many features. Of course, in the area, sulfides such as pyrite, chalcopyrite and sphalerite were reported by \cite{hudsonedwards1999} and \cite{romero2006} and it was possible to observe samples of sulfides laying all over in the south-western branch of the riverbed. 
With such a high resolution, as for the AVIRIS-NG instrument, SAM technique is a good tool to identify minerals based on an endmember list supported by library spectra. It was able to classify pixels and identify interesting concentration trends. MTMF and LSU maps were important in this approach because they helped discard uninteresting spots, such has the vegetated areas, which are immediately visible as they are hybrid coloured and appear like a noisy background. Furthermore, these two techniques next to the SAM helped interpret data and discard false positives or strongly mixed materials, as it can be seen in the slopes heading to the riverbed, where clays are mixed with sulfates and oxides, and a predominant species can’t be identified. 


```{r, label = "fig:24", echo=FALSE, fig.cap="Popcorn-like structure (efflorescence) of bacteria found on steep clayey mounds in PGZ."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/POPCORN.jpg")
```



```{r, label = "fig:25", echo=FALSE, fig.cap="MTMF result mapping the three given classes on the RGB colours. See how the blue distribution of Sulphates recalls the SAM classification for the same class."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/MTMF.png")
```



```{r, label = "fig:26", echo=FALSE, fig.cap="LSU result mapping the three given classes on the RGB colours."}
knitr::include_graphics("F:/User/Documenti/R/thesis/thesis/LSU.png")
```

\newpage

# 4 CONCLUSIONS 
Without a proper ground-truthing it is impossible to have a confirmation of presence or absence of a given material or mineral. However, when speaking of planets others from ours, we can’t rely on a ground-truth itself. Thus, these techniques, so far based on sole remote sensing and spectral libraries, were used to draw a mineralogic map of interest to lead future in-situ missions on Earth for sample collection and ground-truthing. Results agree with literature findings.

## 4.1 Instruments
Even though ASTER radiometer SWIR channels burned out back to 2008, it is still possible to use previous datasets for some observations. The lower spectral resolving power is compensated by a high spatial resolution (considering it was launched in 1999) and a high manageability of its datasets. This makes multispectral instruments a very practical tool for an initial localization and for band-ratio analysis, which in our case showed the distribution of altered rocks and ferric iron and helped us focus on non-anthropized and non-vegetated areas for further studies.
Hyperspectral data from PRISMA carried out high resolution spectra for each pixel, but they were very compromised due to water vapour and were blurred by high overall noise. We had to cut out groups of bands to improve spectra for library comparison. These problems made it harder to spot reliable features at different wavelengths. Furthermore, the dataset is stored in Hierarchical Data Format 5, which many programs cannot open or have issues with keeping data integrity. The use of different steps to retrieve file info and data made the process less efficient and the results more dubious. A future project might skip this data if used for SAM mapping or might use it for a different scope. Panchromatic channel of PRISMA, instead, turned out to be very efficient for detailed observation of the location. It was useful for spotting land features and distinguishing anthropized or vegetated areas from exposed ground or rocky formations.
Hyperspectral airborne instruments such as AVIRIS-NG turned out to be very efficient, for atmosphere and water vapour interfered less, causing minor disturbs in the dataset, and noise was confined to few spectral wavelengths. Its high spatial and spectral resolving power were key for detailed data training and accuracy mapping of smaller areas, and this quality turns out to be very important when trying to spot smaller ecosystems and possible life niches. 

## 4.2 Mapping
It was important to lead a step-by-step analysis following an increase in spectral (and spatial) resolution. This has successfully let us delineate our areas of interest, where we could undisturbedly retrieve a surface mineral composition. Band-ratio map from ASTER data could give us exposed ground information straightforward. Data retrieved from PRISMA could apply an initial mineral characterization to the designated exposed ground. Finally, spectra extracted from airborne surveys of AVIRIS-NG could be compared to library spectra for mineral mapping of a smaller target location. 
Our analysis shows that sulphates are abundant in the whole mining district and particularly exposed in the PGZ, where they lay mainly in the riverbed. Sulphates, clays, and oxides presence in this area is also reported by field missions from Ferrari et al. \cite{ferrari2023} and it is no wonder that these alteration minerals are here copious. Predictions of band-ratios from ASTER dataset in the PGZ area show a mixture of alteration minerals and ferric iron (bluish-green colour in fig. 14), both compositions leading to clays, sulphates, and oxides, and whose results are confirmed by AVIRIS-NG based SAM map. Sulphides, instead, were very hard to detect. Their distribution was often associated to contaminated waters and tailings near ponds. This might be interesting to study with a different approach, for example, by analysing water samples and estimating the sulphide concentrations in the liquid, as well as the ratio of $Fe^{2+}/Fe^{3+}$, and then, for these samples, calculating a new model spectrum to be used in other analog situations. At last, the instruments, used at different spatial and spectral resolution, have turned out useful at different scales: on a smaller scale multispectral tools are skilful and handy to locate preferable targets for later deeper examinations. On a larger scale, hyperspectral data is fundamental for the deep analysis. This means identifying detailed mineral distribution, land cover, and all the critical elements for better addressing a human scaled mission, which is the final goal.

## 4.3 Future implications
All in all, the use of multiple instruments with different spatial and spectral resolutions surely gave us many advantages in this research, as it was seen, but at the same time had a higher cost in terms of time and money. In total, about six months of work was needed, including a field campaign on Rio Tinto. This was necessary to train the staff and acquire information through scientific literature and software and hardware manuals. The biggest part of the cost is due to programs distribution licenses. With the right knowledge this work can be executed quickly given a decent hardware. The increasing spread of open-source software (i.e. QGIS, EnMAP-Box) and open-access science gives these studies a chance to be published widely and to be tinkered to a higher accuracy, fostering a perpetual and constructive debate.

\newpage



\bibliography{mybib}
\bibliographystyle{apalike}