# Seaweeds in two oceans: beta-diversity
Supporting material for the paper: Smit AJ, Bolton JJ, Anderson RJ (in review) 
Seaweeds in two oceans: beta-diversity. Frontiers in Marine Science

## 1. Abstract
Several species assembly mechanisms have been proposed to structure
ecological communities. We assess the biogeography of seaweeds along
2,900 km of South Africa’s coastline in relation to a thermal gradient
produced by the Agulhas Current, and contrast this with the
environmental structure created by the Benguela Current. We subdivided
the coastline into ‘bioregions’ to examine the regional patterning. To
investigate the assembly mechanisms, we decomposed Sørensen’s β-diversity
into ‘turnover’ (β<sub>sim</sub>) and ‘nestedness-resultant’
(β<sub>sne</sub>) dissimilarities, and used distance-based redundancy
analysis (db-RDA) to relate them to the Euclidian thermal difference,
d<sub>E</sub>, and geographical distance. Moran’s eigenvector maps
(MEM) were used as an additional set of spatial constraints. Variation
partitioning was then used to find the relative strengths of thermal and
spatially-structured thermal drivers. Spatial and environmental
predictors explained 97.9% of the total variation in β<sub>sim</sub> and
the thermal gradient accounted for 84.2% of this combined pool.
β<sub>sim</sub> was the major component of overall β-diversity in the
Agulhas Current region, suggesting niche influences (environmental
sorting) as dominant assembly process there. The much weaker thermal
gradient in the Benguela Current-influenced region resulted in a high
amount of β<sub>sne</sub> that could indicate neutral assembly
processes. The intensification of upwelling during the mid-Pliocene
4.6–3.2 Ma (*i.e.* historical factors) were likely responsible for
setting up the strong disjunction between the species-poor west coast
and species-rich south and east coast floras, and this separation
continues to maintain two systems of community structuring mechanisms in
the Atlantic and Indian Ocean influenced sides of South Africa.

**The text below correspond to Appendices B and C of the paper.**

## 2. Appendix B

### 2.1 Spatial analysis background and code
The intention of this section is to show the approach and **R** scripts used to pull apart the spatial scales at which seaweed assemblages are structured around the coast of South Africa. Specifically, I wish to determine if these scales match those expressed by the coastal thermal provinces and the ocean regime underpinned by the Agulhas and Benguela Currents.

#### 2.1.1 The data
I use two data sets. The first, *Y*, comprises distribution records of 846 macroalgal species within each of 58 × 50 km-long sections (Appendix A) of the South African coast (updated from Bolton and Stegenga, 2002). This represents *ca*. 90% of the known seaweed flora of South Africa, but excludes some very small and/or very rare species for which data are insufficient. The data are from verifiable literature sources and John Bolton and Rob Anderson's own collections, assembled from information collected by teams of phycologists over three decades (Bolton, 1986; Bolton and Stegenga, 2002; De Clerck et al., 2005; Stegenga et al., 1997). The second, *E*, is a dataset of *in situ* coastal seawater temperatures (Smit et al., 2013) derived from daily measurements over up to 40 years.

A third data set of explanatory variables --- the spatial variables (*S*) --- is constructed as per the instructions in section *Preparation of spatial variables*, later on.

#### 2.1.2 Setting up the analysis environment
This is **R**, so first I need to find, install and load various packages. Some of the packages will be available on CRAN and can be accessed and installed in the usual way, but others will have to be downloaded from [R Forge](https://r-forge.r-project.org/R/?group_id=195).

``` r
library(betapart)
library(vegan)
library(gridExtra)
library(packfor)
library(grid)
library(gridBase)
library(tidyr)
library(spdep) # for dnearneigh() in PCNM.R
library(AEM) # for moran.I.multi() in PCNM.R
# install.packages("devtools")
# install.packages("AEM", repos = "http://R-Forge.R-project.org")
# install.packages("packfor", repos = "http://R-Forge.R-project.org")
source("functions/pcoa_all.R")
source("functions/PCNM.R")
source("functions/spatial_MEM.R")
```

Now I get to the data. The first step involves the species table (*Y*). First I compute the Sørensen dissimilarity and then I decompose the dissimilarity into the 'turnover' (β) and 'nestedness-resultant' (β) components (Baselga, 2010; Baselga et al., 2013) using the `betapart.core()` and `betapart.pair()` functions of the **betapart** package (Baselga et al., 2013). These are placed into the matrices *Y*1 and *Y*2. Optionally, I can apply a prinipal components analysis (PCA) on *Y* to find the major patterns in the community data. In **vegan** this is done using the `rda()` function and not supplying the constraints (*i.e.* the environment table, *E*, or the spatial table, *S*). The formal analysis will use the species data in distance-based redundancy analyses (db-RDA as per **vegan**'s `capscale()` function) by coupling them with *E* and *S*.

``` r
# Read in the species data (note: on GitHub only the distance
# matrices obtained via 'beta.part' and 'beta.pair' (below) 
# will be provided -- they are read in as 'Y1.Rdata' and 'Y2.Rdata';
# the raw data cannot be shared at this stage, but the distance matrix is provided):
# spp <- read.csv('../stats/seaweeds.csv')
# spp <- dplyr::select(spp, -1)

# Decompose total Sørensen dissimilarity into turnover and 
# nestedness-resultant components:
# Y.core <- betapart.core(spp) 
# Y.pair <- beta.pair(Y.core, index.family = "sor")

# Let Y1 be the turnover component (beta-sim):
# Y1 <- as.matrix(Y.pair<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/cba45b733617df1ac4c2a012b771e0b6.svg?invert_in_darkmode" align=middle width=582.057795pt height=47.64045000000003pt/>beta.sne)
# save(Y2, file = "data/Y2.Rdata")
load("data/Y2.Rdata")
```

It is now necessary to load the environmental data and some setup files that partition the 58 coastal sections (and the species and environmental data that fall within these sections) into bioregions.

The thermal (environmental) data contain various variables, but in the analysis I use only some of them. These data were obtained from many sites along the South African coast, but using interpolation (not included here) I calculated the thermal properties for each of the coastal sections for which seaweed data are available. Consequently we have a data frame with 58 rows and a column for each of the thermal metrics. Before use, I apply **vegan**'s `decostand()` function to scale the data to zero mean and unit variance.

Four bioregions are recognised for South Africa (Bolton and Anderson, 2004), namely the Benguela Marine Province (BMP; coastal sections **1**–**17**), the Benguela-Agulhas Transition Zone (B-ATZ; **18**–**22**), the Agulhas Marine Province (AMP; **19**–**43**/**44**) and the East Coast Transition Zone (ECTZ; **44**/**45**–**58**). My plotting functions partition the data into the bioregions and colour code the figures accordingly so I can see regional patterns in -diversity emerging.

``` r
# Now comes in the in situ temperatures for the 58 coastal sections 
# (interpolated temperaures as per version 2 of the South African Coastal Temperature Network):
load('data/E.RData')
env <- as.data.frame(interpOut)

# I select only some of the thermal vars; the rest
# are collinear with some of the ones I import:
E1 <- dplyr::select(env, febMean, febRange, febSD, augMean,
                    augRange, augSD, annMean, annRange, annSD)

# Calculate z-scores:
E1 <- decostand(E1, method = "standardize")

# Load the coordinates of the coastal sections:
sites <- read.csv("data/sites.csv")
sites <- sites[, c(2, 1)]

# Load the bioregion definition:
bioreg <- read.csv('data/bioregions.csv', header = TRUE)
```

#### 2.1.3 Preparation of spatial variables
I test the niche difference mechanism as the primary species compositional assembly process operating along South African shores. I suggest that the thermal gradient along the coast provides a suite of abiotic (thermal) conditions from which species can select based on their physiological tolerances, and hence this will structure -diversity. For this mechanism to function one would assume that all species have equal access to all sections along this stretch of coast, thus following Beijerinck’s 'Law' that everything is everywhere but the environment selects (Sauer, 1988) (but see main text!).

The basic approach to a spatial analysis structured around a biological response (*e.g.* community structure and composition; *Y*), environmental variables (*E*) and their spatial representation (*S*) involves an analysis of Moran's eigenvector maps (MEM), followed by db-RDA and variance partitioning. Various literature sources discuss principle behind Moran's eigenvector maps (Dray et al., 2006, 2012). Worked examples are also presented in the excellent book *Numerical Ecology with R* (Borcard et al., 2011) in Section 7.4. The method followed here has been adapted from these and other sources.

Obtaining the MEMs to use in the analysis is based on the procedure introduced by Borcard and Legendre (2002), which was later modified by Dray et al. (2006). The basic approach involves:

1.  Set up a geographic or Euclidian distance matrix representing the pairwise distances between the *n* sites (*D* = \[*d*<sub>*i**j*</sub>\]). I already did this when I applied the `decostand` function earlier.

2.  Construct a truncated distance matrix by calculating a Minimum Spanning Tree (*S*<sup>⋆</sup>) and noting the following rules:
    $$S^{\\star} =\\left\\{ \\begin{array}{rl} 0 & \\mbox{if}~i = j \\\\ d\_{ij} & \\mbox{if}~d\_{ij} \\leq t \\\\ 4t & \\mbox{if}~d\_{ij} &gt; t \\end{array} \\right.$$
     Weighting may be applied if desired, resulting in *S*<sub>*w*</sub><sup>⋆</sup>. It is not done here.

3.  Do a Principal Coordinates Analysis (PCoA) of the truncated distance matrix *S*<sup>⋆</sup>.

The spatial properties imprinted on the species and their environment can be specified using a matrix of Euclidian or geographic distances. These coordinates are 'truncated' into a square (section × section) matrix containing non-negative values (*S*<sup>⋆</sup>). By convention the diagonal values are set to zero. A very basic spatial matrix is binary, where 1 codes for pairs of neigbouring sites while 0 denotes non-connected sites according to the chosen network topology. Such matrices are called 'binary connectivity matrices' and relate to graphs made using distance criteria derived from graph theory.

Truncation produced by Minimum Spanning Trees (MST) focuses on the binary relationships between neighbouring sites, discarding any other connections (*i.e.* some sites are considered to be neighbours, while for others the relationships are null). One could also choose a Gabriel graph or another kind of network topology. Such matrix representations show section-to-section connectivities. In the case of South Africa's coastline data, the MST causes sections to be connected only to other sections adjacent to two sides of it: for example, Section **4** is directly connected to *only* Sections **3** and **5**; sections at the termini of the coastal 'string' of sections are each connected to only one other section. The binary connectivity matrices, also called *topology-based connectivity matrices*, can be produced from Euclidian or geographic coordinates using functions in at least two **R** packages (I start with geographic coordinates). One option is to use the **spdep** package's `mst.nb()` function to calculate a MST, but there are also options in the **vegan** package and elsewhere. The neighbours list arrived at from the MST represents the spatial component, *S*<sup>⋆</sup>. The MST results in small connectivity artefacts in the Saldanha Bay region where the closest sections are not necessarily the ones adjacent one-another following along the path around the coast, because sections at opposite sides of the bay may in fact be closer together. This topological inconsistency does not affect the spatial analysis in any way.

Once the truncated distance matrix has been prepared, it is subjected to a PCoA and I keep the eigenvectors that represent positive spatial correlation (positive Moran’s *I*). For the MEM analysis I use the function `PCNM()` that resides in the `functions` folder in the file `PCNM.R` (see notes inside about authorship). PCNM stands for Principal Coordinates Analysis of Neighbourhood Matrices (the neighbourhood matrix in this instance being the MST). This method automatically constructs the spatial variables and calculates the Moran’s I for each. The MEMs are completely orthogonal and represent the spatial structures over the full range of scales from 50 to 2,700 km. The large eigenvectors represent broad spatial scales while smaller ones cover finer features. The *spatial data* will be used as a set of explanatory variables in the multiple regression type analyses applied to a species dissimilarity matrix (*i.e.* the db-RDA; Dray et al., 2012)

The code below reproduces the spatial analysis in the paper. Due to the length of the output I have prevented the script from returning any output here; rather, if the reader is for some odd reason interested in repeating this analysis, s/he may find the data and scripts in my [GitHub](https://github.com/ajsmit/Seaweed-beta) repository, and the full code can be run in its entirety. Well, I hope this will work, but if it doesn't (probably very likely) then write to me at <ajsmit@uwc.ac.za> and I shall assist --- this may depend on if your email has a catchy title that will make it stand out from among all the other emails which other people think are equally important.

``` r
## Auto PCNM:
S.auto <- PCNM(dist(sites), silent = TRUE)
# summary(S.auto)

# The truncation distance:
S.dmin <- S.auto<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/485a8fafd9b477b4153f74faee24ef6d.svg?invert_in_darkmode" align=middle width=346.90474499999993pt height=47.64044999999998pt/>values)

# Expected value of I, no spatial correlation:
S.auto<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/8265d71a0f23ca79144517ba8a094ae5.svg?invert_in_darkmode" align=middle width=525.0519449999999pt height=47.64044999999998pt/>Moran_I<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/c2f5f8a5bacd1481e04896e2e8dca34d.svg?invert_in_darkmode" align=middle width=920.7186449999999pt height=47.64044999999998pt/>vectors)[, S.sel]
```

The code below lets us visualise the configuration of the 58 coastal sections as represented by the minimum spanning tree. Because the sites are constrained by the coast the MST network topology results in a string of coastal sections arranged along the shore between Section **1** and Section **58**. This spatial network therefore also captures the spatial connectivity in the seaweed's dispersal ability along the shore, although no directionality is associated with dispersal. In the paper I discuss the possible influence of ocean currents (*e.g.* Wernberg et al., 2013) and I pointed out that it is tempting to assume that seaweeds would disperse in the direction the major ocean currents. These kinds of networks could conceivably be configured to model dispersal due to currents, but here it is simply used for representing the spatial scale of the study region.

``` r
# The spatial netwwork topology of the coastal sections can be seen by:
plot(S.auto<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/8128718e75f5aea33886cdd54a54da26.svg?invert_in_darkmode" align=middle width=698.16285pt height=363.25707pt/>adj.r.squared

# Variance explained by full model:
sum(S.Y1.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / S.Y1.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/2fd123cea69aa41dcf6dbc816f65b0e1.svg?invert_in_darkmode" align=middle width=698.0968499999999pt height=87.09227999999999pt/>adj.r.squared
sum(S.Y2.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / S.Y2.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/69b24062d2865995197f33a9cc6274be.svg?invert_in_darkmode" align=middle width=745.4271pt height=915.58599pt/>adj.r.squared

# Variance explained by reduced model:
sum(S.Y1.s2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / S.Y1.s2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/31ff202a39c3cf48a2dc89bb2af0cc25.svg?invert_in_darkmode" align=middle width=784.2977999999999pt height=363.25838999999996pt/>adj.r.squared

sum(S.Y2.s2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / S.Y2.s2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/49585c24a3c004e0aafafc03631af8ab.svg?invert_in_darkmode" align=middle width=842.80185pt height=955.0374899999999pt/>constraints)
S.Y1.df_sites<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/1708cda66816858736a4bec1e591c942.svg?invert_in_darkmode" align=middle width=123.566685pt height=22.745910000000016pt/>bolton
S.Y1.df_sites<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/988e7e72c7dc19c17209ba8257a0c445.svg?invert_in_darkmode" align=middle width=618.783495pt height=47.64045000000003pt/>biplot)
S.Y1.bp <- S.Y1.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/e28007f90984ddfba4f3436a0147cc61.svg?invert_in_darkmode" align=middle width=450.74749499999996pt height=24.56552999999997pt/>labels <- rownames(S.Y1.bp)
colnames(S.Y1.bp) <- c("x", "y", "labels")
S.Y1.bp.sign <- S.Y1.bp[S.Y1.bp<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/731465cc44cf4c2c76504e170a0b6eaf.svg?invert_in_darkmode" align=middle width=698.06385pt height=126.54443999999998pt/>biplot)
S.Y1.text <- as.data.frame(S.Y1.text)
S.Y1.text<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/e8c128f8d461f559eea4d41b68fdd1ee.svg?invert_in_darkmode" align=middle width=717.80115pt height=47.64045000000003pt/>labels %in% S.Y1.sign.ax,]

S.Y1.p <- ggplot(data = S.Y1.df_sites, aes(x, y, colour = Bioregion)) + 
  geom_point(size = 4.0) + 
  geom_text(aes(label = Section), size = 3.0, col = "white") + 
  geom_segment(data = S.Y1.bp, 
               aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "red", alpha = 1, size = 0.7) +
  geom_text(data = as.data.frame(S.Y1.text), 
            aes(x, y, label = rownames(S.Y1.text)),
            color = "black") +
  xlab("CAP1") + ylab("CAP2") + 
  ggtitle(expression(paste("Spatial variables and ", beta[sim]))) +
  theme_grey() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 0.8)

# ...then for Y2:
S.Y2.scrs <- scores(S.Y2.s2, display = c("sp","wa","lc","bp","cn"))
S.Y2.df_sites <- data.frame(S.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/cc8ab4b25a89842e09440b65651973d5.svg?invert_in_darkmode" align=middle width=178.58989499999998pt height=24.56552999999997pt/>bioreg <- bioreg<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/eef6a122069e669bd107ed6468a63279.svg?invert_in_darkmode" align=middle width=132.069795pt height=22.745910000000016pt/>section <- seq(1:58)
colnames(S.Y2.df_sites) <- c("x", "y", "Bioregion", "Section")

multiplier <- ordiArrowMul(S.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/c8d302c66b8ed21ed45371ed27d469f0.svg?invert_in_darkmode" align=middle width=291.67429500000003pt height=24.56552999999997pt/>biplot * multiplier
S.Y2.bp <- as.data.frame(S.Y2.bp)
S.Y2.bp<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/b33fa3ceb5070a5c33b3c9c8a8cc5717.svg?invert_in_darkmode" align=middle width=697.9896pt height=47.64045000000003pt/>labels %in% S.Y2.sign.ax,]

S.Y2.text <- text.mult(S.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/cede1dbca00ed3d10ff512244f768fad.svg?invert_in_darkmode" align=middle width=491.297895pt height=24.56552999999997pt/>labels <- rownames(S.Y2.text)
colnames(S.Y2.text) <- c("x", "y", "labels")
S.Y2.text.sign <- S.Y2.text[S.Y2.text<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/e9fea300d44c9353ad60ba286c3f6fee.svg?invert_in_darkmode" align=middle width=844.3825499999999pt height=718.32519pt/>adj.r.squared

# Variance explained by full model:
sum(E.Y1.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / E.Y1.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/33ffe17bfc0c23ab37598e2cb4d91e87.svg?invert_in_darkmode" align=middle width=698.30145pt height=87.09227999999999pt/>adj.r.squared
sum(E.Y2.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / E.Y2.cs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/522912787703bd81421dfb21519c834a.svg?invert_in_darkmode" align=middle width=764.0786999999999pt height=678.8736899999999pt/>adj.r.squared

# Variance explained by reduced (final) model:
sum(E.Y1.cs2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / E.Y1.cs2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/69679b76c171387c44a09a8a7748ba61.svg?invert_in_darkmode" align=middle width=697.9665pt height=599.97069pt/>adj.r.squared

sum(E.Y2.cs2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/987b5bdff7f72fc2ac32ee23c1aa6b8f.svg?invert_in_darkmode" align=middle width=38.044215pt height=22.381919999999983pt/>eig) / E.Y2.cs2<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/b4bbb0e4d311d5da576663f45f0a02b5.svg?invert_in_darkmode" align=middle width=1100.7711pt height=639.4221899999999pt/>constraints)
E.Y1.df_sites<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/1708cda66816858736a4bec1e591c942.svg?invert_in_darkmode" align=middle width=123.566685pt height=22.745910000000016pt/>bolton
E.Y1.df_sites<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/95b6b5f26c76761837cba61b6186e94e.svg?invert_in_darkmode" align=middle width=621.743595pt height=47.64045000000003pt/>biplot)
E.Y1.bp <- E.Y1.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/bd2ed01b3d4bdc2e3e30a72b13573195.svg?invert_in_darkmode" align=middle width=459.62944500000003pt height=24.56552999999997pt/>labels <- rownames(E.Y1.bp)
colnames(E.Y1.bp) <- c("x", "y", "labels")
E.Y1.bp.sign <- E.Y1.bp[E.Y1.bp<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/df6d4831f399828457f5ab0ba23551c9.svg?invert_in_darkmode" align=middle width=222.26539499999998pt height=47.64044999999998pt/>biplot)
E.Y1.text <- as.data.frame(E.Y1.text)
E.Y1.text<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/c1d0bcd08a17a3e176857d250eff23c5.svg?invert_in_darkmode" align=middle width=726.6814499999999pt height=47.64045000000003pt/>labels %in% E.Y1.sign.ax,]

E.Y1.p <- ggplot(data = E.Y1.df_sites, aes(x, y, colour = Bioregion)) + 
  geom_point(size = 4.0) + 
  geom_text(aes(label = Section), size = 3.0, col = "white") + 
  geom_segment(data = E.Y1.bp.sign, 
               aes(x = 0, y = 0, xend = x, yend = y),
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "red", alpha = 1, size = 0.7) +
  geom_text(data = as.data.frame(E.Y1.text.sign), 
            aes(x, y, label = rownames(E.Y1.text.sign)),
            color = "black") +
  xlab("CAP1") + ylab("CAP2") + 
  ggtitle(expression(paste("Thermal variables and ", beta[sim]))) +
  theme_grey() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "none",
        aspect.ratio = 0.8)

E.Y2.scrs <- scores(E.Y2.cs2, display = c("sp","wa","lc","bp","cn"))
E.Y2.df_sites <- data.frame(E.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/e78233f85b1724185526cf50297cd753.svg?invert_in_darkmode" align=middle width=181.549995pt height=24.56552999999997pt/>bioreg <- bioreg<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/53e11ab883482b3480101d2275f2e3b9.svg?invert_in_darkmode" align=middle width=135.03006pt height=22.745910000000016pt/>section <- seq(1:58)
colnames(E.Y2.df_sites) <- c("x", "y", "Bioregion", "Section")

multiplier <- ordiArrowMul(E.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/dd25e1b68c528cf146d88992314a036c.svg?invert_in_darkmode" align=middle width=297.594495pt height=24.56552999999997pt/>biplot * multiplier
E.Y2.bp <- as.data.frame(E.Y2.bp)
E.Y2.bp<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/d5c6e98e51d8e77acffd0a395a92beeb.svg?invert_in_darkmode" align=middle width=697.9665pt height=47.64045000000003pt/>labels %in% E.Y2.sign.ax,]

E.Y2.text <- text.mult(E.Y2.scrs<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/bcfec39573c37379f646485952df2e4c.svg?invert_in_darkmode" align=middle width=500.17984499999994pt height=24.56552999999997pt/>labels <- rownames(E.Y2.text)
colnames(E.Y2.text) <- c("x", "y", "labels")
E.Y2.text.sign <- E.Y2.text[E.Y2.text<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/c4fd1629893e1912c0c1ef2ec06d48ff.svg?invert_in_darkmode" align=middle width=850.3027500000001pt height=3866.3040899999996pt/>Y1 <- as.vector(Y1)

Y1.sl.BMP <- Y1.sl %>%
  dplyr::left_join(., sites, by = c("to" = "site")) %>% 
  dplyr::left_join(., sites, by = c("from" = "site")) %>% 
  dplyr::filter(Y1 <= 0.5 & Y1 != 0) %>% 
  dplyr::filter(from != to & from <= 16)

Y1.sl.BATZ <- Y1.sl %>%
  dplyr::left_join(., sites, by = c("to" = "site")) %>% 
  dplyr::left_join(., sites, by = c("from" = "site")) %>% 
  dplyr::filter(Y1 <= 0.5 & Y1 != 0) %>% 
  dplyr::filter(from != to & from > 16 & from <= 21)

Y1.sl.AMP <- Y1.sl %>%
  dplyr::left_join(., sites, by = c("to" = "site")) %>% 
  dplyr::left_join(., sites, by = c("from" = "site")) %>% 
  dplyr::filter(Y1 <= 0.5 & Y1 != 0) %>% 
  dplyr::filter(from != to & from > 21 & from <= 41)

Y1.sl.ECTZ <- Y1.sl %>%
  dplyr::left_join(., sites, by = c("to" = "site")) %>% 
  dplyr::left_join(., sites, by = c("from" = "site")) %>% 
  dplyr::filter(Y1 <= 0.5 & Y1 != 0) %>% 
  dplyr::filter(from != to & from > 41)

Y1.sl <- rbind(Y1.sl.BMP, Y1.sl.BATZ, Y1.sl.AMP, Y1.sl.ECTZ)

# and then Y2:
Y2.sl<img src="https://rawgit.com/ajsmit/Seaweeds_in_Two_Oceans/master/svgs/b20c789550be056ba9a06bf00b5b0e9c.svg?invert_in_darkmode" align=middle width=732.5868pt height=640.62207pt/>bolton), 
               col = "black", shape = 21) +
    scale_fill_manual(breaks = c("AMP", "B-ATZ", "BMP", "ECTZ"),
                      values = col.seq, name = "bioregion", guide = FALSE) +
    scale_colour_gradient(name = expression(paste(beta[sim])), 
                          low = "black", high = "red") +
    coord_fixed(ratio = 1, expand = TRUE) +
    scale_x_continuous(breaks = seq(15, 35, 5),
                       labels = scales::unit_format("°E", sep = "")) +
    scale_y_continuous(breaks = seq(-35, -25, 5),
                       labels = c("35°S", "30°S", "25°S")) +
    scale_alpha_continuous(guide = FALSE) +
    theme_grey() + xlab(NULL) + ylab(NULL) +
    theme(panel.grid.minor = element_blank()) +
    ggtitle(title)}

a <- net.plot.Y1(Y1.sl.BMP, "Benguela Marine Province", 
                 col.seq = c("black", "black", "white", "black")) + # alphabetical
  theme(legend.direction = "horizontal",
        legend.position = c(x = 0.5, y = 0.8),
        legend.key.height = unit(0.3, "cm")) 
b <- net.plot.Y1(Y1.sl.BATZ, "Benguela-Agulhas Transition Zone", 
                 col.seq = c("black", "white", "black", "black")) +
  theme(legend.position = "none")
c <- net.plot.Y1(Y1.sl.AMP, "Agulhas Marine Province", 
                 col.seq = c("white", "black", "black", "black")) +
  theme(legend.position = "none")
d <- net.plot.Y1(Y1.sl.ECTZ, "East Coast Transition Zone", 
                 col.seq = c("black", "black", "black", "white")) +
  theme(legend.position = "none")

# pdf("Fig4.pdf", width = 7, height = 3.8)
# grid::grid.newpage()
# grid::pushViewport(grid::viewport(layout = grid::grid.layout(2,2)))
# vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)
# print(a, vp = vplayout(1,1))
# print(b, vp = vplayout(1,2))
# print(c, vp = vplayout(2,1))
# print(d, vp = vplayout(2,2))
# dev.off()
```

And that's it, folks. You'll notice that I haven't reproduced Fig. 5 here. I'll leave that up to you... or ask me and I'll send the code. All of the matrices have (mostly) been calculated above and they can be used together with some **dplyr** and **ggplot2** know-how in the creation of that graph.

Legalise seaweed!

## 4. References
Anderson, M. J., Tolimieri, N., and Millar, R. B. (2013). Beta diversity of demersal fish assemblages in the North-Eastern Pacific: interactions of latitude and depth. *PLOS ONE* 8, e57918.

Baselga, A. (2010). Partitioning the turnover and nestedness components of beta diversity. *Global Ecology and Biogeography* 19, 134–143.

Baselga, A., Orme, D., Villeger, S., Bortoli, J. D., and Leprieur, F. (2013). *betapart: Partitioning beta diversity into turnover and nestedness components*. Available at: <http://CRAN.R-project.org/package=betapart>.

Blanchet, F. G., Legendre, P., and Borcard, D. (2008). Forward selection of explanatory variables. *Ecology* 89, 2623–2632.

Bolton, J. J. (1986). Marine phytogeography of the Benguela upwelling region on the west coast of southern Africa: A temperature dependent approach. *Botanica Marina* 29, 251–256.

Bolton, J. J., and Anderson, R. J. (2004). “Marine Vegetation,” in *Vegetation of southern africa*, eds. R. M. Cowling, D. M. Richardson, and S. M. Pierce (Cambridge University Press), 348–370.

Bolton, J. J., and Stegenga, H. (2002). Seaweed species diversity in South Africa. *South African Journal of Marine Science* 24, 9–18.

Borcard, D., Gillet, F., and Legendre, P. (2011). *Numerical Ecology with R*. Springer New York Available at: <https://books.google.co.za/books?id=dtQNxsH4Y2wC>.

Davidar, P., Rajagopal, B., Mohandass, D., Puyravaud, J.-P., Condit, R., Wright, S., et al. (2007). The effect of climatic gradients, topographic variation and species traits on the beta diversity of rain forest trees. *Global Ecology and Biogeography* 16, 510–518.

De Clerck, O., Bolton, J. J., Anderson, R. J., and Coppejans, E. (2005). Guide to the seaweeds of KwaZulu-Natal. *Scripta Botanica Belgica* 33, 294 pp.

Dray, S., Legendre, P., and Peres-Neto, P. R. (2006). Spatial modelling: a comprehensive framework for principal coordinate analysis of neighbour matrices (PCNM). *Ecological Modelling* 196, 483–493.

Dray, S., Pélissier, R., Couteron, P., Fortin, M. J., Legendre, P., Peres-Neto, P. R., et al. (2012). Community ecology in the age of multivariate multiscale spatial analysis. *Ecological Monographs* 82, 257–275.

Nekola, J. C., White, P. S., Carolina, N., Nekola, C., and Curriculum, P. S. W. (1999). The distance decay of similarity in biogeography and ecology. *Journal of Biogeography* 26, 867–878.

Peres-Neto, P. R., and Legendre, P. (2010). Estimating and controlling for spatial structure in the study of ecological communities. *Global Ecology and Biogeography* 19, 174–184.

Peres-Neto, P. R., Legendre, P., Dray, S., and Borcard, D. (2006). Variation partitioning of species data matrices: estimation and comparison of fractions. *Ecology* 87, 2614–2625.

Sauer, J. D. (1988). *Plant migration: The dynamics of geographic patterning in seed plant species*. University of California Press.

Smit, A. J., Roberts, M., Anderson, R. J., Dufois, F., Dudley, S. F. J., Bornman, T. G., et al. (2013). A coastal seawater temperature dataset for biogeographical studies: large biases between *in situ* and remotely-sensed data sets around the coast of South Africa. *PLOS ONE* 8, e81944.

Stegenga, H., Bolton, J. J., and Anderson, R. J. (1997). Seaweeds of the South African west coast. *Contributions of the Bolus Herbarium* 18, 3–637.

Wernberg, T., Thomsen, M. S., Connell, S. D., Russell, B. D., Waters, J. M., Zuccarello, G. C., et al. (2013). The footprint of continental-scale ocean currents on the biogeography of seaweeds. *PLOS ONE* 8, e80168.

