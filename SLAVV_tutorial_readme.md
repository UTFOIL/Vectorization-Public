# **Interactive Curator Tutorial ReadMe**
## Segmentation-Less, Automated, Vascular Vectorization (SLAVV)
#### Software ([source in MATLAB](https://github.com/UTFOIL/Vectorization-Public/blob/master/source)) Tutorial
Demonstrating the use of the graphical curator interface (built-into the SLAVV software) on three large (~1 mm<sup>3</sup>), images of living adult mouse brain microvasculature. The three images (**Image A**, **B**, and **C**) are all from the same mouse, approximately the same field of view, and timed 2 weeks apart.

#### In the screenshots to follow in this tutorial, depictions of: 
**Image A** are in the left column | **Image B** the center | and **Image C** the right
:--------------------------------- | :--------------------: | ------------------------:

## SLAVV Methodology Overview

The SLAVV method consists of four workflow steps:
1. linear "energy" filtering
2. vertex extraction
3. edge extraction
4. network organization

### Table of Contents:
- [SLAVV Methodology Overview](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#slavv-methodology-overview)
1. [Inputs and Parameters](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#Inputs-and-Parameters)
  1. [Wrapper Script](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#slavv-Wrapper-script)
2. [Vertex Extraction](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#slavv-Vertex-Extraction)
2. [Vertex Curator](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#slavv-Vertex-Curator)
  1. [Global Threshold Selection](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md#Global-threshold-selection)
 
<!--#### [This folder](https://github.com/UTFOIL/Vectorization-Public/blob/master/tutorial/) contains the screenshots demonstrating the use of the SLAVV software on **Images A**, **B**, and **C**.-->

## Inputs and Parameters

### Input Images

<!--Images 1-2:-->
### Wrapper Script
Image A Input Selection | Starting Workflow Step Selection
:---: | :---:
![](tutorial/1.png)|![](tutorial/2.png)-->

<!--Images 10's:-->
## Vertex Curator
### Global threshold selection
Image A | B | C
:---: | :---: | :---:
some_text | some_text | some_text
![](tutorial/10.png)||
![](tutorial/11.png)|![](tutorial/10b.png)|![](tutorial/10c.png)
![](tutorial/12.png)|![](tutorial/11b.png)|![](tutorial/11c.png)
![](tutorial/14.png)|![](tutorial/12b.png)|![](tutorial/12c.png)
![](tutorial/15.png)|![](tutorial/13b.png)|
![](tutorial/16.png)||
![](tutorial/17.png)||

<!--Images 20's:--> 
### Local  Threshold Selection
A | B | C
:---: | :---: | :---:
![](tutorial/20.png)|![](tutorial/20b.png)|![](tutorial/20c.png)
![](tutorial/21.png)|![](tutorial/21b.png)|
![](tutorial/22.png)|![](tutorial/22b.png)|

### Individual Vertex Toggling
| A | B | C
| :---: | :---: | :---:
||![](tutorial/23b.png)|
||![](tutorial/24b.png)|
||![](tutorial/25b.png)|
||![](tutorial/26b.png)|

<!--Image  30's: -->
### Vertex Curator: Final curated Vertex set for edge tracing  
3XA
![](tutorial/30.png)
3XB
![](tutorial/30b.png)
![](tutorial/31b.png)
### Images 40's:   Edge Curator: Volume navigation, exploration, and (for "B" and "C" images) local thresholding.  
4XA
![](tutorial/40.png)
![](tutorial/41.png)
![](tutorial/42.png)
![](tutorial/43.png)
![](tutorial/44.png)
![](tutorial/45.png)
![](tutorial/46.png)
![](tutorial/47.png)
4XB
![](tutorial/40b.png)
![](tutorial/41b.png)
![](tutorial/42b.png)
![](tutorial/43b.png)
![](tutorial/44b.png)
![](tutorial/45b.png)
![](tutorial/46b.png)
![](tutorial/47b.png)
4XC
![](tutorial/40c.png)
![](tutorial/41c.png)
![](tutorial/42c.png)
![](tutorial/43c.png)
![](tutorial/44c.png)
![](tutorial/45c.png)
![](tutorial/46c.png)
### Images 50's:   Edge Curator: Orthogonal Views  
5XA
![](tutorial/50.png)
![](tutorial/51.png)
![](tutorial/52.png)
![](tutorial/53.png)
![](tutorial/54.png)
### Images 60's:   Edge Curator: Addition Tool  
6XA
![](tutorial/60.png)
![](tutorial/61.png)
![](tutorial/62.png)
![](tutorial/63.png)
### Images 70's:   Edge Curator: Toggling Tool  
7XA
![](tutorial/70.png)
![](tutorial/71.png)
### Image  79b :   Edge Curator: Final output MIP over the whole volume with graded edge color weighting.  
79B
![](tutorial/79b.png)

### Image    80: Maximum Intensity Projection Outputs (Coloring: strand uniques, depth, and direction) from the middle (1/3 in all 3 dimensions) of the vectors overlaying the image.  
8XA
![](tutorial/80.png)
8XB
![](tutorial/80b.png)
8XC
![](tutorial/80c.png)

### Image    81: Histograms of various statistics of interest from the vectors in the image.  
81A  
![](tutorial/81.png)
81B  
![](tutorial/81b.png)
81C
![](tutorial/81c.png)

### Images  90+: Visualizations of the output vectors using a .vmv file output to the VessMorphoVis plugin to Blender.  
9XA and 10XA  
![](tutorial/90.png)
![](tutorial/91.png)
![](tutorial/92.png)
![](tutorial/93.png)
![](tutorial/94.png)
![](tutorial/95.png)
![](tutorial/96.png)
![](tutorial/97.png)
![](tutorial/98.png)
![](tutorial/99.png)
![](tutorial/100.png)
![](tutorial/101.png)
![](tutorial/102.png)
![](tutorial/103.png)
9XB  
![](tutorial/90b.png)
![](tutorial/91b.png)
![](tutorial/92b.png)
![](tutorial/93b.png)
![](tutorial/94b.png)
![](tutorial/95b.png)
![](tutorial/96b.png)
9XC  
![](tutorial/90c.png)
![](tutorial/91c.png)
![](tutorial/92c.png)
![](tutorial/93c.png)
![](tutorial/94c.png)
![](tutorial/95c.png)
