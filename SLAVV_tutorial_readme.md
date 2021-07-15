# **Graphical Curation Interface Tutorial ReadMe**
## Segmentation-Less, Automated, Vascular Vectorization (SLAVV)
The SLAVV software ([supporting functions](https://github.com/UTFOIL/Vectorization-Public/blob/master/source)) gives the user the ability to curate the vectorization output with live visual feedback. This graphical curator interface automatically opens at the appropriate times during the default execution of the [main function](https://github.com/UTFOIL/Vectorization-Public/blob/master/vectorize_V200.m). 

### SLAVV Methodology Overview
The SLAVV method consists of four workflow steps as described in the [Methodology Manuscript](https://github.com/UTFOIL/Vectorization-Public#Methodology-Manuscript) section: Automated Vascular Vectorization:
1. linear "energy" filtering
2. vertex extraction
3. edge extraction
4. network organization

These steps must be executed sequentially, starting with the "energy" step. The input parameters for each specific workflow are described in the [Documentation](https://github.com/UTFOIL/Vectorization-Public#Documentation) section: _Workflow-Specific Parameters_.

The graphical curator interface opens after the vertex and edge steps by default as described in the [Documentation](https://github.com/UTFOIL/Vectorization-Public#Documentation) section: Logistical Parameters.

## Tutorial Contents
This tutorial demonstrates the use of the graphical curator interface on three large (~1 mm<sup>3</sup>), images of living adult mouse brain microvasculature. The three images (**Image A**, **B**, and **C**) are all from the same mouse, approximately the same field of view, and timed 2 weeks apart. 
### Table of Contents
- [Input Images](#Input-images)
- [Energy Images](#Energy-images)
- [Vertex Curator](#Vertex-Curator)
  - [Global Thresholding](#Vertex-Global-thresholding)
  - [Local Thresholding](#Vertex-Local-Thresholding)
  - [Vertex Toggling](#Vertex-Toggling)
  - [Final Curated Vertex Set](Final-curated-Vertex-set)
- [Edge Curator](Edge-Curator)
  - [Global Thresholding](#edge-Global-Thresholding)
  - [Local Thresholding](#edge-Local-Thresholding)
  - [Orthogonal Views](#edge-orthogonal-views)
  - [Edge Toggling](#edge-toggling)
  - [Final Curated Edge Set](Final-curated-edge-set)
- [Visualizations](#visualizations)
<!---- [SLAVV Methodology Overview](#methodology-overview)
- [Inputs and Parameters](#Inputs-and-Parameters)
  - [Wrapper Script](#Wrapper-script)
- [Vertex Extraction](#Vertex-Extraction)-->

<!--The SLAVV methodology is described in the [Methodology Manuscript](https://github.com/UTFOIL/Vectorization-Public#Methodology-Manuscript). -->

### In the following screenshots, depictions of: 
**Image A** are in the left column | **Image B** the center | and **Image C** the right
:--------------------------------- | :--------------------: | ------------------------:

<!--#### [This folder](https://github.com/UTFOIL/Vectorization-Public/blob/master/tutorial/) contains the screenshots demonstrating the use of the SLAVV software on **Images A**, **B**, and **C**.-->
<!--
## Inputs and Parameters
-->
## Input Images
There is no pre-processing (interpolation, filtering, etc.) required. Simply input a raw TIF from a file location or a matrix from the workspace, as described in the [Documentation](https://github.com/UTFOIL/Vectorization-Public#Documentation), section: Optional Input.

#### Maximum intensity projections of the first 100 microns of each of the three **input** image stacks:
A|B|C
:---: | :---: | :---:
![](tutorial/000.png)|![](tutorial/000b.png)|![](tutorial/000c.png)

## Energy Images
#### _Minimum_ intensity projection of the first 100 microns of each the three **energy-filtered** images:
A|B|C
:---: | :---: | :---:
![](tutorial/001.png)|![](tutorial/001b.png)|![](tutorial/001c.png)

<!--Images 1-2:
### Wrapper Script
Image A Input Selection | Starting Workflow Step Selection
:---: | :---:
![](tutorial/1.png)|![](tutorial/2.png)-->

<!--Images 10's:-->
## Vertex Curator
Vertices have both point-location and radial-sizing components, and should be considered true (blue color in curator) when both the location and size match a vessel in the underlying original image. The user has the ability to curate the vertices using local thresholding as well as individual/group toggling. 
### Vertex Global thresholding
A|B|C
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
### Vertex Local Thresholding
A|B|C
:---: | :---: | :---:
![](tutorial/20.png)|![](tutorial/20b.png)|![](tutorial/20c.png)
![](tutorial/21.png)|![](tutorial/21b.png)|
![](tutorial/22.png)|![](tutorial/22b.png)|

### Vertex Toggling
|A|B|C
| :---: | :---: | :---:
||![](tutorial/23b.png)|
||![](tutorial/24b.png)|
||![](tutorial/25b.png)|
||![](tutorial/26b.png)|

<!--Image  30's: -->
### Final Curated Vertex Set
To be passed to the Edge Tracing step.
|A|B|C
| :---: | :---: | :---:
|![](tutorial/30.png)|![](tutorial/30b.png)|
||![](tutorial/31b.png)|

<!--Images 40's:-->
## Edge Curator
### Edge Global Thresholding
|A|B|C
| :---: | :---: | :---:
|![](tutorial/40.png)|![](tutorial/40b.png)|![](tutorial/40c.png)
|![](tutorial/41.png)|![](tutorial/41b.png)|![](tutorial/41c.png)
|![](tutorial/42.png)|![](tutorial/42b.png)|![](tutorial/42c.png)
|![](tutorial/43.png)|![](tutorial/43b.png)|![](tutorial/43c.png)
|![](tutorial/44.png)|![](tutorial/44b.png)|![](tutorial/45c.png)
|![](tutorial/45.png)|![](tutorial/45b.png)|![](tutorial/46c.png)
|![](tutorial/46.png)|![](tutorial/46b.png)|
|![](tutorial/47.png)|![](tutorial/47b.png)|

### Edge Local Thresholding
|A|B|C
| :---: | :---: | :---:
|||![](tutorial/44c.png)

### Orthogonal Views
A|B|C
:---: | :---: | :---:
![](tutorial/50.png)||
![](tutorial/51.png)||
![](tutorial/52.png)||
![](tutorial/53.png)||
![](tutorial/54.png)||
### Edge Addition  
A|B|C
:---: | :---: | :---:
![](tutorial/60.png)||
![](tutorial/61.png)||
![](tutorial/62.png)||
![](tutorial/63.png)||
### Edge Toggling
![](tutorial/70.png)
![](tutorial/71.png)
### Final Curated Edge Set 
#### MIP over the whole volume with graded edge color weighting.  
|A|B|C
| :---: | :---: | :---:
||![](tutorial/79b.png)|

### Maximum Intensity Projection Outputs 
(Coloring: strand uniques, depth, and direction) from the middle (1/3 in all 3 dimensions) of the vectors overlaying the image.  
A|B|C
:---: | :---: | :---:
![](tutorial/80.png)|![](tutorial/80b.png)|![](tutorial/80c.png)

### Histograms 
of various statistics of interest from the vectors in the image.  
A|B|C
:---: | :---: | :---:
![](tutorial/81.png)|![](tutorial/81b.png)|![](tutorial/81c.png)

### Visualizations 
of the output vectors using a .vmv file output to the VessMorphoVis plugin to Blender.  

![](tutorial/90.png)

A|B|C
:---: | :---: | :---:
![](tutorial/91.png)|![](tutorial/90b.png)|![](tutorial/90c.png)
![](tutorial/92.png)|![](tutorial/91b.png)|![](tutorial/91c.png)
![](tutorial/93.png)|![](tutorial/92b.png)|![](tutorial/92c.png)
![](tutorial/94.png)|![](tutorial/93b.png)|![](tutorial/93c.png)
![](tutorial/95.png)|![](tutorial/94b.png)|![](tutorial/94c.png)
![](tutorial/96.png)|![](tutorial/95b.png)|![](tutorial/95c.png)

A|B|C
:---: | :---: | :---:
![](tutorial/97.png)|![](tutorial/96b.png)|
![](tutorial/98.png)||
![](tutorial/99.png)||
![](tutorial/100.png)||
![](tutorial/101.png)||
![](tutorial/102.png)||
![](tutorial/103.png)||
