#Xwalk
(version 0.6)

##INTRODUCTION
Chemical cross-linking of proteins or protein complexes and the mass spectrometry based localization of the cross-linked amino acids is a powerful method for generating distance information on the substrate topology. Here we introduce the algorithm Xwalk for predicting and validating these cross-links on existing protein structures. Xwalk calculates and displays non-linear distances between chemically cross-linked amino acids on protein surfaces, while mimicking the flexibility and non-linearity of cross-linker molecules. It returns a **Solvent Accessible Surface (SAS)** distance, which corresponds to the length of
the shortest path between two amino acids, where the path leads through solvent
occupied space without penetrating the protein surface.


#TEST
Test examples to execute Xwalk can be found in the test subdirectory.


#CLASSPATH
You might want to consider adding the bin directory into the CLASSPATH
environment of your SHELL, which avoids the usage of the ```-cp``` flag when executing
Xwalk.


#COMMANDLINE PARAMETERS
A list of all commandline parameters can be retrieved by typing ```-help``` before
a Xwalk execution.  


#OUTPUT
Distance information will be printed out to the STDOUT channel or via ```-out``` to a
local file in the following tab delimeted format:

```
1	1brs.pdb	LYS-1-D-CB	LYS-2-D-CB	1	5.9	6.6	0.0449	0.0685	KKAVINGEQIR-KAVINGEQIR
```
Column | Information       |
-------|-------------------|
1      | Index     
2      | Filename     
3      | PDB information of the 1st amino acid in the format: PDBThreeLetterCode-ResidueId-ChainId-AtomName
4      | PDB information of the 2nd amino acid in the format: PDBThreeLetterCode-ResidueId-ChainId-AtomName
5      | Sequence distance within the PDB file as calculated from the residue identifiers of the 1st and 2nd amino acid.
6      | Euclidean distance between 1st and 2nd amino acid.
7      | SAS distance between  1st and 2nd amino acid. 
8      | Probability for observing the Euclidean distance in a XL experiment using DSS or BS3. 
9      | Probability for observing the SAS distance in a XL experiment using DSS or BS3. 
10     | Shortest tryptic peptide sequence for both amino acids.

Setting ```-pymol -out xxx.pml``` on the commandline will write a PyMOL script into
the file ```xxx.pml```, which can subsequently be loaded into the molecular viewer PyMOL to visualise
SAS distance paths (see [NOTES](#notes)).


#ALGORITHM
The SAS distance is calculated using a grid and the breadth-first search algorithm
to search within the grid for the shortest path between two points on the protein
surface using following algorithm:

1. Read in input data
   1. ```xyz.pdb``` spatial coordinates of a protein or protein complex in PDB
       format.
   2. ```maxdist``` maximum distance of the path (i.e. the length of the
       cross-linker + AA side chain length)
   3. ```listXL``` a list of experimentally determined cross-linked lysine residues.
2. Remove all non-protein atoms in ```xyz.pdb``` and assign protein atoms a van der Waals radius sum of SURFNET atom radii + solvent radius
<a name="step3"></a>
3. Select a random lysine pair *(K<sub>i</sub>,K<sub>j</sub>)* from ```listXL```, 
4. Check Euclidean distance *(Euc)* of *(K<sub>i</sub>,K<sub>j</sub>)*. Continue, if ```Euc > maxdist```, disregard otherwise and go back to [3.](#step3) 
5. Generate a grid of size maxdist and grid spacing 1 Angstroem centered at AAa
6. Set ```Integer.MAX_VALUE``` as distance for all grid cells and label grid cells as residing in the
  1. protein
  2. solvent
  3. boundary between protein and solvent.
7. Label grid cells residing in *(K<sub>i</sub>,K<sub>j</sub>)* as solvent
8. Set distance ```dist = 0.0``` for central grid cell of *K<sub>i</sub>* and store grid cell in the active list listactive
    <a name="step9"></a>
9. Start breadth-first search. Iterate through listactive 
  1. Check that grid cell *i* is labeled as solvent 
  2. Find all immediate neighbors listneighbour
  3. Iterate through listneighbour
    1. Check that grid cell *j* is labeled as solvent
    <a name="step9.3.2"></a>
    2. Compute new distance for grid cell *j* as the sum of the distance in grid cell *i* and the Euclidean distance between grid cell *i* and *j* 
    3. If distance sum in [9.3.2](#9.3.2) is smaller than the current distance in grid cell *j*, store the distance sum as new distance for grid cell *j* and add grid cell *j* to the new active list ```listnew_active```,
10. Go back to step [9.](#step9) with ```listactive = listnew_active```.



<a name="notes"></a>
#NOTES
- As the SAS distance is based on a grid calculation, the default heap size of
  the Java VM with 64MB is likely to be too small. You can increase the heap
  size with ```java -Xmx512m```
- You can obtain PyMOL [here](https://sourceforge.net/projects/pymol/) for free.
- Beware that in order for PyMOL to recognize the script file, the file must
  have the name ending ```.pml```.
- You can load the script directly at the startup of PyMOL, i.e. with the
  command 
  ```
  pymol 1brs.pml.
  ```


#CONTACT
<abdullah.kahraman@uzh.ch>


#LICENCE
Xwalk executable and libraries are available under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License via the Xwalk website.

* Anyone is free:
     * to copy, modify, distribute the software;

* Under the following conditions:
    * the original authors must be given credit by citing the original Xwalk paper: 
        * Kahraman A., Malmström L., Aebersold R. (2011). Xwalk: Computing and
        Visualizing Distances in Cross-linking Experiments. Bioinformatics,
        doi:10.1093/bioinformatics/btr348.
    * the software or derivate works must not be used for commercial purposes
    * derivate works must be licenced under the same or similar licence

* Any of the above conditions can be waived if the authors give permission.


VERSION HISTORY
---------------
### 0.6
* Added parameter ```-bfactor```, which adds the uncertainty of the atom coordinates
       as expressed by their B-factor/temperature factor to the maximum distance
       threshold.
* Distance probabilities were updated based on the latest update of XLdb with
       506 intra-protein and 62 inter-protein cross-links.
* Added ```Electro.java``` to the utils package for calculating electrostatic
       potentials on molecular surfaces. Requires [APBS, PDB2PQR](https://github.com/Electrostatics/apbs-pdb2pqr) and SIMS.
* Added ```Compactness.java``` to the utils package for calculating the number
       of atoms within a certain radius of any atom in a protein.
* Bug fixes.

###0.5
* Added parameter ```-xSC```, which is similar to ```-bb``` except that it leaves all side
       chain atoms of non-modified amino acid untouched, while removing only side
       chain atoms of cross-linked amino acids (except their CB atoms). In addition
       the solvent probe sphere is kept at the conventional radius of 1.4 Angstroem
       prior to SAS distance calculations. This parameter might be of value when
       side chain conformations of cross-linked residues are unknown.
* The premature halt of SAS distance calculation has been removed, which in the
       previous versions ensured a faster execution time, however, at the cost of
       the distance calculation accuracy.
* Distance probabilities have been recalculated using a larger data set of
       published cross-links.
* Bug fixes.

###0.4
* Maximum distance threshold for SAS distance calculation has been increased to 100
       Angstroem.
* Added option to read in ```.gzip``` files
* Warning messages are output for redundant XL entries in ```-dist``` files. 
* Added parameter ```-keepName``` to keep protein names from ```-dist``` files in the output.
* white background colour command in PyMOL scripts has been removed.

###0.3
* Improved performance and memory usage. Xwalk runs now 3-5 times
      faster, while using a third of the previous memory.

### 0.2.1
* Fixed a bug that caused Xwalk to crash on cases where cross-linked atoms
        were located outside the grid.
* Added a further number code to the SAS distance column for atom pairs that
        are located in closed cavities prohibiting a shortest path calculation.
* Bug fixes   

###0.2
* Probability calculation with option -prob
* Number code in SAS distance column (see ```-help```)
* Euclidean distance for all residue pairs in a distance file (```-dist```)
* PyMol option changed from ```-p``` to ```-pymol```
* Removed ```-grid```, ```-global``` and ```-all``` commandline options
* Bug fixes

###0.1
* First version of Xwalk
