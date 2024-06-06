import sys
import os
# import XCLASS packages
from xclass import task_LineIdentification
path = "/home/parkerwise/Research/Cloud-C/spectral-analysis/A.spw29/"
LocalPath="/home/parkerwise/Research/xclass/xclass_pip_off/examples/LineID/"
os.chdir(path)


###########################################################################
# TO MODIFY BY THE USER


# define path and name of default molfit file
DefaultMolfitFile = "molecules.molfit"


# define path and name of obs. xml file
ObsXMLFileName = "A.spw29.0.xml"



# define upper limit of overestimation
MaxOverestimationHeight = 500.0


# define tolerance
Tolerance = 65.0


# define path and name of algorithm xml files
AlgorithmXMLFileSMF = LocalPath + "files/my_algorithm-settings.xml"
AlgorithmXMLFileOverAll = LocalPath + "files/my_algorithm-settings.xml"


# define lower limit for column density of core components
MinColumnDensityEmis = 0.0


# define lower limit for column density of foreground components
MinColumnDensityAbs = 0.0


# define source name
SourceName = ""

SelectedMolecules = []

# define list of so-called strong molecules
StrongMoleculeList = []


## define path and name of cluster file
clusterdef = ""
###########################################################################


# call LineID function
IdentifiedLines, JobDir = task_LineIdentification.LineIdentificationCore(
                        MaxOverestimationHeight = MaxOverestimationHeight,
                        SourceName = SourceName,
                        DefaultMolfitFile = DefaultMolfitFile,
                        Tolerance = Tolerance,
                        SelectedMolecules = SelectedMolecules,
                        StrongMoleculeList = StrongMoleculeList,
                        MinColumnDensityEmis = MinColumnDensityEmis,
                        MinColumnDensityAbs = MinColumnDensityAbs,
                        AlgorithmXMLFileSMF = AlgorithmXMLFileSMF,
                        AlgorithmXMLFileOverAll = AlgorithmXMLFileOverAll,
                        experimentalData = ObsXMLFileName,
                        clusterdef = clusterdef)
