# Shock-Tube-Bifurcation-Prediction
Shock Tube Bifurcation Prediction

"MlappPredictBifurcation" - Mitchell D. Hageman October 2024
PURPOSE:
  *Determine liklihood and extent of shock bifurcation in a shock tube
PREREQUISITES
  * Experimental Pressure Trace
  * Experimental laser absorbtion trace with strong schlieren spikes
  * MlappPredictBifurcationHelp.pdf is needed in the PWD (present working directory) if you want the help button to work.
REFERENCES
Ref 1: THE INTERACTION OF A REFLECTED SHOCK WAVE WITH THE BOUNDARY LAYER IN A SHOCK TUBE
      Herman Mark (NACA TM 1418) 1958
Ref 2: Measurement of reflected-shock bifurcation over a wide range of gas composition and pressure
      E. L. Petersen · R. K. Hanson, Shock Waves (2006) 15:333–340DOI 10.1007/s00193-006-0032-3
Ref 3: Influence of Reflected Shock and Boundary‐Layer Interaction on Shock‐Tube Flows
      L. Davies; J. L. Wilson Phys. Fluids 12, I-37–I-43 (1969)
INPUTS:
  * VacuumFilePath - full file path and name of csv or excel file with your Vacuum sample data
      -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_VacuumData.csv'
      -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
  * ShockFilePath - full file path and name of csv or excel file with your sample data
      -example:  'C:\Users\mitchell.hageman\Desktop\Data\20240923_001_ShockData.csv'
      -assumes that any dark (zero) signal offset correction has already been applied to the recorded voltages.
      -Assumes the sample data file and vacuum data file are organized identically.
  * NumHeaderLines - number of header lines in csv file before data begins.
      -example: My data has two header lines.  Row 1 is data labels, and Row 2 is the offset voltage. So my data  begins in row 3.
  *timeColumn - column of csv where time trace is.  Mine is Column A so Timecolumn=1.
  *PressureColumn - column of csv where detector trace is.  Mine is column B so PressureColumn=2.
  *PitchColumn - column of csv where LAS reference (Pitch) detector trace is.  Mine is column C so PitchColumn=3.
  *CatchColumn - column of csv where LAS signal (Catch) detector trace is.  Mine is column D so CatchColumn=4.
  *D [m] - Diameter of the pressure transducer.
  **Note: The following would likely come from a normal shock equation solver, such as the FROzen chemistry SHock solver (FROSH)**
  *Mis [-] - Incident shock Mach number
  *Mrs [-] - Reflected shock Mach number
   *gamma [-] - specific heat ratio of the test gas.  NEED TO CHECK WHETHER TO USE STATE 1 or 2 FOR THESE CALCS
  *MW_Mix [kg/kmol] - Molecular weight of the test gas.
OUTPUTS:
  *PointData
      -BifurcationHeight [m]                - Value - See Fig1 & Eq(5) Ref 1 (BifurcationHeight=l)
      -Initial Reflected Shock Spike time, tI [s]   - Value - Determined from Pressure trace
      -Reflected Shock arrival time, tA [s] - Value - Determined from Schlieren spike
  * PointLables - Data lables for PointData - Character Matrix
  *"Write Plot Data to File" Button writes PlotData and Point Data to .csv in user-selected folder, with lables in the first row.
VERSION NUMBER:
  * 1.0: January 2025 - initial release, Mitchell D. Hageman
