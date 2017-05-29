// This function is executed outside the loop over all q's - the idea is to compute the geometry, scattering properties and other constraints here and pass the to the Model-function via the variable Constraints - keep in mind that the Model-function will be executed in parallel, if this has been enabled. Depending on the selected fitting routine, this function is also executed in parallel.

void ComputeConstraints(double * Parameters, double * VolumesOfMolecules, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /*
     * Parameters is an array of the parameters from the interface.
     * VoluesOfMolecules is an array the values from the sample information file.
     * ScatteringLengths is an array of the values from the sample information flie - computed so they match the contrast of the specific dataset.
     * Contrast is the contrast of the current dataset.
     * Constraints is the array, this function is meant to populate.
     * Concentration is value of the concentration read from the .card-file - rescaled from mM to 1/cm^3.
     * ProteinStructure contains all the information of the desired PDB-file in the format of the structure defined in Structs.h.
     * UserDefinedStructure is intended to be a "free" space to put whichever information one might need in order compute the scattering from a given model. The structure is defined in Structs.h.
     */
}


// This function recieves a q-value and should return the model intensity for this q-value - keep in mind that this function will be executed in parallel, if this has been enabled.

double Model(double q, double * Parameters, double * Constraints, double Contrast, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
    /*
     * q is the current value of the momentum transfer.
     * Parameters is an array of the parameters from the interface.
     * Constraints is the array populated by the ComputeConstraints-functions.
     * Contrast is the number provided in the .card-file describing the current dataset.
     * ProteinStructure contains all the information of the desired PDB-file in the format of the structure defined in Structs.h.
     * UserDefinedStructure is intended to be a "free" space to put whichever information one might need in order compute the scattering from a given model. The structure is defined in Structs.h.
     */

	double Intensity = 0.0;

	return Intensity;
}


// This function is executed after the fitting routine has run and is used to output any desired parameter to the screen or a given file.

void OutputData(double ChiSquare, double QMin, double QMax, struct Parameter * Parameters, int NumberOfParameters,
                struct Dataset * Data, int NumberOfSpectra, char cardfilename[128], struct Protein ProteinStructure, 
                struct UserDefined UserDefinedStructure, char SampleFilename[256])
{
    /*
     * ChiSquare is the final chisquare.
     * QMin and QMax are the fitting ranges passed by the GUI.
     * Parameters is an struct of the parameter values, their limits and names.
     * NumberOfParameters is the number of parameters passed from the interface.
     * Dataset is a struct containing all datasets along with the constraints computed for these datasets.
     * NumberOfSpectra is the number of datasets imported from the .card-file.
     * cardfilename is the location of the .card-file assigned by the graphical interface.
     * ProteinStructure is the struct describing the .pdb-file assigned by the graphical interface.
     * UserDefinedStructure is intended to be a "free" space to put whichever information one might need in order compute the scattering from a given model. The structure is defined in Structs.h.
     * SampleFilename is the location of the sample information file assigned by the graphical interface.
     */
}
