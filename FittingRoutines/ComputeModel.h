int ComputeModel(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, double *ChiXX, int NumberOfSmearingFolds,
                 double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int TotalNumberOfDatapoints, int NumberOfFreeParameters){
    /// Declarations
    double ChiSquare;

    /// Comptutation
    printf("Computing model function for given parameters...\n");

    ChiSquare = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                 TotalNumberOfDatapoints, NumberOfFreeParameters);

    *ChiXX = ChiSquare;

    return 0;
}
