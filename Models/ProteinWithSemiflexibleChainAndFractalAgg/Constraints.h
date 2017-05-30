enum CONSTRAINTS {
    SLDWater,
    Conc,
    X,
    Y,
    Z,
    PolymerChainExcessScatteringLenght,
    I0,
    RZero
};
enum PARAMETERS {
    BACKN100,
    SCALEN100,
    BACKN0,
    SCALEN0,
    BACKX,
    SCALEX,
    SCALECONC,
    HYDR,
    GLYCV
};
enum SCATTERINGDATA {
    SOLVENT,
    MOLECULE
};

void ComputeConstraints(double * Parameters, double * Volumes, double * ScatteringLengths, double Contrast,
                        double Concentration, double * Constraints, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure)
{
  // Declare dummy variables needed in function
  double ScatteringLengthDensityOfWater;
  double ScatteringLengthOfWater;
  double VolumeOfWater;

  // Volume of water
  VolumeOfWater = Volumes[0];

  // Obtain scattering length densities
  ScatteringLengthOfWater  = ScatteringLengths[0];

  // Derive scattering lengths
  ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;

  // Assign values of the constraints
  Constraints[SLDWater] = ScatteringLengthDensityOfWater;
  Constraints[Conc] = 1.0 ;//Concentration;

}
