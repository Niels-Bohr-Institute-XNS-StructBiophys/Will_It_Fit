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
    ChainX,
    ChainY,
    ChainZ,
    ChainDr,
    ContourLength,
    KuhnLength,
    ChainContrast,
    BACKN100,
    SCALEN100,
    BACKN0,
    SCALEN0,
    BACKX,
    SCALEX,
    SCALECONC,
    PROTSCALE,
    FracDim,
    AggRg,
    RZeroPar,
    Xagg,
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
  double ChainDistanceFromSurface;
  double alphaSQ;
  double RgSQ ;

  double Sagg;
  double Dm = Parameters[FracDim];
  double RgAgg = Parameters[AggRg];
  double r0 = Parameters[RZeroPar];
  double XiAgg = sqrt(2 * pow(RgAgg, 2) / (Dm * (Dm + 1)));
  // Volume of water
  VolumeOfWater = Volumes[0];

  // Obtain scattering length densities
  ScatteringLengthOfWater  = ScatteringLengths[0];

  // Derive scattering lengths
  ScatteringLengthDensityOfWater = ScatteringLengthOfWater / VolumeOfWater;

  // Assign values of the constraints
  Constraints[SLDWater] = ScatteringLengthDensityOfWater;
  Constraints[Conc] = Concentration;
  Constraints[X] = Parameters[ChainDr]*Parameters[ChainX];
  Constraints[Y] = Parameters[ChainDr]*Parameters[ChainY];
  Constraints[Z] = Parameters[ChainDr]*Parameters[ChainZ];
  Constraints[PolymerChainExcessScatteringLenght] = Parameters[ChainContrast]*ScatteringLengths[1];
  Constraints[I0] = 0.0;
  alphaSQ =pow(1
    +pow((Parameters[ContourLength]/Parameters[KuhnLength]/3.12),2)
    +pow((Parameters[ContourLength]/Parameters[KuhnLength]/8.67),3),(0.17/3));
  RgSQ = (Parameters[ContourLength]*Parameters[KuhnLength]/6
    -pow(Parameters[KuhnLength],2)/4
    +Parameters[KuhnLength]/Parameters[ContourLength]/4
    -1/8/pow(Parameters[ContourLength],2)*pow(Parameters[KuhnLength],2)*(1-exp(-2*Parameters[ContourLength]/Parameters[KuhnLength]))
  ) ;
  Constraints[RZero] = pow(alphaSQ*RgSQ,0.5);
//  printf("Aggregate aggregation number: %f\n",tgamma(Dm+1)*pow(XiAgg/r0,Dm));

//  printf("Alpha^2 %f\n", alphaSQ);
//  printf("Rg of aggregate subunit %f\n", Constraints[RZero]);
}
