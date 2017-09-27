/// Function for computing chisqaure for a given set of parameters
double ComputeChiSquare(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules,
                        struct Protein * Ensemble, int NumberOfProteins, double * ProteinWeights, struct UserDefined * UserDefinedStructure, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Declarations
    double ChiSquare = 0.0;

    int i;
    int j;
    int k;
    int p;

    double q;
    double DummyQ;
    double DifferenceInIntensity;
    double Stepsize;
    double Intensity;
    double DummyParameters[NumberOfParameters];
    double SigmaOfQ;
    double Weight;

    struct UserDefined UserDefinedCopy;
    struct Protein ProteinStructureCopy;
    struct Protein ProteinStructure;

    // Generate parameter array
    for (i = 0; i < NumberOfParameters; ++i) {
        DummyParameters[i] = Parameters[i].Value;
    }

    // Compute constraints and elliptic distribution for each dataset
    
    for (i = 0; i < NumberOfSpectra; ++i) {
        //Loop over proteins in ensemble here
        printf("Dataset has %d points\n", Data[i].NumberOfDatapoints);
        printf("Name of dataset: %s \n", Data[i].Filename);
        for(p = 0; p < NumberOfProteins; ++p){
            printf("P is % d \n", p);
            ProteinStructure = Ensemble[p];
            printf("Residues: %f Atoms: %f \n", ProteinStructure.NumberOfResidues, ProteinStructure.NumberOfAtoms);
            Weight = ProteinWeights[i];
            ComputeConstraints(DummyParameters, VolumesOfMolecules, Data[i].ScatteringLengths, Data[i].Contrast,
                           Data[i].Concentration, Data[i].Constraints, ProteinStructure, &*UserDefinedStructure);
            //printf("Constraints computed\n");
            #pragma omp parallel for schedule(dynamic) private(q, j, k, Intensity, Stepsize, DummyQ, DifferenceInIntensity, ProteinStructureCopy, UserDefinedCopy, SigmaOfQ)
            for (j = Data[i].NMin; j < Data[i].NMax; ++j) {
                q = Data[i].QValues[j];

                InitializeUserDefinedStructure(&UserDefinedCopy);
                CopyUserDefined(&*UserDefinedStructure, &UserDefinedCopy);
                //This block should maybe be before the q loop?   
    		  	if (ProteinStructure.NumberOfAtoms != 0) {
                    ProteinStructureCopy.NumberOfAtoms    = ProteinStructure.NumberOfAtoms;
                    ProteinStructureCopy.NumberOfResidues = ProteinStructure.NumberOfResidues;
                    sprintf(ProteinStructureCopy.PDBFileLocation, "%s", ProteinStructure.PDBFileLocation);
                    //printf("The error is here:\n");
                    AllocateProteinStructure(&ProteinStructureCopy, ProteinStructureCopy.NumberOfResidues, ProteinStructureCopy.NumberOfAtoms);
                    //printf("Allocation done");
                    for (k = 0; k < ProteinStructureCopy.NumberOfResidues; ++k) {
                        ProteinStructureCopy.Residues[k].xVolume = ProteinStructure.Residues[k].xVolume;
                        ProteinStructureCopy.Residues[k].yVolume = ProteinStructure.Residues[k].yVolume;
                        ProteinStructureCopy.Residues[k].zVolume = ProteinStructure.Residues[k].zVolume;

                        ProteinStructureCopy.Residues[k].xXRayScattering = ProteinStructure.Residues[k].xXRayScattering;
                        ProteinStructureCopy.Residues[k].yXRayScattering = ProteinStructure.Residues[k].yXRayScattering;
                        ProteinStructureCopy.Residues[k].zXRayScattering = ProteinStructure.Residues[k].zXRayScattering;

                        ProteinStructureCopy.Residues[k].xNeutronScattering = ProteinStructure.Residues[k].xNeutronScattering;
                        ProteinStructureCopy.Residues[k].yNeutronScattering = ProteinStructure.Residues[k].yNeutronScattering;
                        ProteinStructureCopy.Residues[k].zNeutronScattering = ProteinStructure.Residues[k].zNeutronScattering;

                        ProteinStructureCopy.Residues[k].XRayScatteringLength    = ProteinStructure.Residues[k].XRayScatteringLength;
                        ProteinStructureCopy.Residues[k].NeutronScatteringLength = ProteinStructure.Residues[k].NeutronScatteringLength;
                        ProteinStructureCopy.Residues[k].Volume                  = ProteinStructure.Residues[k].Volume;
                        ProteinStructureCopy.Residues[k].ResidueID               = ProteinStructure.Residues[k].ResidueID;

                        strcpy(ProteinStructureCopy.Residues[k].Name, ProteinStructure.Residues[k].Name);
                    }

                    for (k = 0; k < ProteinStructureCopy.NumberOfAtoms; ++k) {
                        ProteinStructureCopy.Atoms[k].x = ProteinStructure.Atoms[k].x;
                        ProteinStructureCopy.Atoms[k].y = ProteinStructure.Atoms[k].y;
                        ProteinStructureCopy.Atoms[k].z = ProteinStructure.Atoms[k].z;

                        ProteinStructureCopy.Atoms[k].XRayScatteringLength    = ProteinStructure.Atoms[k].XRayScatteringLength;
                        ProteinStructureCopy.Atoms[k].NeutronScatteringLength = ProteinStructure.Atoms[k].NeutronScatteringLength;
                        ProteinStructureCopy.Atoms[k].Volume                  = ProteinStructure.Atoms[k].Volume;
                        ProteinStructureCopy.Atoms[k].Name                    = ProteinStructure.Atoms[k].Name;
                    }
    		  	}
                if (Data[i].IncludeResolutionEffects == true) {
                    Intensity = 0.0;
                    SigmaOfQ  = Data[i].SigmaQValues[j];
                    Stepsize  = 6.0 * SigmaOfQ / (1.0 * NumberOfSmearingFolds);

                    for (k = 0; k < NumberOfSmearingFolds; ++k) {
                        DummyQ = q + (k + 0.5 - NumberOfSmearingFolds / 2.0) * Stepsize;
    
                        if (DummyQ < 1e-5) {
                            DummyQ = 1e-5;
                        }
    
                        Intensity += Weight*Model(DummyQ, DummyParameters, Data[i].Constraints, Data[i].Contrast, ProteinStructureCopy, &UserDefinedCopy) * Data[i].ResolutionWeights[j][k];
                    }
                } else {
                    Intensity += Weight*Model(q, DummyParameters, Data[i].Constraints, Data[i].Contrast, ProteinStructureCopy, &UserDefinedCopy);
                }
                Data[i].FitValues[j] = Intensity;
                DifferenceInIntensity = Data[i].IValues[j] - Intensity;
    
                #pragma omp atomic
                ChiSquare += Weight*(DifferenceInIntensity * DifferenceInIntensity * Data[i].SigmaValues[j] / (1.0 * (TotalNumberOfDatapoints - NumberOfFreeParameters)));
    
                FreeUserDefined(&UserDefinedCopy);
    
                if (ProteinStructure.NumberOfAtoms != 0) {
                    free(ProteinStructureCopy.Residues);
                    free(ProteinStructureCopy.Atoms);
                }
            }
        }    
    }

    return ChiSquare;
}

/// Compute the gradient for a given point
void ComputeGradient(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds, double * VolumesOfMolecules,
                    struct Protein * Ensemble, int NumberOfProteins, double * ProteinWeights, struct UserDefined * UserDefinedStructure, double * Gradient,
                    double DeltaForDifferentiations, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Declarations
    int i;
    double dParameter;
    double InitialValueOfParameter;
    const double Delta = DeltaForDifferentiations;

    double ChiSquarePlus;
    double ChiSquareMinus;

    // Computation
    for (i = 0; i < NumberOfParameters; ++i) {
        Gradient[i] = 0;

        if (Parameters[i].iParameter == true) {
            InitialValueOfParameter = Parameters[i].Value;

            if (fabs(Parameters[i].Value) < 1e-10) {
                Parameters[i].Value = 1e-10;
            }

            dParameter = Parameters[i].Value * Delta;

            Parameters[i].Value += dParameter;

            ChiSquarePlus  = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                            ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

            Parameters[i].Value = InitialValueOfParameter;
            Parameters[i].Value -= dParameter;

            ChiSquareMinus = ComputeChiSquare(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                            ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

            Parameters[i].Value = InitialValueOfParameter;
            Gradient[i] = (ChiSquarePlus - ChiSquareMinus) / (2.0 * dParameter);
        }
    }
}

/// Perform linesearch (lnsrch from Numerical Recipes)
double PerformLinesearch(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, int NumberOfSmearingFolds,
                         double * VolumesOfMolecules, struct Protein * Ensemble, int NumberOfProteins, double * ProteinWeights, struct UserDefined * UserDefinedStructure, double OriginalValue,
                         double * Gradient, double * Direction, struct Parameter * NewParameters, double MaxStepSize, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Convergence declarations
    const double Alpha = 1e-4;
    const double ToleranceLevel = 1e-15;

    // Dummy declarations
    int i;
    double Dummy1 = 0.0;
    double Dummy2;
    double CurrentValue;

    // Supporting declarations
    double Lambda;
    double LambdaMin;
    double PreviousLambda = 0.0;
    double TemporaryLambda;
    double PreviousValue = 0.0;

    double AParameter;
    double BParameter;
    double Discriminant;
    double RHS1;
    double RHS2;

    double Sum = 0.0;
    double Slope = 0.0;

    // Computation
    for (i = 0; i < NumberOfParameters; ++i) {
        Sum += pow(Direction[i], 2);
    }

    Sum = sqrt(Sum);

    // Rescale, if step is too big
    if (Sum > MaxStepSize) {

        for (i = 0; i < NumberOfParameters; ++i) {
            Direction[i] *= MaxStepSize / Sum;
        }
    }

    // Compute slope in the given direction
    for (i = 0; i < NumberOfParameters; ++i) {
        Slope += Gradient[i] * Direction[i];
    }

    if (Slope >= 0.0) {
        // New value is to close to previous value
        printf("Roundoff problem occured, when performing linesearch. \n");
        return OriginalValue;
    }

    // Compute min lambda
    for (i = 0; i < NumberOfParameters; ++i) {

        if (fabs(Parameters[i].Value) > 1.0) {
            Dummy2 = fabs(Direction[i]) / fabs(Parameters[i].Value);

        } else {
            Dummy2 = fabs(Direction[i]);
        }

        if (Dummy2 > Dummy1) {
            Dummy1 = Dummy2;
        }
    }

    LambdaMin = ToleranceLevel / Dummy1;
    Lambda = 1.0;

    // Initialize (infinite) loop
    while(true){

        for (i = 0; i < NumberOfParameters; ++i) {
            NewParameters[i].Value = Parameters[i].Value + Lambda * Direction[i];
        }

        CurrentValue = ComputeChiSquare(Data, NumberOfSpectra, NewParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, Ensemble, NumberOfProteins,
                                    ProteinWeights, &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

        if (Lambda < LambdaMin) {

            for (i = 0; i < NumberOfParameters; ++i) {
                NewParameters[i].Value = Parameters[i].Value;
            }

            // No minimum found - convergence reached
            return OriginalValue;

        } else if (CurrentValue <= OriginalValue + Alpha * Lambda * Slope) {

            // Sufficient function decrease found - conclude search
            return CurrentValue;

        } else {
            // Rescale lambda and try again
            if (Lambda == 1.0) {
                TemporaryLambda = - Slope / (2.0 * (CurrentValue - OriginalValue - Slope));

            } else {
                RHS1 = CurrentValue - OriginalValue - Lambda * Slope;
                RHS2 = PreviousValue - OriginalValue - PreviousLambda * Slope;

                AParameter = (RHS1 / pow(Lambda, 2) - RHS2 / pow(PreviousLambda, 2)) / (Lambda - PreviousLambda);
                BParameter = (- PreviousLambda * RHS1 / pow(Lambda, 2) + Lambda * RHS2 / pow(PreviousLambda, 2)) / (Lambda - PreviousLambda);

                if (AParameter == 0.0) {
                    TemporaryLambda = - Slope / (2.0 * BParameter);

                } else {
                    Discriminant = pow(BParameter, 2) - 3.0 * AParameter * Slope;

                    if (Discriminant < 0.0) {
                        TemporaryLambda = 0.5 * Lambda;

                    } else if (BParameter <= 0.0) {
                        TemporaryLambda = (- BParameter + sqrt(Discriminant)) / (3.0 * AParameter);

                    } else {
                        TemporaryLambda = - Slope / (BParameter + sqrt(Discriminant));
                    }
                }

                if (TemporaryLambda > 0.5 * Lambda) {
                    TemporaryLambda = 0.5 * Lambda;
                }
            }
        }

        PreviousLambda = Lambda;
        PreviousValue = CurrentValue;

        if (TemporaryLambda > 0.1 * Lambda) {
            Lambda = TemporaryLambda;
        } else {
            Lambda = 0.1 * Lambda;
        }
    }
}
