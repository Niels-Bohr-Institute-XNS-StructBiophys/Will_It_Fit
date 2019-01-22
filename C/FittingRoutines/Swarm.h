int Swarm(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, double * Chisquare, int NumberOfSmearingFolds,
          double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int NumberOfAgents, int NumberOfFlights, int MaximumNumberOfIterations,
          bool ProfileLikelihood, int HighestNumberOfDatapoints, int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Declaration of dummy variables
    int i;
    int j;
    int k;
    int l;
    int m;
    int n;

    // Declarations of supporting variables
    double MaxPosVelocity[NumberOfParameters];
    double MaxNegVelocity[NumberOfParameters];
    struct Parameter * DummyParameters;
    struct Dataset * DummyData;
    struct UserDefined UserDefinedCopy;
    struct Protein ProteinStructureCopy;
    double VelocityRange;
    MaximumNumberOfIterations = 30;

    const double CognitiveAccelerationCoefficient = 2.05;
    const double SocialAccelerationCoefficient = 2.05;
    const double CombinedAccelerationCoefficient = CognitiveAccelerationCoefficient + SocialAccelerationCoefficient;
    const double K = 1.0;
    const double ConstrictionFactor = 2.0 * K / (fabs(2.0 - CombinedAccelerationCoefficient -
                                                      sqrt(pow(CombinedAccelerationCoefficient, 2) - 4.0 * CombinedAccelerationCoefficient)));

    // Properties of agents
    struct Agent Population[NumberOfAgents];
    struct Agent BestCandidate;

    // Seed RNG
    srand(time(0));

    // Initialize agents
    for (i = 0; i < NumberOfAgents; ++i) {
        InitializeAgent(&Population[i], NumberOfParameters);
    }

    InitializeAgent(&BestCandidate, NumberOfParameters);

    for (i = 0; i < NumberOfParameters; ++i) {
        BestCandidate.Parameters[i] = Parameters[i].Value;
    }

    BestCandidate.Chisquare = ComputeChiSquareSwarm(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure,
                                                    &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Initialize main loop
    i = 0;

    while (i < NumberOfFlights) {
        // Update flight number
        ++i;

        // Initialize agents
        for (j = 0; j < NumberOfAgents; ++j) {

            for (k = 0; k < NumberOfParameters; ++k) {

                if (Parameters[k].iParameter == true) {
                    Population[j].Parameters[k] = BestCandidate.Parameters[k];

                    MaxPosVelocity[k] = 0.05 * (Parameters[k].MaxValue - Parameters[k].Value);
                    MaxNegVelocity[k] = 0.05 * (Parameters[k].MinValue - Parameters[k].Value);
                    VelocityRange = fabs(MaxPosVelocity[k] - MaxNegVelocity[k]);

                    Population[j].Velocity[k] = MaxNegVelocity[k] + RandomFraction() * VelocityRange;

                    Population[j].PersonalBestParameters[k] = Parameters[k].Value;
                    Population[j].PersonalBestChisquare = BestCandidate.Chisquare;
                } else {
                    Population[j].Parameters[k] = Parameters[k].Value;
                    MaxPosVelocity[k] = 0.0;
                    MaxNegVelocity[k] = 0.0;
                    Population[j].Velocity[k] = 0.0;
                }
            }
        }

        // Initialize i'th flight
        j = 0;

        while (j < MaximumNumberOfIterations) {
            // Update iteration
            ++j;

            // Evaluate initial conditions
            if (j != 1) {

                #pragma omp parallel for schedule(dynamic) private(k, l, m, DummyParameters, DummyData, ProteinStructureCopy, UserDefinedCopy)
                for (k = 0; k < NumberOfAgents; ++k) {
                    // Initialization
                    AllocateParameters(&DummyParameters, NumberOfParameters);
                    AllocateData(&DummyData, NumberOfSpectra);

                    InitializeUserDefinedStructure(&UserDefinedCopy);
                    CopyUserDefined(&*UserDefinedStructure, &UserDefinedCopy);

                    for (l = 0; l < NumberOfParameters; ++l) {
                        DummyParameters[l].Value = Population[k].Parameters[l];
                    }

                    for (l = 0; l < NumberOfSpectra; ++l) {
                        Initialize1DArray(&DummyData[l].QValues,           HighestNumberOfDatapoints);
                        Initialize1DArray(&DummyData[l].IValues,           HighestNumberOfDatapoints);
                        Initialize1DArray(&DummyData[l].FitValues,         HighestNumberOfDatapoints);
                        Initialize1DArray(&DummyData[l].SigmaValues,       HighestNumberOfDatapoints);
                        Initialize1DArray(&DummyData[l].SigmaQValues,      HighestNumberOfDatapoints);
                        Initialize1DArray(&DummyData[l].Constraints,       MaxNumberOfConstraints);
                        Initialize1DArray(&DummyData[l].ScatteringLengths, NumberOfSampleInformations);
                        Initialize2DArray(&DummyData[l].ResolutionWeights, HighestNumberOfDatapoints, NumberOfSmearingFolds);

                        for (m = 0; m < HighestNumberOfDatapoints; ++m) {
                            DummyData[l].QValues[m]      = Data[l].QValues[m];
                            DummyData[l].IValues[m]      = Data[l].IValues[m];
                            DummyData[l].SigmaValues[m]  = Data[l].SigmaValues[m];
                            DummyData[l].SigmaQValues[m] = Data[l].SigmaQValues[m];

                            for (n = 0; n < NumberOfSmearingFolds; ++n) {
                                DummyData[l].ResolutionWeights[m][n] = Data[l].ResolutionWeights[m][n];
                            }
                        }

                        for (m = 0; m < NumberOfSampleInformations; ++m) {
                            DummyData[l].ScatteringLengths[m] = Data[l].ScatteringLengths[m];
                        }

                        DummyData[l].IncludeResolutionEffects = Data[l].IncludeResolutionEffects;
                        DummyData[l].Concentration            = Data[l].Concentration;
                        DummyData[l].Contrast                 = Data[l].Contrast;
                        DummyData[l].NMax                     = Data[l].NMax;
                        DummyData[l].NMin                     = Data[l].NMin;
                    }

			// Initialization and full copy of ProteinStructure (ProteinStructureCopy)
			if ( ProteinStructure.NumberOfAtoms != 0 )
			{
				ProteinStructureCopy.NumberOfAtoms    = ProteinStructure.NumberOfAtoms ;
				ProteinStructureCopy.NumberOfResidues = ProteinStructure.NumberOfResidues ;

				AllocateProteinStructure( &ProteinStructureCopy, ProteinStructureCopy.NumberOfResidues, ProteinStructureCopy.NumberOfAtoms) ;

				CopyProteinStructure( &ProteinStructureCopy, &ProteinStructure) ;
			}
/*
                    if (ProteinStructure.NumberOfAtoms != 0) {
                        ProteinStructureCopy.NumberOfAtoms    = ProteinStructure.NumberOfAtoms;
                        ProteinStructureCopy.NumberOfResidues = ProteinStructure.NumberOfResidues;
                        sprintf(ProteinStructureCopy.PDBFileLocation, "%s", ProteinStructure.PDBFileLocation);

                        AllocateProteinStructure(&ProteinStructureCopy, ProteinStructureCopy.NumberOfResidues, ProteinStructureCopy.NumberOfAtoms);

                        for (l = 0; l < ProteinStructureCopy.NumberOfResidues; ++l) {
                            ProteinStructureCopy.Residues[l].xVolume = ProteinStructure.Residues[l].xVolume;
                            ProteinStructureCopy.Residues[l].yVolume = ProteinStructure.Residues[l].yVolume;
                            ProteinStructureCopy.Residues[l].zVolume = ProteinStructure.Residues[l].zVolume;

                            ProteinStructureCopy.Residues[l].xXRayScattering = ProteinStructure.Residues[l].xXRayScattering;
                            ProteinStructureCopy.Residues[l].yXRayScattering = ProteinStructure.Residues[l].yXRayScattering;
                            ProteinStructureCopy.Residues[l].zXRayScattering = ProteinStructure.Residues[l].zXRayScattering;

                            ProteinStructureCopy.Residues[l].xNeutronScattering = ProteinStructure.Residues[l].xNeutronScattering;
                            ProteinStructureCopy.Residues[l].yNeutronScattering = ProteinStructure.Residues[l].yNeutronScattering;
                            ProteinStructureCopy.Residues[l].zNeutronScattering = ProteinStructure.Residues[l].zNeutronScattering;

                            ProteinStructureCopy.Residues[l].XRayScatteringLength    = ProteinStructure.Residues[l].XRayScatteringLength;
                            ProteinStructureCopy.Residues[l].NeutronScatteringLength = ProteinStructure.Residues[l].NeutronScatteringLength;
                            ProteinStructureCopy.Residues[l].Volume                  = ProteinStructure.Residues[l].Volume;
                            ProteinStructureCopy.Residues[l].ResidueID               = ProteinStructure.Residues[l].ResidueID;

                            strcpy(ProteinStructureCopy.Residues[l].Name, ProteinStructure.Residues[l].Name);
                        }

                        for (l = 0; l < ProteinStructureCopy.NumberOfAtoms; ++l) {
                            ProteinStructureCopy.Atoms[l].x = ProteinStructure.Atoms[l].x;
                            ProteinStructureCopy.Atoms[l].y = ProteinStructure.Atoms[l].y;
                            ProteinStructureCopy.Atoms[l].z = ProteinStructure.Atoms[l].z;

                            ProteinStructureCopy.Atoms[l].XRayScatteringLength    = ProteinStructure.Atoms[l].XRayScatteringLength;
                            ProteinStructureCopy.Atoms[l].NeutronScatteringLength = ProteinStructure.Atoms[l].NeutronScatteringLength;
                            ProteinStructureCopy.Atoms[l].Volume                  = ProteinStructure.Atoms[l].Volume;
                            ProteinStructureCopy.Atoms[l].Name                    = ProteinStructure.Atoms[l].Name;
                        }
                    }
*/


                    // Computation
                    Population[k].Chisquare = ComputeChiSquareSwarm(DummyData, NumberOfSpectra, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules,
                                                                    ProteinStructureCopy, &UserDefinedCopy, TotalNumberOfDatapoints, NumberOfFreeParameters);

                    // Free memory
                    free(DummyParameters);

                    for (l = 0; l < NumberOfSpectra; ++l) {
                        free(DummyData[l].QValues);
                        free(DummyData[l].IValues);
                        free(DummyData[l].FitValues);
                        free(DummyData[l].SigmaValues);
                        free(DummyData[l].SigmaQValues);
                        free(DummyData[l].Constraints);
                        free(DummyData[l].ScatteringLengths);

                        for (m = 0; m < HighestNumberOfDatapoints; ++m) {
                            free(DummyData[l].ResolutionWeights[m]);
                        }

                        free(DummyData[l].ResolutionWeights);
                    }

                    FreeUserDefined(&UserDefinedCopy);
                    free(DummyData);

                    if (ProteinStructure.NumberOfAtoms != 0) {
                        free(ProteinStructureCopy.Residues);
                        free(ProteinStructureCopy.Atoms);
                    }
                }

                for (k = 0; k < NumberOfAgents; ++k) {

                    // Check for personal improvement
                    if (Population[k].PersonalBestChisquare > Population[k].Chisquare) {
                        Population[k].PersonalBestChisquare = Population[k].Chisquare;

                        for (l = 0; l < NumberOfParameters; ++l) {
                            Population[k].PersonalBestParameters[l] = Population[k].Parameters[l];
                        }
                    }

                    // Check for global improvement
                    if (BestCandidate.Chisquare > Population[k].Chisquare) {
                        BestCandidate.Chisquare = Population[k].Chisquare;

                        for (l = 0; l < NumberOfParameters; ++l) {
                            BestCandidate.Parameters[l] = Population[k].Parameters[l];
                        }
                    }
                }
            }

            if (j != MaximumNumberOfIterations) {

                // Update velocities
                for (k = 0; k < NumberOfAgents; ++k) {

                    for (l = 0; l < NumberOfParameters; ++l) {

                        if (Parameters[l].iParameter == true) {
                            Population[k].Velocity[l] = ConstrictionFactor * (Population[k].Velocity[l] +
                                                        CognitiveAccelerationCoefficient * RandomFraction() *
                                                                              (Population[k].PersonalBestParameters[l] - Population[k].Parameters[l]) +
                                                        SocialAccelerationCoefficient * RandomFraction() *
                                                                              (BestCandidate.Parameters[l] - Population[k].Parameters[l]));

                            if (Population[k].Velocity[l] > MaxPosVelocity[l]) {
                                Population[k].Velocity[l] = MaxPosVelocity[l];
                            } else if (Population[k].Velocity[l] < MaxNegVelocity[l]) {
                                Population[k].Velocity[l] = MaxNegVelocity[l];
                            }
                        }
                    }
                }

                // Update coordinates
                for (k = 0; k < NumberOfAgents; ++k) {

                    for (l = 0; l < NumberOfParameters; ++l) {

                        if (Parameters[l].iParameter == true) {
                            Population[k].Parameters[l] += Population[k].Velocity[l];

                            if (Population[k].Parameters[l] > Parameters[l].MaxValue) {
                                Population[k].Parameters[l] = Parameters[l].MaxValue;
                            } else if (Population[k].Parameters[l] < Parameters[l].MinValue) {
                                Population[k].Parameters[l] = Parameters[l].MinValue;
                            }
                        }
                    }
                }
            }

            // Print best candidate
            if (j % 5 == 0 && ProfileLikelihood == false) {
                printf("Best solution - after %d iterations in flight number %d: \n", j, i);
                printf("Chisquare = %g \n", BestCandidate.Chisquare);
                printf("\n");
                printf("Parameters: \n");

                for (k = 0; k < NumberOfParameters; ++k) {
                    printf("%-30s = %g \n", Parameters[k].Name, BestCandidate.Parameters[k]);
                }

                printf("\n");
                ClearScreen( stdout );
            }
        }
    }

    // Conclusion
    for (i = 0; i < NumberOfParameters; ++i) {
        Parameters[i].Value = BestCandidate.Parameters[i];
    }

    *Chisquare = ComputeChiSquareSwarm(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                       TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Free agents from memory
    free(BestCandidate.Parameters);
    free(BestCandidate.Velocity);
    free(BestCandidate.PersonalBestParameters);

    for (i = 0; i < NumberOfAgents; ++i) {
        free(Population[i].Parameters);
        free(Population[i].Velocity);
        free(Population[i].PersonalBestParameters);
    }

    return 0;
}
