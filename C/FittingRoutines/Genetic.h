int Genetic(struct Dataset * Data, int NumberOfSpectra, struct Parameter * Parameters, int NumberOfParameters, double *Chisquare, int NumberOfSmearingFolds,
            double * VolumesOfMolecules, struct Protein ProteinStructure, struct UserDefined * UserDefinedStructure, int SizeOfPopulation, int MaxNumberOfGenerations,
            int HighestNumberOfDatapoints, int NumberOfSampleInformations, int TotalNumberOfDatapoints, int NumberOfFreeParameters)
{
    // Declaration of dummy variables
    int i;
    int j;
    int k;
    int l;

    // Declare populations
    struct Individual Population[SizeOfPopulation];
    struct Individual NewPopulation[SizeOfPopulation];
    struct Individual ParentOne;
    struct Individual ParentTwo;
    struct Individual BestCandidate;

    double FitnessRanking[SizeOfPopulation];
    double ProbabilityOfMutation;

    // More parameters
    int Generation = 0;
    double Criteria;
    double SumOfFitness;
    bool Stop = false;
    const bool DoExcessivePrinting = false;
    struct Parameter * DummyParameters;
    struct Dataset * DummyData;
    struct UserDefined UserDefinedCopy;
    struct Protein ProteinStructureCopy;
    FILE *Outfile;

    // Allocate memonry for subjects
    AllocateSubject(&ParentOne, NumberOfParameters);
    AllocateSubject(&ParentTwo, NumberOfParameters);
    AllocateSubject(&BestCandidate, NumberOfParameters);

    for (i = 0; i < SizeOfPopulation; ++i) {
        AllocateSubject(&Population[i], NumberOfParameters);
        AllocateSubject(&NewPopulation[i], NumberOfParameters);
    }

    // Seed RNG
    srand(time(0));

    // Populate
    for (i = 0; i < SizeOfPopulation; ++i) {
        InitializeSubject(&Population[i], Parameters, NumberOfParameters);
    }

    // Choose arbitrary subject as initial best candidate
    for (i = 0; i < NumberOfParameters; ++i) {
        Population[0].Parameters[i] = Parameters[i].Value;
    }

    Population[0].Chisquare = ComputeChiSquareGenetic(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure,
                                                      &*UserDefinedStructure, TotalNumberOfDatapoints, NumberOfFreeParameters);

    BestCandidate.Chisquare = Population[0].Chisquare;

    for (i = 0; i < NumberOfParameters; ++i) {
        BestCandidate.Parameters[i] = Population[0].Parameters[i];
    }

    // Do overly excessive printing
    if (DoExcessivePrinting == true) {
        Outfile = fopen("Subjects.mcp", "w+");

        fprintf(Outfile, "Chisquare         ");

        for (i = 0; i < NumberOfParameters; ++i) {
            fprintf(Outfile, "%15s   ", Parameters[i].Name);
        }

        fprintf(Outfile, "\n");
    }

    // Begin main loop
    while (Stop == false) {

        // Compute chisquare
        #pragma omp parallel for schedule(dynamic) private(i, j, k, DummyParameters, DummyData, ProteinStructureCopy, UserDefinedCopy)
        for (i = 0; i < SizeOfPopulation; ++i) {
            AllocateParameters(&DummyParameters, NumberOfParameters);
            AllocateData(&DummyData, NumberOfSpectra);

            InitializeUserDefinedStructure(&UserDefinedCopy);
            CopyUserDefined(&*UserDefinedStructure, &UserDefinedCopy);

            for (j = 0; j < NumberOfParameters; ++j) {
                DummyParameters[j].Value = Population[i].Parameters[j];
            }

            for (j = 0; j < NumberOfSpectra; ++j) {
                Initialize1DArray(&DummyData[j].QValues,           HighestNumberOfDatapoints);
                Initialize1DArray(&DummyData[j].IValues,           HighestNumberOfDatapoints);
                Initialize1DArray(&DummyData[j].FitValues,         HighestNumberOfDatapoints);
                Initialize1DArray(&DummyData[j].SigmaValues,       HighestNumberOfDatapoints);
                Initialize1DArray(&DummyData[j].SigmaQValues,      HighestNumberOfDatapoints);
                Initialize1DArray(&DummyData[j].Constraints,       MaxNumberOfConstraints);
                Initialize1DArray(&DummyData[j].ScatteringLengths, NumberOfSampleInformations);
                Initialize2DArray(&DummyData[j].ResolutionWeights, HighestNumberOfDatapoints, NumberOfSmearingFolds);


                for (k = 0; k < HighestNumberOfDatapoints; ++k) {
                    DummyData[j].QValues[k]      = Data[j].QValues[k];
                    DummyData[j].IValues[k]      = Data[j].IValues[k];
                    DummyData[j].SigmaValues[k]  = Data[j].SigmaValues[k];
                    DummyData[j].SigmaQValues[k] = Data[j].SigmaQValues[k];

                    for (l = 0; l < NumberOfSmearingFolds; ++l) {
                        DummyData[j].ResolutionWeights[k][l] = Data[j].ResolutionWeights[k][l];
                    }
                }

                for (k = 0; k < NumberOfSampleInformations; ++k) {
                    DummyData[j].ScatteringLengths[k] = Data[j].ScatteringLengths[k];
                }

                DummyData[j].IncludeResolutionEffects = Data[j].IncludeResolutionEffects;
                DummyData[j].Concentration            = Data[j].Concentration;
                DummyData[j].Contrast                 = Data[j].Contrast;
                DummyData[j].NMax                     = Data[j].NMax;
                DummyData[j].NMin                     = Data[j].NMin;
            }

            if (ProteinStructure.NumberOfAtoms != 0) {
                ProteinStructureCopy.NumberOfAtoms    = ProteinStructure.NumberOfAtoms;
                ProteinStructureCopy.NumberOfResidues = ProteinStructure.NumberOfResidues;
                sprintf(ProteinStructureCopy.PDBFileLocation, "%s", ProteinStructure.PDBFileLocation);

                AllocateProteinStructure(&ProteinStructureCopy, ProteinStructureCopy.NumberOfResidues, ProteinStructureCopy.NumberOfAtoms);

                for (j = 0; j < ProteinStructureCopy.NumberOfResidues; ++j) {
                    ProteinStructureCopy.Residues[j].xVolume = ProteinStructure.Residues[j].xVolume;
                    ProteinStructureCopy.Residues[j].yVolume = ProteinStructure.Residues[j].yVolume;
                    ProteinStructureCopy.Residues[j].zVolume = ProteinStructure.Residues[j].zVolume;

                    ProteinStructureCopy.Residues[j].xXRayScattering = ProteinStructure.Residues[j].xXRayScattering;
                    ProteinStructureCopy.Residues[j].yXRayScattering = ProteinStructure.Residues[j].yXRayScattering;
                    ProteinStructureCopy.Residues[j].zXRayScattering = ProteinStructure.Residues[j].zXRayScattering;

                    ProteinStructureCopy.Residues[j].xNeutronScattering = ProteinStructure.Residues[j].xNeutronScattering;
                    ProteinStructureCopy.Residues[j].yNeutronScattering = ProteinStructure.Residues[j].yNeutronScattering;
                    ProteinStructureCopy.Residues[j].zNeutronScattering = ProteinStructure.Residues[j].zNeutronScattering;

                    ProteinStructureCopy.Residues[j].XRayScatteringLength    = ProteinStructure.Residues[j].XRayScatteringLength;
                    ProteinStructureCopy.Residues[j].NeutronScatteringLength = ProteinStructure.Residues[j].NeutronScatteringLength;
                    ProteinStructureCopy.Residues[j].Volume                  = ProteinStructure.Residues[j].Volume;
                    ProteinStructureCopy.Residues[j].ResidueID               = ProteinStructure.Residues[j].ResidueID;

                    strcpy(ProteinStructureCopy.Residues[j].Name, ProteinStructure.Residues[j].Name);
                }

                for (j = 0; j < ProteinStructureCopy.NumberOfAtoms; ++j) {
                    ProteinStructureCopy.Atoms[j].x = ProteinStructure.Atoms[j].x;
                    ProteinStructureCopy.Atoms[j].y = ProteinStructure.Atoms[j].y;
                    ProteinStructureCopy.Atoms[j].z = ProteinStructure.Atoms[j].z;

                    ProteinStructureCopy.Atoms[j].XRayScatteringLength    = ProteinStructure.Atoms[j].XRayScatteringLength;
                    ProteinStructureCopy.Atoms[j].NeutronScatteringLength = ProteinStructure.Atoms[j].NeutronScatteringLength;
                    ProteinStructureCopy.Atoms[j].Volume                  = ProteinStructure.Atoms[j].Volume;
					ProteinStructureCopy.Atoms[j].Name                    = ProteinStructure.Atoms[j].Name;
                }
			}

            Population[i].Chisquare = ComputeChiSquareGenetic(Data, NumberOfSpectra, DummyParameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules,
                                                              ProteinStructureCopy, &UserDefinedCopy, TotalNumberOfDatapoints, NumberOfFreeParameters);

            free(DummyParameters);

            for (j = 0; j < NumberOfSpectra; ++j) {
                free(DummyData[j].QValues);
                free(DummyData[j].IValues);
                free(DummyData[j].FitValues);
                free(DummyData[j].SigmaValues);
                free(DummyData[j].SigmaQValues);
                free(DummyData[j].Constraints);
                free(DummyData[j].ScatteringLengths);

                for (k = 0; k < HighestNumberOfDatapoints; ++k) {
                    free(DummyData[j].ResolutionWeights[k]);
                }

                free(DummyData[j].ResolutionWeights);
            }

            FreeUserDefined(&UserDefinedCopy);
            free(DummyData);

            if (ProteinStructure.NumberOfAtoms != 0) {
                free(ProteinStructureCopy.Residues);
                free(ProteinStructureCopy.Atoms);
            }
        }

        if (DoExcessivePrinting == true) {

            for (i = 0; i < SizeOfPopulation; ++i) {
                fprintf(Outfile, "%15g,  ", Population[i].Chisquare);

                for (j = 0; j < NumberOfParameters; ++j) {
                    fprintf(Outfile, "%15g,  ", Population[i].Parameters[j]);
                }

                fprintf(Outfile, "\n");
            }
        }

        // Compare chisquares
        for (i = 0; i < SizeOfPopulation; ++i) {

            if (BestCandidate.Chisquare > Population[i].Chisquare) {
                BestCandidate.Chisquare = Population[i].Chisquare;

                for (j = 0; j < NumberOfParameters; ++j) {
                    BestCandidate.Parameters[j] = Population[i].Parameters[j];
                }
            }
        }

        // Print current best candidate
        if (Generation % 5 == 0) {
            PrintBestCandidate(BestCandidate, Generation, NumberOfParameters, Parameters);
        }

        // Evaluate total fitness of current population
        SumOfFitness = 0.0;

        for (i = 0; i < SizeOfPopulation; ++i) {
            SumOfFitness += EvaluateFitness(Population[i].Chisquare);
            FitnessRanking[i] = SumOfFitness;
        }

        // Assign normalized fitness for population
        for (i = 0; i < SizeOfPopulation; ++i) {
            Population[i].Fitness = EvaluateFitness(Population[i].Chisquare) / SumOfFitness;
            FitnessRanking[i] /= SumOfFitness;
        }

        // Breed new population
        for (i = 0; i < SizeOfPopulation; ++i) {

            // Search for first parent
            j = 0;

            Criteria = RandomFraction();

            while (FitnessRanking[j] < Criteria) {
                ++j;
            }

            for (k = 0; k < NumberOfParameters; ++k) {
                ParentOne.Parameters[k] = Population[j].Parameters[k];
            }

            // Search for second parent
            j = 0;

            Criteria = RandomFraction();

            while (FitnessRanking[j] < Criteria) {
                ++j;
            }

            for (k = 0; k < NumberOfParameters; ++k) {
                ParentTwo.Parameters[k] = Population[j].Parameters[k];
            }

            // Create new individual
            CreateIndividual(&NewPopulation[i], ParentOne, ParentTwo, NumberOfParameters);
        }

        // Mutate new population
        ProbabilityOfMutation = 0.8 * Generation / (1.0 * MaxNumberOfGenerations) + 0.1;

        for (i = 0; i < SizeOfPopulation; ++i) {
            MutateSubject(&NewPopulation[i], Parameters, NumberOfParameters, ProbabilityOfMutation);
        }

        // Replace the old population
        for (i = 0; i < SizeOfPopulation; ++i) {
            Population[i].Chisquare = NewPopulation[i].Chisquare;
            Population[i].Fitness   = NewPopulation[i].Fitness;

            for (j = 0; j < NumberOfParameters; ++j) {
                Population[i].Parameters[j] = NewPopulation[i].Parameters[j];
            }
        }

        // Conclusion of loop
        ++Generation;

        if (Generation > MaxNumberOfGenerations) {
            Stop = true;
        }
    }

    // Conclusion
    if (DoExcessivePrinting == true) {
        fclose(Outfile);
    }

    for (i = 0; i < NumberOfParameters; ++i) {
        Parameters[i].Value = BestCandidate.Parameters[i];
    }

    *Chisquare = ComputeChiSquareGenetic(Data, NumberOfSpectra, Parameters, NumberOfParameters, NumberOfSmearingFolds, VolumesOfMolecules, ProteinStructure, &*UserDefinedStructure,
                                         TotalNumberOfDatapoints, NumberOfFreeParameters);

    // Free subjects from memory
    free(ParentOne.Parameters);
    free(ParentTwo.Parameters);
    free(BestCandidate.Parameters);

    for (i = 0; i < SizeOfPopulation; ++i) {
        free(NewPopulation[i].Parameters);
        free(Population[i].Parameters);
    }

    return 0;
}
