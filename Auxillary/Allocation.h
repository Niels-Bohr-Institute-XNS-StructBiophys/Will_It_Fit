/// Function used to allocate the data struct
void AllocateData(struct Dataset ** Data, int NumberOfSpectra)
{
    // Declarations
    struct Dataset * DummyData;

    // Initialization
    DummyData = (struct Dataset *) calloc(NumberOfSpectra, sizeof(struct Dataset));

    *Data = DummyData;
}

/// Function used to allocate the parameter struct
void AllocateParameters(struct Parameter ** Parameters, int NumberOfParameters)
{
    // Declarations
    struct Parameter * ParametersDummy;

    // Initialization
    ParametersDummy = (struct Parameter *) calloc(NumberOfParameters, sizeof(struct Parameter));

    *Parameters = ParametersDummy;
}

/// Function used to allocate the Protein struct
void AllocateProteinStructure(struct Protein * ProteinStructure, int NumberOfResidues, int NumberOfAtoms)
{
    ProteinStructure->Residues = (struct Residue *) calloc(NumberOfResidues, sizeof(struct Residue));
    ProteinStructure->Atoms    = (struct Atom *)    calloc(NumberOfAtoms,    sizeof(struct Atom));
}

/// Functions used to allocate memory for arrays
void Initialize1DArray(double ** Array, int Length)
{
    // Declarations
    int i;
    double * DummyArray;

    // Initialization
    DummyArray = (double *) calloc(Length, sizeof(double));

    for (i = 0; i < Length; ++i) {
        DummyArray[i] = 0.0;
    }

    *Array = DummyArray;
}

void Initialize1DIntegerArray(int ** Array, int Length)
{
    // Declarations
    int i;
    int * DummyArray;

    // Initialization
    DummyArray = (int *) calloc(Length, sizeof(int));

    for (i = 0; i < Length; ++i) {
        DummyArray[i] = 0.0;
    }

    *Array = DummyArray;
}

void Initialize2DArray(double *** Array, int Length, int Width)
{
    // Declarations
    int i;
    int j;
    double ** DummyArray;

    // Initialization
    DummyArray = (double **) calloc(Length, sizeof(double));

    for (i = 0; i < Length; ++i) {
        DummyArray[i] = (double *) calloc(Width, sizeof(double));

        for (j = 0; j < Width; ++j) {
            DummyArray[i][j] = 0.0;
        }
    }

    *Array = DummyArray;
}

void Free2DArray(double ** Array, int Length) {
    // Declaration
    int i;

    // Emptying memory
    for (i = 0; i < Length; ++i) {
        free(Array[i]);
    }

    free(Array);
}
