/*
 * This file offers a model the option to transport a structure through
 * the entire architecture of the program.
 * 
 * The file must include a definitoin of the UserDefined-struct as well
 * as the three associated functions, which are needed to correctly
 * allocate the structs correctly.
 * 
 * In order to include this header in the compilation of the program,
 * please put in next to the ModelInfo.h-header (and do not write any
 * #include-statements anywhere)
 */

/// Structure used to transport user-defined variables to the model
struct UserDefined {
    bool DummyBool;
    int DummyInt;
    double DummyDouble;
    double ** DummyDoublePointer;
};

/// Function initializing the structure above - this function is executed once during the initialisation part of the program
void InitializeUserDefinedStructure(struct UserDefined * UserDefinedStruct)
{
    // Declarations
    int i;
    int j;
    const double Length = 9;
    const double Width  = 5;

    // Allocation and initialization
    UserDefinedStruct->DummyBool   = false;
    UserDefinedStruct->DummyInt    = 0;
    UserDefinedStruct->DummyDouble = 0.0;

    UserDefinedStruct->DummyDoublePointer = (double **) calloc(Length, sizeof(double));

    for (i = 0; i < Length; ++i) {

        UserDefinedStruct->DummyDoublePointer[i] = (double *) calloc(Width, sizeof(double));

        for (j = 0; j < Width; ++j) {
            UserDefinedStruct->DummyDoublePointer[i][j] = 0.0;
        }
    }
}

/// Function creating an identical (already allocated) copy of a structure - used for parallelisation
void CopyUserDefined(struct UserDefined * OriginalStruct, struct UserDefined * CopyOfStruct)
{
    // Declarations
    int i;
    int j;
    const double Length = 9;
    const double Width  = 5;

    // Copy
    CopyOfStruct->DummyBool   = OriginalStruct->DummyBool;
    CopyOfStruct->DummyInt    = OriginalStruct->DummyInt;
    CopyOfStruct->DummyDouble = OriginalStruct->DummyDouble;

    for (i = 0; i < Length; ++i) {

        for (j = 0; j < Width; ++j) {
            CopyOfStruct->DummyDoublePointer[i][j] = OriginalStruct->DummyDoublePointer[i][j];
        }
    }
}

/// Function freeing eventual pointers in a structure - this must be supplied, if pointers are present in the structure to prevent memory overflows
void FreeUserDefined(struct UserDefined * UserDefinedStruct)
{
    // Declarations
    int i;
    int Length = 9;

    // Freeubg
    for (i = 0; i < Length; ++i) {
        free(UserDefinedStruct->DummyDoublePointer[i]);
    }

    free(UserDefinedStruct->DummyDoublePointer);
}
