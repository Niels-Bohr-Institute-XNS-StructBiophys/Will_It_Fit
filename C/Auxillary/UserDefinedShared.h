/*
 * This file offers a model the option to transport a structure through
 * the entire architecture of the program.
 * 
 * The file must include a definition of the UserDefinedSheared-structure as well
 * as the two associated functions, which are needed to correctly
 * (de)allocate the structs correctly.
 * 
 * By contrast to UserDefined this structure will be shared among
 * parallel threads, thus changes should not be applied within parallelized sections
 * 
 * In order to include this header in the compilation of the program,
 * please put in next to the ModelInfo.h-header (and do not write any
 * #include-statements anywhere)
 */

/// Structure used to transport user-defined variables to the model
struct UserDefinedShared {
	bool DummyBool;
	int DummyInt;
	double DummyDouble;
	double ** DummyDoublePointer;
};

/// Function initializing the structure above - this function is executed once during the initialisation part of the program
void InitializeUserDefinedSharedStructure(struct UserDefined * UserDefinedSharedStruct)
{
	// Declarations
	int i;
	int j;
	const double Length = 9;
	const double Width  = 5;

	// Allocation and initialization
	UserDefinedSharedStruct->DummyBool   = false;
	UserDefinedSharedStruct->DummyInt    = 0;
	UserDefinedSharedStruct->DummyDouble = 0.0;

	UserDefinedSharedStruct->DummyDoublePointer = (double **) calloc(Length, sizeof(double));

	for (i = 0; i < Length; ++i)
	{
		UserDefinedSharedStruct->DummyDoublePointer[i] = (double *) calloc(Width, sizeof(double));

		for (j = 0; j < Width; ++j)
		{
			UserDefinedSharedStruct->DummyDoublePointer[i][j] = 0.0;
		}
	}
}

/// Function freeing eventual pointers in a structure - this must be supplied, if pointers are present in the structure to prevent memory overflows
void FreeUserDefinedShared(struct UserDefinedShared * UserDefinedSharedStruct)
{
	// Declarations
	int i;
	int Length = 9;

	// Free
	for (i = 0; i < Length; ++i) {
	free(UserDefinedSharedStruct->DummyDoublePointer[i]);
	}

	free(UserDefinedSharedStruct->DummyDoublePointer);
}
