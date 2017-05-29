/// Declaration of Dataset-structure
struct Dataset {
    double * QValues;
    double * IValues;
    double * SigmaValues;
    double * SigmaQValues;
    double ** ResolutionWeights;
    double * FitValues;
    double Concentration;
    double Contrast;
    double * Constraints;
    double * ScatteringLengths;
    bool IncludeResolutionEffects;
    int NMin;
    int NMax;
    int NumberOfDatapoints;
    char Filename[256];
};

/// Declaration of parameter-structure
struct Parameter {
    double Value;
    bool iParameter;
    double MinValue;
    double MaxValue;
    double Error;
    char Name[64];
};

/// Declaration of structures used for PDB-files
struct Atom {
    // Coordinates
    double x;
    double y;
    double z;

    // Physical properties
    double XRayScatteringLength;
    double NeutronScatteringLength;
    double Volume;

    // Book-keeping
    char Name;
};

struct Residue {
    // Coordinates of center of volume
    double xVolume;
    double yVolume;
    double zVolume;

    double xXRayScattering;
    double yXRayScattering;
    double zXRayScattering;

    double xNeutronScattering;
    double yNeutronScattering;
    double zNeutronScattering;

    // Physical properties
    double XRayScatteringLength;
    double NeutronScatteringLength;
    double Volume;

    // Book-keeping
    char Name[5];
    int ResidueID;
};

struct Protein {
    struct Atom * Atoms;
    struct Residue * Residues;
    int NumberOfAtoms;
    int NumberOfResidues;
    char PDBFileLocation[256];
};
