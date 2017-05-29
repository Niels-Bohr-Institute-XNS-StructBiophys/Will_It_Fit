int CheckNumberOfResiduesInPDBFile(char Filename[256])
{
    // Declarations
    FILE *PointerToFile;
    char Linebuffer[256];
    double xDummy;
    double yDummy;
    double zDummy;
    int NumberOfResidues = 0;
    int ResidueID = 0;
    int PreviousResidueID = 0;

    // I/O
    PointerToFile = fopen(Filename, "r");

    if(PointerToFile == 0){
        printf("An error occured when attempting to open PDB-file. \n");
        return -1;
    }

    while(fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {
        ResidueID = 0;

        if(sscanf(Linebuffer, "ATOM%*18c%d%*4c%lf%lf%lf", &ResidueID, &xDummy, &yDummy, &zDummy) == 4){

            if(ResidueID != PreviousResidueID && ResidueID != 0) {
                ++NumberOfResidues;
            }

            PreviousResidueID = ResidueID;
        }
        if(sscanf(Linebuffer, "HETATM%d%*18c%*4c%lf%lf%lf", &ResidueID, &xDummy, &yDummy, &zDummy) == 4){
            if(ResidueID != PreviousResidueID && ResidueID != 0) {
                ++NumberOfResidues;
            }

            PreviousResidueID = ResidueID;
        }
    }

    fclose(PointerToFile);

    return NumberOfResidues;
}

int CheckNumberOfAtomsInPDBFile(char Filename[256])
{
    // Declarations
    FILE *PointerToFile;
    char Linebuffer[256];
    double xDummy;
    double yDummy;
    double zDummy;
    int NumberOfAtoms = 0;
    int DummyID;

    // I/O
    PointerToFile = fopen(Filename, "r");

    if(PointerToFile == 0){
        printf("An error occured when attempting to open PDB-file. \n");
        return -1;
    }

    while(fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {

        if(sscanf(Linebuffer, "ATOM%*18c%d%*4c%lf%lf%lf", &DummyID, &xDummy, &yDummy, &zDummy) == 4){
            ++NumberOfAtoms;
        }
        if(sscanf(Linebuffer, "HETATM%*18c%d%*4c%lf%lf%lf", &DummyID, &xDummy, &yDummy, &zDummy) == 4){
            ++NumberOfAtoms;
        }
    }

    fclose(PointerToFile);

    return NumberOfAtoms;
}

void ImportAtomsFromPDBFile(char Filename[256], struct Protein ProteinStruct, int NumberOfAtoms)
{
    // Declarations
    FILE *PointerToFile;
    char Linebuffer[256];
    double xDummy;
    double yDummy;
    double zDummy;
    char Dummychar;
    char AtomName;
    int DummyResidueID;
    int IDOfCurrentAtom = 0;

    // Physical values
    const double HVolume = 5.15;
    const double DVolume = 5.15;
    const double CVolume = 16.44;
    const double NVolume = 2.49;
    const double OVolume = 9.13;
    const double PVolume = 5.73;
    const double SVolume = 19.86;
    const double H2OVolume =  30.0;
    const double CLVolume = 28.81;
    const double ZNVolume = 9.85;


    const double HXRayScatteringLength =  1.0 * 2.82e-13;
    const double DXRayScatteringLength =  1.0 * 2.82e-13;
    const double CXRayScatteringLength =  6.0 * 2.82e-13;
    const double NXRayScatteringLength =  7.0 * 2.82e-13;
    const double OXRayScatteringLength =  8.0 * 2.82e-13;
    const double PXRayScatteringLength = 15.0 * 2.82e-13;
    const double SXRayScatteringLength = 16.0 * 2.82e-13;
    const double H2OXRayScatteringLength =  10.0 * 2.82e-13;
    const double CLXRayScatteringLength = 17.0 * 2.82e-13;
    const double ZNXRayScatteringLength = 30.0 * 2.82e-13;


    const double HNeutronScatteringLength = -3.742e-13;
    const double DNeutronScatteringLength = 6.674e-13;
    const double CNeutronScatteringLength = 6.6484e-13;
    const double NNeutronScatteringLength = 9.36e-13;
    const double ONeutronScatteringLength = 5.803e-13;
    const double PNeutronScatteringLength = 5.13E-13;
    const double SNeutronScatteringLength = 2.847e-13;
    const double CLNeutronScatteringLength = 5.13E-13;
    const double ZNNeutronScatteringLength = 2.847e-13;
    const double GlycanScale = 1.0;
    // I/O
    PointerToFile = fopen(Filename, "r");

    while (fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {

        if (sscanf(Linebuffer, "ATOM%*9c%c%*8c%d%*4c%lf%lf%lf%*23c%c", &Dummychar, &DummyResidueID, &xDummy, &yDummy, &zDummy, &AtomName) == 6) {
            ProteinStruct.Atoms[IDOfCurrentAtom].x    = xDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].y    = yDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].z    = zDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].Name = AtomName;

            switch (AtomName) {
                case 'H':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = HXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = HNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = HVolume;
                break;

                case 'D':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = DXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = DNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = DVolume;
                break;

                case 'C':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = CXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = CNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = CVolume;
                break;

                case 'N':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = NXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = NNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = NVolume;
                break;

                case 'O':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = OXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = ONeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = OVolume;
                break;

                case 'P':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = PXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = PNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = PVolume;
                break;

                case 'S':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = SXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = SNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = SVolume;
                break;

                case 'Q':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength = 4.092*H2OXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = SNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume  = 4.092*(H2OVolume);
                break;

                default:
                    printf("Encountered unknown atom in PDB-file - labelled %c at site %d.\n", AtomName, IDOfCurrentAtom);
                break;
            }

            ++IDOfCurrentAtom;
        }
        if (sscanf(Linebuffer, "HETATM%*9c%c%*8c%d%*4c%lf%lf%lf%*23c%c", &Dummychar, &DummyResidueID, &xDummy, &yDummy, &zDummy, &AtomName) == 6) {
            ProteinStruct.Atoms[IDOfCurrentAtom].x    = xDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].y    = yDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].z    = zDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].Name = AtomName;

            switch (AtomName) {
                case 'H':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = HXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = HNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*HVolume;
                break;

                case 'D':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = DXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = DNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*DVolume;
                break;

                case 'C':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = CXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = CNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*CVolume;
                break;

                case 'N':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = NXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = NNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*NVolume;
                break;

                case 'O':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = OXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = ONeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*OVolume;
                break;

                case 'P':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = PXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = PNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*PVolume;
                break;

                case 'S':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength    = SXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = SNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume                  = GlycanScale*SVolume;
                break;

                case 'Q':
                    ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength = 4.092*H2OXRayScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength = SNeutronScatteringLength;
                    ProteinStruct.Atoms[IDOfCurrentAtom].Volume  = 4.092*(H2OVolume);
                break;

                default:
                    printf("Encountered unknown atom in PDB-file - labelled %c at site %d.\n", AtomName, IDOfCurrentAtom);
                break;
            }

            ++IDOfCurrentAtom;
        }
    }

    fclose(PointerToFile);
}

void ImportResiduesFromPDBFile(char Filename[256], struct Protein ProteinStruct, int NumberOfResidues)
{
    // Declarations
    FILE *PointerToFile;
    char Linebuffer[256];
    char AtomName[2];
    char Dummychar;
    int PreviousResidueID = 0;
    int ResidueID = 0;
    int IDOfCurrentResidue = 0;
    int AtomRecg = 0; // value for checking that atom is found in database

    double VolumeOfResidue = 0.0;
    double XRayScatteringLengthOfResidue = 0.0;
    double NeutronScatteringLengthOfResidue = 0.0;

    double xCenterOfXRayScattering = 0.0;
    double yCenterOfXRayScattering = 0.0;
    double zCenterOfXRayScattering = 0.0;

    double xCenterOfNeutronScattering = 0.0;
    double yCenterOfNeutronScattering = 0.0;
    double zCenterOfNeutronScattering = 0.0;

    double xCenterOfVolume = 0.0;
    double yCenterOfVolume = 0.0;
    double zCenterOfVolume = 0.0;

    double VolumeOfCurrentAtom;
    double XRayScatteringLengthOfCurrentAtom;
    double NeutronScatteringLengthOfCurrentAtom;

    double xDummy;
    double yDummy;
    double zDummy;

    // Physical values
    const double HVolume = 5.15;
    const double DVolume = 5.15;
    const double CVolume = 16.44;
    const double NVolume = 2.49;
    const double OVolume = 9.13;
    const double PVolume = 5.73;
    const double SVolume = 19.86;
    const double H2OVolume =  30.0;
    const double CLVolume = 28.81;
    const double ZNVolume = 9.85;


    const double HXRayScatteringLength =  1.0 * 2.82e-13;
    const double DXRayScatteringLength =  1.0 * 2.82e-13;
    const double CXRayScatteringLength =  6.0 * 2.82e-13;
    const double NXRayScatteringLength =  7.0 * 2.82e-13;
    const double OXRayScatteringLength =  8.0 * 2.82e-13;
    const double PXRayScatteringLength = 15.0 * 2.82e-13;
    const double SXRayScatteringLength = 16.0 * 2.82e-13;
    const double H2OXRayScatteringLength =  10.0 * 2.82e-13;
    const double CLXRayScatteringLength = 17.0 * 2.82e-13;
    const double ZNXRayScatteringLength = 30.0 * 2.82e-13;


    const double HNeutronScatteringLength = -3.742e-13;
    const double DNeutronScatteringLength = 6.674e-13;
    const double CNeutronScatteringLength = 6.6484e-13;
    const double NNeutronScatteringLength = 9.36e-13;
    const double ONeutronScatteringLength = 5.803e-13;
    const double PNeutronScatteringLength = 5.13E-13;
    const double SNeutronScatteringLength = 2.847e-13;
    const double CLNeutronScatteringLength = 5.13E-13;
    const double ZNNeutronScatteringLength = 2.847e-13;
    const double GlycanScale = 1.0;
//    const double SNeutronScatteringLength = 10.0 * 2.82e-13;


    // I/O
    PointerToFile = fopen(Filename, "r");

    while (fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) { //Read File, defined above, line by line and store it in "Linebuffer"
        ResidueID = 0;

        // **** These if statements are for testing what gets read
	if (sscanf(Linebuffer, "ATOM%*9c%c%*8c%d%*4c%lf%lf%lf%*22c%2c", &Dummychar, &ResidueID, &xDummy, &yDummy, &zDummy, &AtomName) == 6) {
        //Check that we read
        //printf("%s",Linebuffer);
        //printf("%i %s %c %i\n",ResidueID,AtomName, Dummychar,IDOfCurrentResidue);
        }
        if (sscanf(Linebuffer, "HETATM%d%*3c%c%*8c%*4c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &xDummy, &yDummy, &zDummy, &AtomName) == 6) {
        //Check that we read
        //printf("%s",Linebuffer);
        //printf("%i %s %c %i\n",ResidueID,AtomName, Dummychar,IDOfCurrentResidue);
	}
	// ****
	
        if (ResidueID != PreviousResidueID) {

            if (PreviousResidueID != 0) {

                ProteinStruct.Residues[IDOfCurrentResidue].Volume = VolumeOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].XRayScatteringLength = XRayScatteringLengthOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].NeutronScatteringLength = NeutronScatteringLengthOfResidue;

                ProteinStruct.Residues[IDOfCurrentResidue].xVolume = xCenterOfVolume / VolumeOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].yVolume = yCenterOfVolume / VolumeOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].zVolume = zCenterOfVolume / VolumeOfResidue;

                ProteinStruct.Residues[IDOfCurrentResidue].xXRayScattering = xCenterOfXRayScattering / XRayScatteringLengthOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].yXRayScattering = yCenterOfXRayScattering / XRayScatteringLengthOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].zXRayScattering = zCenterOfXRayScattering / XRayScatteringLengthOfResidue;

                ProteinStruct.Residues[IDOfCurrentResidue].xNeutronScattering = xCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].yNeutronScattering = yCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;
                ProteinStruct.Residues[IDOfCurrentResidue].zNeutronScattering = zCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;

                xCenterOfVolume = 0.0;
                yCenterOfVolume = 0.0;
                zCenterOfVolume = 0.0;
                VolumeOfResidue = 0.0;

                xCenterOfXRayScattering = 0.0;
                yCenterOfXRayScattering = 0.0;
                zCenterOfXRayScattering = 0.0;
                XRayScatteringLengthOfResidue = 0.0;

                xCenterOfNeutronScattering = 0.0;
                yCenterOfNeutronScattering = 0.0;
                zCenterOfNeutronScattering = 0.0;
                NeutronScatteringLengthOfResidue = 0.0;

                ++IDOfCurrentResidue;
            }

            PreviousResidueID = ResidueID;
        }

        if (ResidueID == NumberOfResidues - 1) {
            ProteinStruct.Residues[IDOfCurrentResidue].Volume = VolumeOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].XRayScatteringLength = XRayScatteringLengthOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].NeutronScatteringLength = NeutronScatteringLengthOfResidue;

            ProteinStruct.Residues[IDOfCurrentResidue].xVolume = xCenterOfVolume / VolumeOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].yVolume = yCenterOfVolume / VolumeOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].zVolume = zCenterOfVolume / VolumeOfResidue;

            ProteinStruct.Residues[IDOfCurrentResidue].xXRayScattering = xCenterOfXRayScattering / XRayScatteringLengthOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].yXRayScattering = yCenterOfXRayScattering / XRayScatteringLengthOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].zXRayScattering = zCenterOfXRayScattering / XRayScatteringLengthOfResidue;

            ProteinStruct.Residues[IDOfCurrentResidue].xNeutronScattering = xCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].yNeutronScattering = yCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;
            ProteinStruct.Residues[IDOfCurrentResidue].zNeutronScattering = zCenterOfNeutronScattering / NeutronScatteringLengthOfResidue;
        }

        if (sscanf(Linebuffer, "ATOM%*13c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name)== 1) {
//        printf("Atom %i ",sscanf(Linebuffer, "ATOM%*13c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name) );
//        printf("%c \n",AtomName );

        // the following if statements are to read the atoms. Note this was done in a swich statement before
        // but that does not work for multi character atom names such as "Zn"
        if (AtomName[0] == ' ' && AtomName[1] == 'H' ){
          XRayScatteringLengthOfCurrentAtom = HXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = HNeutronScatteringLength;
          VolumeOfCurrentAtom = HVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'D' ){
          XRayScatteringLengthOfCurrentAtom = DXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = DNeutronScatteringLength;
          VolumeOfCurrentAtom = DVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'C' ){
          XRayScatteringLengthOfCurrentAtom = CXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = CNeutronScatteringLength;
          VolumeOfCurrentAtom = CVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'N' ){
          XRayScatteringLengthOfCurrentAtom = NXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = NNeutronScatteringLength;
          VolumeOfCurrentAtom = NVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'O' ){
          XRayScatteringLengthOfCurrentAtom = OXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = ONeutronScatteringLength;
          VolumeOfCurrentAtom = OVolume;
          AtomRecg = 1;
            }
        if (AtomName[0] == ' ' && AtomName[1] == 'P' ){
          XRayScatteringLengthOfCurrentAtom = PXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = PNeutronScatteringLength;
          VolumeOfCurrentAtom = PVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'S' ){
          XRayScatteringLengthOfCurrentAtom = SXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          VolumeOfCurrentAtom = SVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'Q' ){
          XRayScatteringLengthOfCurrentAtom = 4.092*H2OXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          VolumeOfCurrentAtom = 4.092*(H2OVolume);
          AtomRecg = 1;
        }
        if (AtomName[0] == 'Z' && AtomName[1] == 'N' ){
          XRayScatteringLengthOfCurrentAtom = ZNXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = ZNNeutronScatteringLength;
          VolumeOfCurrentAtom = ZNVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == 'C' && AtomName[1] == 'L' ){
          XRayScatteringLengthOfCurrentAtom = CLXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = CLNeutronScatteringLength;
          VolumeOfCurrentAtom = CLVolume;
          AtomRecg = 1;
        }
        if (AtomRecg == 0){
          printf("Atom name %c%c not found in database\n", AtomName[0],AtomName[1]);
        }
        AtomRecg = 0;
      }
        if (sscanf(Linebuffer, "HETATM%*11c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name)== 1) {
//        printf("hetatm %i ",sscanf(Linebuffer, "HETATM%*13c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name) );
          if (ProteinStruct.Residues[IDOfCurrentResidue].Name[2] != '9' ){

         ProteinStruct.Residues[IDOfCurrentResidue].Name[0] = 'X';
         }
        // the following if statements are to read te atoms. Note this was done in a swich statement before
        // but that does not work for multi character atom names such as "Zn"
        if (AtomName[0] == ' ' && AtomName[1] == 'H' ){
          XRayScatteringLengthOfCurrentAtom = HXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = HNeutronScatteringLength;
          VolumeOfCurrentAtom = HVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'D' ){
          XRayScatteringLengthOfCurrentAtom = DXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = DNeutronScatteringLength;
          VolumeOfCurrentAtom = DVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'C' ){
          XRayScatteringLengthOfCurrentAtom = CXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = CNeutronScatteringLength;
          VolumeOfCurrentAtom = CVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'N' ){
          XRayScatteringLengthOfCurrentAtom = NXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = NNeutronScatteringLength;
          VolumeOfCurrentAtom = NVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'O' ){
          XRayScatteringLengthOfCurrentAtom = OXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = ONeutronScatteringLength;
          VolumeOfCurrentAtom = OVolume;
          AtomRecg = 1;
            }
        if (AtomName[0] == ' ' && AtomName[1] == 'P'){
          XRayScatteringLengthOfCurrentAtom = PXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = PNeutronScatteringLength;
          VolumeOfCurrentAtom = PVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'S' ){
          XRayScatteringLengthOfCurrentAtom = SXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          VolumeOfCurrentAtom = SVolume;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'Q' ){
          XRayScatteringLengthOfCurrentAtom = 4.092*H2OXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          VolumeOfCurrentAtom = 4.092*(H2OVolume);
          AtomRecg = 1;
        }

        if (AtomName[0] == 'Z' && AtomName[1] == 'N' ){
          XRayScatteringLengthOfCurrentAtom = ZNXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = ZNNeutronScatteringLength;
          VolumeOfCurrentAtom = ZNVolume;
          AtomRecg = 1;
        }

        if (AtomName[0] == 'C' && AtomName[1] == 'L' ){
          XRayScatteringLengthOfCurrentAtom = CLXRayScatteringLength;
          NeutronScatteringLengthOfCurrentAtom = CLNeutronScatteringLength;
          VolumeOfCurrentAtom = CLVolume;
          AtomRecg = 1;
        }

        if (AtomRecg == 0){
          printf("Atom name %c%c not found in database\n", AtomName[0],AtomName[1]);
          printf("%s \n", Linebuffer);
        }
        AtomRecg = 0;
      }

        VolumeOfResidue += VolumeOfCurrentAtom;
        XRayScatteringLengthOfResidue += XRayScatteringLengthOfCurrentAtom;
        NeutronScatteringLengthOfResidue += NeutronScatteringLengthOfCurrentAtom;

        xCenterOfVolume += VolumeOfCurrentAtom * xDummy;
        yCenterOfVolume += VolumeOfCurrentAtom * yDummy;
        zCenterOfVolume += VolumeOfCurrentAtom * zDummy;

        xCenterOfXRayScattering += XRayScatteringLengthOfCurrentAtom * xDummy;
        yCenterOfXRayScattering += XRayScatteringLengthOfCurrentAtom * yDummy;
        zCenterOfXRayScattering += XRayScatteringLengthOfCurrentAtom * zDummy;

        xCenterOfNeutronScattering += NeutronScatteringLengthOfCurrentAtom * xDummy;
        yCenterOfNeutronScattering += NeutronScatteringLengthOfCurrentAtom * yDummy;
        zCenterOfNeutronScattering += NeutronScatteringLengthOfCurrentAtom * zDummy;
    }

    fclose(PointerToFile);
}
