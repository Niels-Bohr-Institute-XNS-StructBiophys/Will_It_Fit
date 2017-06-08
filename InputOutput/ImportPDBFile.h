int CheckNumberOfResiduesInPDBFile(char Filename[256])
{
    // Declarations
    FILE *PointerToFile;
    char Linebuffer[256];
    double xDummy;
    double yDummy;
    double zDummy;
    int NumberOfResidues = 0;
    int PreviousResidueID = 0;
    int ResidueID;
    char AtomName[2];
    char ResidueName[3];
    char Dummychar;
   // I/O
    PointerToFile = fopen(Filename, "r");

    if(PointerToFile == 0){
        printf("An error occured when attempting to open PDB-file. \n");
        return -1;
    }

    while(fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {
        ResidueID = 0;
	//if (sscanf(Linebuffer, "ATOM%*9c%c%*3c%3c%*2c%d%*c%lf%lf%lf%*22c%2c", &Dummychar, &ResidueName,&ResidueID, &xDummy, &yDummy, &zDummy, &AtomName) == 7) { //This Line for reading amino acids
        
	if (sscanf(Linebuffer, "ATOM%d%*3c%c%*2c%3c%*6c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &ResidueName[0], &xDummy, &yDummy, &zDummy, &AtomName[0]) == 7) { // This line for reading individual atoms

            if(ResidueID != PreviousResidueID && ResidueID != 0) {
                ++NumberOfResidues;
            }

            PreviousResidueID = ResidueID;
        }
	if (sscanf(Linebuffer, "HETATM%d%*3c%c%*12c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &xDummy, &yDummy, &zDummy, &AtomName[0]) == 6) {

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
    char AtomName[2];
    char ResidueName[3];
    char Dummychar;
    int ResidueID;

    // I/O
    PointerToFile = fopen(Filename, "r");

    if(PointerToFile == 0){
        printf("An error occured when attempting to open PDB-file. \n");
        return -1;
    }

    while(fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {
        
	if (sscanf(Linebuffer, "ATOM%5d%*2c%c%*3c%3c%*6c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &ResidueName[0], &xDummy, &yDummy, &zDummy, &AtomName[0]) == 7) { // This line for reading individual atoms
//	printf("%d %s %s \n",ResidueID,AtomName,ResidueName);
            ++NumberOfAtoms;
        }
	if (sscanf(Linebuffer, "HETATM%5d%*2c%c%*3c%3c%*6c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &ResidueName[0], &xDummy, &yDummy, &zDummy, &AtomName[0]) == 7) { // This line for reading individual atoms
//	printf("%d %c%c %s \n",ResidueID,AtomName[0],AtomName[1],ResidueName);
            ++NumberOfAtoms;
        }
    }

    fclose(PointerToFile);

    return NumberOfAtoms;
}

void AssignAtom(char AtomName[2], double *XRayScatteringLengthOfCurrentAtom,double  *NeutronScatteringLengthOfCurrentAtom,double *VolumeOfCurrentAtom, double *WeightOfCurrentAtom)
{  // function for assigning the scattering length to the correct atoms
   int AtomRecg = 0;

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

    const double HWeight = 1.008;
    const double CWeight = 12.011;
    const double NWeight = 14.007;  
    const double OWeight = 15.999;
    const double PWeight = 30.974;
    const double SWeight = 32.06 ;
    const double CLWeight = 35.45;
    const double ZNWeight = 65.38;

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

// the following if statements are to read the atoms. Note this was done in a swich statement before
        // but that does not work for multi character atom names such as "Zn"
        if (AtomName[0] ==  ' ' && AtomName[1] == 'H' ){
          *XRayScatteringLengthOfCurrentAtom = HXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = HNeutronScatteringLength;
          *VolumeOfCurrentAtom = HVolume;
	  *WeightOfCurrentAtom = HWeight;
          AtomRecg = 1;
	}
        if (strcmp(AtomName, " D")== 0 ){
          *XRayScatteringLengthOfCurrentAtom = DXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = DNeutronScatteringLength;
          *VolumeOfCurrentAtom = DVolume;
          *WeightOfCurrentAtom = HWeight;
	  AtomRecg = 1;

        }
        if (AtomName[0] == ' ' && AtomName[1] == 'C' ){
         *XRayScatteringLengthOfCurrentAtom = CXRayScatteringLength;
         *NeutronScatteringLengthOfCurrentAtom = CNeutronScatteringLength;
         *VolumeOfCurrentAtom = CVolume;
         *WeightOfCurrentAtom = CWeight;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'N' ){
          *XRayScatteringLengthOfCurrentAtom = NXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = NNeutronScatteringLength;
          *VolumeOfCurrentAtom = NVolume;
          *WeightOfCurrentAtom = NWeight;
          AtomRecg = 1;
        }

        if (AtomName[0] == ' ' && AtomName[1] == 'O' ){
          *XRayScatteringLengthOfCurrentAtom = OXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = ONeutronScatteringLength;
          *VolumeOfCurrentAtom = OVolume;
          *WeightOfCurrentAtom = OWeight;
         AtomRecg = 1;
            }
        if (AtomName[0] == ' ' && AtomName[1] == 'P' ){
          *XRayScatteringLengthOfCurrentAtom = PXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = PNeutronScatteringLength;
          *VolumeOfCurrentAtom = PVolume;
          *WeightOfCurrentAtom = PWeight;
          AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'S' ){
          *XRayScatteringLengthOfCurrentAtom = SXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          *VolumeOfCurrentAtom = SVolume;
          *WeightOfCurrentAtom = SWeight;
 	  AtomRecg = 1;
        }
        if (AtomName[0] == ' ' && AtomName[1] == 'Q' ){
          *XRayScatteringLengthOfCurrentAtom = 4.092*H2OXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = SNeutronScatteringLength;
          *VolumeOfCurrentAtom = 4.092*(H2OVolume);
	  *WeightOfCurrentAtom = 4.092*(2*HWeight+OWeight);
          AtomRecg = 1;
        }
        if (AtomName[0] == 'Z' && AtomName[1] == 'N' ){
          *XRayScatteringLengthOfCurrentAtom = ZNXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = ZNNeutronScatteringLength;
          *VolumeOfCurrentAtom = ZNVolume;
          *WeightOfCurrentAtom = ZNWeight;
           AtomRecg = 1;
        }

        if (AtomName[0] == 'C' && AtomName[1] == 'L' ){
          *XRayScatteringLengthOfCurrentAtom = CLXRayScatteringLength;
          *NeutronScatteringLengthOfCurrentAtom = CLNeutronScatteringLength;
          *VolumeOfCurrentAtom = CLVolume;
          *WeightOfCurrentAtom = CLWeight;
	   AtomRecg = 1;
        }
        if (AtomRecg == 0){
         printf("Atom name %s not found in database\n", AtomName);
        }
        AtomRecg = 0;
return;
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
    char AtomName[2];
    int DummyResidueID;
    int IDOfCurrentAtom = 0;
    int LinesToPrint = 10;
    int CountLines = 0 ;
    // I/O
    PointerToFile = fopen(Filename, "r");
printf("Printing first 10 atoms read\nType, x,y,z, ResidueNo\n");
    while (fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) {
	//                     "ATOM%*9c%c%*8c%d%*4c%lf%lf%lf%*22c%2c", &Dummychar, &ResidueID,      &xDummy, &yDummy, &zDummy, &AtomName
        if (sscanf(Linebuffer, "ATOM%*9c%c%*8c%d%*4c%lf%lf%lf%*22c%2c", &Dummychar, &DummyResidueID, &xDummy, &yDummy, &zDummy, &AtomName[0]) == 6) {
            
	    ProteinStruct.Atoms[IDOfCurrentAtom].x    = xDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].y    = yDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].z    = zDummy;
            //ProteinStruct.Atoms[IDOfCurrentAtom].Name = AtomName;
	    AssignAtom(AtomName, &ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength,&ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength,&ProteinStruct.Atoms[IDOfCurrentAtom].Volume,&ProteinStruct.Atoms[IDOfCurrentAtom].Weight);
        if ( CountLines < LinesToPrint) { // Print first 10 atom read. can be changed by LinesToPrint
	    printf("%c%c %lf %lf %lf %d \n", AtomName[0], AtomName[1], xDummy, yDummy, zDummy,DummyResidueID );
	    ++CountLines  ;
	}
            ++IDOfCurrentAtom;
        }
	//                     "HETATM%d%*3c%c%*8c%*4c%lf%lf%lf%*22c%2c", &ResidueID,      &Dummychar, &xDummy, &yDummy, &zDummy, &AtomName
        if (sscanf(Linebuffer, "HETATM%d%*3c%c%*8c%*4c%lf%lf%lf%*22c%2c", &DummyResidueID, &Dummychar, &xDummy, &yDummy, &zDummy, &AtomName[0]) == 6) {

	      
	    ProteinStruct.Atoms[IDOfCurrentAtom].x    = xDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].y    = yDummy;
            ProteinStruct.Atoms[IDOfCurrentAtom].z    = zDummy;
            //ProteinStruct.Atoms[IDOfCurrentAtom].Name = AtomName;


           
          AssignAtom(AtomName, &ProteinStruct.Atoms[IDOfCurrentAtom].XRayScatteringLength,&ProteinStruct.Atoms[IDOfCurrentAtom].NeutronScatteringLength,&ProteinStruct.Atoms[IDOfCurrentAtom].Volume,&ProteinStruct.Atoms[IDOfCurrentAtom].Weight);
        if ( CountLines < LinesToPrint) {// Print first 10 atom read. can be changed by LinesToPrint
	     printf("%c%c %lf %lf %lf %d \n", AtomName[0], AtomName[1], xDummy, yDummy, zDummy,DummyResidueID );

++CountLines ; 
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
    FILE *PointerToDummyPDBFile;
    char Linebuffer[256];
    char AtomName[2];
    char ResidueName[3];
    char Dummychar;
    int PreviousResidueID = 0;
    int ResidueID = 0;
    int IDOfCurrentResidue = 0;
	
    double VolumeOfResidue = 0.0;
    double WeightOfResidue = 0.0;
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
    double WeightOfCurrentAtom;
    double xDummy;
    double yDummy;
    double zDummy;
	    

    int LinesToPrint = 10;
    int CountLines = 0;
    // I/O
    PointerToFile = fopen(Filename, "r");

    // Reading of Amino acids:
    // We read each line that starts with either "ATOM" or "HETATM"
    // For the "ATOM" lines residues are assigned based on ResidueID.
    // Each residue has a weight, a volume and a scattering length density (in xrays and neutrons)
    // As the lines are read these properties are updated according the atom on the current line until the ResidueID changes
   
    PointerToDummyPDBFile = fopen("LinesNotRead.txt", "w+");
    while (fgets(Linebuffer, sizeof(Linebuffer), PointerToFile) != NULL) { //Read File, defined above, line by line and store it in "Linebuffer"
        ResidueID = 0;

        // **** These if statements are for testing what gets read
	//                     "ATOM%*9c%c%*8c%d%*4c%lf%lf%lf%*22c%2c", &Dummychar, &ResidueID,      &xDummy, &yDummy, &zDummy, &AtomName
	//if (sscanf(Linebuffer, "ATOM%*9c%c%*3c%3c%*2c%d%*c%lf%lf%lf%*22c%2c", &Dummychar, &ResidueName,&ResidueID, &xDummy, &yDummy, &zDummy, &AtomName) == 7) { //This Line for reading amino acids
        if (sscanf(Linebuffer, "ATOM%5d%*3c%c%*3c%3c%*6c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &ResidueName[0], &xDummy, &yDummy, &zDummy, &AtomName[0]) == 7) { // This line for reading individual atoms
	//Check that we read
        //printf("%s",Linebuffer);

	
	}

	else if (sscanf(Linebuffer, "HETATM%5d%*3c%c%*3c%3c%*6c%lf%lf%lf%*22c%2c", &ResidueID,&Dummychar,  &ResidueName[0], &xDummy, &yDummy, &zDummy, &AtomName[0]) == 7) {
		//Check that we read
        //printf("%s",Linebuffer);
        //printf("%i %s %c %i\n",ResidueID,AtomName, Dummychar,IDOfCurrentResidue);
	}
	else {
              fprintf(PointerToDummyPDBFile,"%s", Linebuffer);

	}
	
 	// ****
	
        if (ResidueID != PreviousResidueID) {
        	if (CountLines < LinesToPrint) {
        printf("%i, %s\n",ResidueID, ResidueName);
        ++CountLines;
	}
	    	if (PreviousResidueID != 0) {

                ProteinStruct.Residues[IDOfCurrentResidue].Volume = VolumeOfResidue;
		ProteinStruct.Residues[IDOfCurrentResidue].Weight = WeightOfResidue;

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
		WeightOfResidue = 0.0;
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
            ProteinStruct.Residues[IDOfCurrentResidue].Weight = WeightOfResidue;

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
//
        AssignAtom(AtomName, &XRayScatteringLengthOfCurrentAtom,&NeutronScatteringLengthOfCurrentAtom,&VolumeOfCurrentAtom, &WeightOfCurrentAtom);
      }

        if (sscanf(Linebuffer, "HETATM%*11c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name)== 1) {
//        printf("hetatm %i ",sscanf(Linebuffer, "HETATM%*13c%3c", ProteinStruct.Residues[IDOfCurrentResidue].Name) );
          if (strcmp(ProteinStruct.Residues[IDOfCurrentResidue].Name, "UZ9") == 0 ){ //Does Residue name match "UZ9"
              strcpy(ProteinStruct.Residues[IDOfCurrentResidue].Name, "X"); //Change name to "X" 

	  }
	 AssignAtom(AtomName, &XRayScatteringLengthOfCurrentAtom,&NeutronScatteringLengthOfCurrentAtom,&VolumeOfCurrentAtom, &WeightOfCurrentAtom);
  
      }

        VolumeOfResidue += VolumeOfCurrentAtom;
        WeightOfResidue += WeightOfCurrentAtom;
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
    fclose(PointerToDummyPDBFile);

    fclose(PointerToFile);
}
