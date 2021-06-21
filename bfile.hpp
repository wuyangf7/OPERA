//
//  bfile.hpp
//  SMR_CPP
//
//  Created by Futao Zhang on 5/07/2018.
//  Copyright © 2018 Futao Zhang. All rights reserved.
//

#ifndef bfile_hpp
#define bfile_hpp

#include "SMR_data.h"

namespace SMRDATA
{
    void ld_report(char* outFileName, char* bFileName,char* indilstName, char* indilst2remove,char* snplstName, char* snplst2exclde,int chr, double maf, bool ldr, bool ldr2, int ldWind);
}

#endif /* bfile_hpp */
