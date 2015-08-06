//
//  VariantEntry.h
//  GetBaseCounts
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#ifndef __GetBaseCounts__VariantEntry__
#define __GetBaseCounts__VariantEntry__

#include <stdio.h>
#include <string>
using namespace std;

class VcfEntry
{
public:
    
    VcfEntry() {}
    
    ~VcfEntry() {}
    
    string chrom;
    int pos;
    string id;
    string ref;
    string alt;
    string output;
};

#endif /* defined(__GetBaseCounts__VariantEntry__) */
