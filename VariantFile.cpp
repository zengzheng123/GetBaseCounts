//
//  VariantFile.cpp
//  GetBaseCounts
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 7/31/15.
//  Copyright (c) 2015 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//

#include "VariantFile.h"


VcfFile::VcfFile(string input_vcf_file)
{
    if((input_vcf_file.length() > 3 && input_vcf_file.substr(input_vcf_file.length() - 3) == ".gz") || (input_vcf_file.length() > 5 && input_vcf_file.substr(input_vcf_file.length() - 5) == ".gzip"))
    {
        cout << "[INFO] Input vcf is gzip file" << endl;
        vcf_gz_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_gz_fs;
    }
    else if(input_vcf_file.length() > 4 && input_vcf_file.substr(input_vcf_file.length() - 4) == ".vcf")
    {
        cout << "[INFO] Input vcf is plain text file" << endl;
        vcf_txt_fs.open(input_vcf_file.c_str());
        vcf_fs = &vcf_txt_fs;
    }
    else
    {
        cerr << "[ERROR] Input vcf is in unrecognized format" << endl;
        exit(1);
    }
    if(!(*vcf_fs))
    {
        cerr << "[ERROR] Fail to open input target file: " << input_vcf_file << endl;
        exit(1);
    }
}

VcfFile::~VcfFile()
{
    close();
}

bool VcfFile::get_next(string &line)
{
    if(cur_line.empty())
    {
        if(!getline((*vcf_fs), line))
            return false;
    }
    else
    {
        line = cur_line;
        cur_line.clear();
    }
    return true;
}

void VcfFile::roll_back(string &line)
{
    cur_line = line;
}

void VcfFile::close()
{
    if(vcf_txt_fs)
        vcf_txt_fs.close();
    if(vcf_gz_fs)
        vcf_gz_fs.close();
}

bool VcfFile::eof()
{
    return vcf_fs->eof();
}

