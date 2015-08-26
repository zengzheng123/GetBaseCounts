//
//  GetBaseCounts.cpp
//  GetBaseCounts
//
//  Created by Zeng, Zheng/Sloan-Kettering Institute on 9/5/14.
//  Copyright (c) 2014 Zeng, Zheng/Sloan-Kettering Institute. All rights reserved.
//


#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string>
#include <stdio.h>
#include <sstream>
#include <fstream>
#include <getopt.h>
#include <vector>
#include <algorithm>
#include "api/BamReader.h"
#include "gzstream.h"
#include "omp.h"
#include "VariantFile.h"
#include "VariantEntry.h"
#include "Global.h"

//#define _DEBUG
//#define _PARSING_DEBUG
//#define _FASTA_DEBUG

using namespace std;
using namespace BamTools;

const string VERSION = "GetBaseCounts 1.4.0";

const char bases[] = {'A', 'C', 'G', 'T', 'a', 'c', 'g', 't', '+', '-'};
const size_t base_size = 10;

map<char, size_t> base_index;
string input_fasta_file;
string input_bam_file;
string input_vcf_file;
string output_file;
int mapping_quality_threshold = 15;
int base_quality_threshold = 20;
int minimum_coverage_threshold = 0;
int quality_scale = 33;
int filter_duplicate = 1;
int filter_improper_pair = 1;
int filter_qc_failed = 0;
int filter_indel = 0;
int filter_non_primary = 1;
int maximum_vcf_block_size = 10000;
int maximum_vcf_block_distance = 1000000;
int num_thread = 1;
bool sort_output = false;
bool compress_output = false;
int max_warning_per_type = 3;
int warning_vcf_unsorted = 0;
int warning_inconsistent_ref_allele = 0;




void printUsage(string msg = "")
{
    cout << endl;
    cout << VERSION << endl;
    cout << "Usage: " << endl;
    cout << "[REQUIRED ARGUMENTS]" << endl;
    cout << "\t--fasta                 <string>                        Input reference sequence file" << endl;
    cout << "\t--bam                   <string>                        Input bam file" << endl;
    cout << "\t--vcf                   <string>                        Input vcf file, it needs to be sorted within each chromosome to optimize the running speed, gzipped vcf file is supported" << endl;
    cout << "\t--output                <string>                        Output file" << endl;
    cout << endl;
    cout << "[OPTIONAL ARGUMENTS]" << endl;
    cout << "\t--thread                <int>                           Number of thread. Default " << num_thread << endl;
    cout << "\t--sort_output                                           Sort output file by genomic position, this option requires addtional memory" << endl;
    cout << "\t--compress_output                                       Compress the output and write gzipped file directly" << endl;
    cout << "\t--maq                   <int>                           Mapping quality threshold. Default " << mapping_quality_threshold << endl;
    cout << "\t--baq                   <int>                           Base quality threshold, Default " << base_quality_threshold << endl;
    cout << "\t--cov                   <int>                           Minimum coverage applied to BASEQ_depth. Default " << minimum_coverage_threshold << endl;
    cout << "\t--filter_duplicate      [0, 1]                          Whether to filter reads that are marked as duplicate. 0=off, 1=on. Default " << filter_duplicate << endl;
    cout << "\t--filter_improper_pair  [0, 1]                          Whether to filter reads that are marked as improperly paired. 0=off, 1=on. Default " <<  filter_improper_pair << endl;
    cout << "\t--filter_qc_failed      [0, 1]                          Whether to filter reads that are marked as failed quality control. 0=off, 1=on. Default " << filter_qc_failed << endl;
    cout << "\t--filter_indel          [0, 1]                          Whether to filter reads that contain indels. 0=off, 1=on. Default " << filter_indel << endl;
    cout << "\t--filter_non_primary    [0, 1]                          Whether to filter reads that are marked as non primary alignment. Default " << filter_non_primary << endl;
    cout << "\t--suppress_warning      <int>                           Only print a limit number of warnings for each type. Default " << max_warning_per_type << endl;
    cout << "\t--help                                                  Print command line usage" << endl;
    cout << endl;
    cout << "[ADVANCED ARGUMENTS, CHANGING THESE ARGUMENTS MAY SIGNIFICANTLY AFFECT MEMORY USAGE AND RUNNING TIME. USE WITH CAUTION]" << endl;
    cout << "\t--max_block_size        <int>                           The maximum number of vcf entries that can be processed at once per thread. Default " << maximum_vcf_block_size << endl;
    cout << "\t--max_block_dist        <int>                           The longest spanning region (bp) of vcf chunks that can be processed at once per thread. Default " << maximum_vcf_block_distance << endl;
    cout << endl;
    if(!msg.empty())
        cerr << msg << endl;
    exit(1);
}


static struct option long_options[] =
{
    {"fasta",                   required_argument,      0,     'f'},
    {"bam",                     required_argument,      0,     'b'},
    {"vcf",                     required_argument,      0,     'v'},
    {"output",                  required_argument,      0,     'o'},
    {"sort_output",             no_argument,            0,     's'},
    {"compress_output",         no_argument,            0,     'C'},
    {"thread",                  required_argument,      0,     't'},
    {"maq",                     required_argument,      0,     'Q'},
    {"baq",                     required_argument,      0,     'q'},
    {"cov",                     required_argument,      0,     'c'},
    {"filter_duplicate",        required_argument,      0,     'd'},
    {"filter_improper_pair",    required_argument,      0,     'p'},
    {"filter_qc_failed",        required_argument,      0,     'l'},
    {"filter_indel",            required_argument,      0,     'i'},
    {"filter_non_primary",      required_argument,      0,     'n'},
    {"suppress_warning",        required_argument,      0,     'w'},
    {"max_block_size",          required_argument,      0,     'M'},
    {"max_block_dist",          required_argument,      0,     'm'},
    {"help",                    no_argument,            0,     'h'},
    {0, 0, 0, 0}
};


void parseOption(int argc, const char* argv[])
{
    if(argc == 1)
        printUsage();
    int next_option;
    int option_index = 0;
    do
    {
        next_option = getopt_long(argc, const_cast<char**>(argv), "f:b:v:o:sCt:Q:q:c:d:p:l:i:n:w:M:m:h", long_options, &option_index);
        switch(next_option)
        {
            case 'f':
                input_fasta_file = optarg;
                break;
            case 'b':
                input_bam_file = optarg;
                break;
            case 'v':
                input_vcf_file = optarg;
                break;
            case 'o':
                output_file = optarg;
                break;
            case 's':
                sort_output = true;
                break;
            case 'C':
                compress_output = true;
                break;
            case 't':
                num_thread = atoi(optarg);
                break;
            case 'Q':
                mapping_quality_threshold = atoi(optarg);
                break;
            case 'q':
                base_quality_threshold = atoi(optarg);
                break;
            case 'c':
                minimum_coverage_threshold = atoi(optarg);
                break;
            case 'd':
                filter_duplicate = atoi(optarg);
                break;
            case 'p':
                filter_improper_pair = atoi(optarg);
                break;
            case 'l':
                filter_qc_failed = atoi(optarg);
                break;
            case 'i':
                filter_indel = atoi(optarg);
                break;
            case 'n':
                filter_non_primary = atoi(optarg);
                break;
            case 'w':
                max_warning_per_type = atoi(optarg);
                break;
            case 'M':
                maximum_vcf_block_size = atoi(optarg);
                break;
            case 'm':
                maximum_vcf_block_distance = atoi(optarg);
                break;
            case 'h':
                printUsage();
                break;
            case -1:
                break; //parsed all options
            default:
                printUsage(string("[ERROR] Argument error: ") + argv[optind - 1]);
        }
    } while(next_option != -1);

    if(input_fasta_file.empty())
        printUsage("[ERROR] Please specify input fasta file");
    if(input_bam_file.empty())
        printUsage("[ERROR] Please specify input bam file");
    if(input_vcf_file.empty())
        printUsage("[ERROR] Please specify input vcf file");
    if(output_file.empty())
        printUsage("[ERROR] Please specify output file");
    if(num_thread <= 0)
        printUsage("[ERROR] Invalid number of threads");
    if(maximum_vcf_block_size <= 0)
        printUsage("[ERROR] Invalid max_block_size");
    if(maximum_vcf_block_distance <= 0)
        printUsage("[ERROR] Invalid max_block_dist");
    if(filter_duplicate != 0 && filter_duplicate != 1)
        printUsage("[ERROR] --filter_duplicate should be 0 or 1");
    if(filter_improper_pair != 0 && filter_improper_pair != 1)
        printUsage("[ERROR] --filter_improper_pair should be 0 or 1");
    if(filter_qc_failed != 0 && filter_qc_failed != 1)
        printUsage("[ERROR] --filter_qc_failed should be 0 or 1");
    if(filter_indel != 0 && filter_indel != 1)
        printUsage("[ERROR] --filter_indel should be 0 or 1");
    if(filter_non_primary != 0 && filter_non_primary != 1)
        printUsage("[ERROR] --filter_non_primary should be 0 or 1");
    base_quality_threshold += quality_scale;
    if(compress_output && (output_file.length() < 3 || output_file.substr(output_file.length() - 3) != ".gz")) // fix output file extension
    {
        output_file.append(".gz");
    }
#ifdef _DEBUG
    cout << "[DEBUG] Parsing options complete." << endl;
#endif
}

void build_base_index()
{
    for(size_t i = 0; i < base_size; i++)
        base_index[bases[i]] = i + 1;
    base_index['N'] = base_size;
}

void reset_hash_entry(map<char, int>& base_count)
{
    for(size_t i = 0; i < base_size; i++)
        base_count[bases[i]] = 0;
}

bool sort_vcf_by_pos(const VcfEntry& lhs, const VcfEntry& rhs)
{
    if(lhs.chrom != rhs.chrom)
    {
        string chrom1 = lhs.chrom;
        string chrom2 = rhs.chrom;
        if(chrom1.substr(0, 3) == "chr")
            chrom1 = chrom1.substr(3);
        if(chrom2.substr(0, 3) == "chr")
            chrom2 = chrom2.substr(3);
        bool chrom1_num = is_number(chrom1);
        bool chrom2_num = is_number(chrom2);
        if(chrom1_num && !chrom2_num)
            return true;
        if(!chrom1_num && chrom2_num)
            return false;
        if(chrom1_num && chrom2_num)
            return atoi(chrom1.c_str()) < atoi(chrom2.c_str());
        else //both string
            return chrom1 < chrom2;
    }
    else
        return lhs.pos < rhs.pos;
}


bool AlignmentHasIndel(BamAlignment& my_bam_alignment)
{
    vector<CigarOp>::const_iterator cigarIter;
    for(cigarIter = my_bam_alignment.CigarData.begin() ; cigarIter != my_bam_alignment.CigarData.end(); cigarIter++)
    {
        if(cigarIter->Type == 'I' || cigarIter->Type == 'D')
            return true;
    }
    return false;
}


void GetBaseCounts()
{
    vector<VcfEntry> output_sort_vec;
    build_base_index();
    
    // open VCF file
    VcfFile my_vcf(input_vcf_file);
    
    // open output file
    ostream *output_fs;
    ogzstream output_gz_fs;
    ofstream output_txt_fs;
    if(compress_output)
    {
        output_gz_fs.open(output_file.c_str());
        output_fs = &output_gz_fs;
    }
    else
    {
        output_txt_fs.open(output_file.c_str());
        output_fs = &output_txt_fs;
    }
    if(!(*output_fs))
    {
        cerr << "[ERROR] Fail to open output file: " << output_file << endl;
        exit(1);
    }
    (*output_fs) << "Chrom\tPos\tRef\tAlt\tRefidx\tTOTAL_depth\tMAPQ_depth\tBASEQ_depth\tA\tC\tG\tT\ta\tc\tg\tt\tINS\tDEL\tID" << endl;
    
    // load reference sequence
    map<string, string> reference_sequence;
    load_reference_sequence_speedup(input_fasta_file, reference_sequence, num_thread);
    
    cout << "[INFO] Processing bam file: " << input_bam_file << endl;
#pragma omp parallel num_threads(num_thread)
    {
        //int thread_num = omp_get_thread_num();
        
        // open BAM file
        BamReader my_bam_reader;
        if(!my_bam_reader.Open(input_bam_file))
        {
#pragma omp critical(output_stderr)
            {
                cerr << "[ERROR] Fail to open input bam file: " << input_bam_file << endl;
          	}
            exit(1);
        }
//         string input_bam_index_file1 = input_bam_file.substr(0, input_bam_file.length() - 3) + "bai";
//         string input_bam_index_file2 = input_bam_file + ".bai";
//         if(!my_bam_reader.OpenIndex(input_bam_index_file1))
//         {
//             if(!my_bam_reader.OpenIndex(input_bam_index_file2))
//             {
// #pragma omp critical(output_stderr)
//                 {
//                 	cerr << "[ERROR] Fail to open input bam index file: " << input_bam_index_file1 << ", or " << input_bam_index_file2 << endl;
//                 }
//                 exit(1);
//             }
//         }
        
        vector<VcfEntry> vcf_block;
        map<char, int> base_count;
        while(!my_vcf.eof())
        {
#pragma omp critical(read_vcf)
            {// read a chunk of vcf entries
                string line;
                while(my_vcf.get_next(line))
                {
                    if(line[0] == '#')
                    	continue;
                    VcfEntry new_vcf_entry;
                    
                    //split version
                    /*vector<string> vcf_item;
                    split(line, '\t', vcf_item);
                    new_vcf_entry.chrom = vcf_item[0];
                    new_vcf_entry.pos = atoi(vcf_item[1].c_str()) - 1; //convert to 0-based
                    new_vcf_entry.ref = vcf_item[3];
                    new_vcf_entry.alt = vcf_item[4];
                    new_vcf_entry.id = vcf_item[2];
                     */
                    
                    //sscanf version
                    char chrom_tmp[5000], ref_tmp[5000], alt_tmp[5000], id_tmp[5000];
                    sscanf(line.c_str(), "%s\t%d\t%s\t%s\t%s", chrom_tmp, &(new_vcf_entry.pos), id_tmp, ref_tmp, alt_tmp);
                    new_vcf_entry.chrom = chrom_tmp;
                    new_vcf_entry.pos--;
                    new_vcf_entry.id = id_tmp;
                    new_vcf_entry.ref = ref_tmp;
                    new_vcf_entry.alt = alt_tmp;
                    
                    if(((int)vcf_block.size() >= maximum_vcf_block_size) || ((int)vcf_block.size() != 0 && (new_vcf_entry.chrom != vcf_block[0].chrom || new_vcf_entry.pos - vcf_block[0].pos < 0 || new_vcf_entry.pos - vcf_block[0].pos > maximum_vcf_block_distance)))
                    { //reach chunk limit, save the current line back to vcf file handler
	                  	if(warning_vcf_unsorted < max_warning_per_type && new_vcf_entry.chrom == vcf_block[0].chrom && new_vcf_entry.pos - vcf_block[0].pos < 0)
                        {
                            cout << "[WARNING] The vcf file is not sorted within each chromosome, which might slow down the running performace significantly" << endl;
                            warning_vcf_unsorted++;
                        }
                        my_vcf.roll_back(line);
	                  	break;
                    }
                    else
                    {
	                  	vcf_block.push_back(new_vcf_entry);
                    }
              	}
            }
            
            if(vcf_block.size() > 0)
            {
                /*#pragma omp critical(output_debug)
                 {
                 cout << "thread " << thread_num << " setting block(" << vcf_block.size() << "):\t" << vcf_block[0].chrom << "\t" << vcf_block[0].pos - 1 << "\t" << vcf_block[vcf_block.size() - 1].chrom << "\t" << vcf_block[vcf_block.size() - 1].pos - 1 << endl;
                 }*/
                int refid1 = my_bam_reader.GetReferenceID(vcf_block[0].chrom);
                int refid2 = my_bam_reader.GetReferenceID(vcf_block[vcf_block.size() - 1].chrom);
                if(refid1 == -1)
                {
                    cerr << "[ERROR] Could not find vcf chrom: " << vcf_block[0].chrom << " in the bam file" << endl;
                    exit(1);
                }
                if(refid2 == -1)
                {
                    cerr << "[ERROR] Could not find vcf chrom: " << vcf_block[vcf_block.size() - 1].chrom << " in the bam file" << endl;
                    exit(1);
                }
                
                // load bam alignments for vcf chunk
                my_bam_reader.SetRegion(refid1, vcf_block[0].pos, refid2, vcf_block[vcf_block.size() - 1].pos + 1);
                vector<BamAlignment> bam_vec;
                BamAlignment new_bam_alignment;
                while(my_bam_reader.GetNextAlignment(new_bam_alignment))
                {
                    // filter duplicates, improper pairs, QC failed and indels
                    if((filter_duplicate && new_bam_alignment.IsDuplicate()) || (filter_improper_pair && !new_bam_alignment.IsProperPair()) || (filter_indel && AlignmentHasIndel(new_bam_alignment)) || (filter_qc_failed && new_bam_alignment.IsFailedQC()) || (filter_non_primary && !(new_bam_alignment.IsPrimaryAlignment())) )
                        continue;
                    bam_vec.push_back(new_bam_alignment);
                }
                
                // get base counts
                size_t start_bam_index = 0;
                bool start_bam_index_set = false;
                for(size_t vcf_index = 0; vcf_index < vcf_block.size(); vcf_index++)
                {
	                VcfEntry& my_vcf_entry = vcf_block[vcf_index];
	                string refseq_ref;
	                if(reference_sequence.find(my_vcf_entry.chrom) != reference_sequence.end()) // get ref base on reference sequence
	                	refseq_ref = reference_sequence[my_vcf_entry.chrom].substr(my_vcf_entry.pos, my_vcf_entry.ref.length());
	                else
	                {
#pragma omp critical(output_stderr)
	                	{
	                		cerr << "[ERROR] Could not find vcf chrom: " << my_vcf_entry.chrom << " in reference sequence fasta" << endl;
	                	}
                        exit(1);
	                }
	                if(refseq_ref != my_vcf_entry.ref)
	                {
#pragma omp critical(output_stderr)
	                	{
	                		if(warning_inconsistent_ref_allele < max_warning_per_type)
                            {
                                cout << "[WARNING] Inconsistency found, refseq ref base: " << refseq_ref << ", vcf ref base: " << my_vcf_entry.ref << " at position: " << my_vcf_entry.chrom << "\t" << my_vcf_entry.pos + 1 << endl;
                                warning_inconsistent_ref_allele++;
                            }
                        }
	                }
	                reset_hash_entry(base_count);  // reset base count to 0
	                //map<string, int> event_counter;
	                start_bam_index_set = false;
                    int tdepth = 0; // number of any reads cover the site
                	int depth = 0; // number of reads satisfy MAPQ that cover the site
                	int qdepth = 0; // number of reads satisfy both MAPQ and BASEQ that cover the site
                    
                    //iterate though bam entries
	                for(size_t bam_index = start_bam_index; bam_index < bam_vec.size(); bam_index++)
	                {
	                    BamAlignment &my_bam_alignment = bam_vec[bam_index];
	                    if(my_bam_alignment.Position > my_vcf_entry.pos)
	                    	break;
	                    if(my_bam_alignment.GetEndPosition(false, true) < my_vcf_entry.pos)  // GetEndPosition(bool padding, bool closed_interval)
	                    	continue;
                        if(!start_bam_index_set) // save roll back point
	                    {
	                    	start_bam_index = bam_index;
	                    	start_bam_index_set = true;
	                    }
	                   	tdepth ++;
                        if(my_bam_alignment.MapQuality < mapping_quality_threshold)
                            continue;
                        depth ++;
	                    size_t ref_pos = my_bam_alignment.Position; // the position on the reference sequence for the next base to be processed
	                    size_t read_pos = 0; // the position on the read for the next base to be processed
	                    size_t vcf_pos = 0; // target vcf position in the read
	                    bool vcf_ins = false; // there is an insertion between current base and next base
	                    bool vcf_del = false; // there is a deletion in the current base
	                    vector<CigarOp>::const_iterator cigarIter, cigarIterNext;
                        
                        // iterate through cigar strings
                        for(cigarIter = my_bam_alignment.CigarData.begin() ; cigarIter != my_bam_alignment.CigarData.end(); cigarIter++)
                        {
                            // check if the next postition to be processed is already > the target position of vcf
                            if((int)ref_pos > my_vcf_entry.pos)
                                break;
                            const CigarOp& op = (*cigarIter);
                            switch(op.Type)
                            {
                                case 'M':
                                {
                                    if((int)(ref_pos + op.Length) > my_vcf_entry.pos) // target vcf postion fall in this cigar region
                                    {
                                        vcf_pos = my_vcf_entry.pos - ref_pos + read_pos; //calculate the index of target vcf position in the read
                                        if((int)(ref_pos + op.Length) - 1 == my_vcf_entry.pos) // current postion is the last base of a M cigar region
                                        {
                                            // check if there is an indel right next to current position
                                            cigarIterNext = cigarIter + 1;
                                            if(cigarIterNext != my_bam_alignment.CigarData.end())
                                            {
                                                const CigarOp& opNext = (*cigarIterNext);
                                                ostringstream event_buffer;
                                                if(opNext.Type == 'I') // there is insertion right after the target position
                                                {
                                                    vcf_ins = true;
                                                }
                                                else if(opNext.Type == 'D') // there is deletion right after the target position
                                                {
                                                    // do nothing for now
                                                }
                                            }
                                        }
                                    }
                                    ref_pos += op.Length;
                                    read_pos += op.Length;
                                    break;
                                }
                                case 'I':
                                {
                                    read_pos += op.Length;
                                    break;
                                }
                                case 'S' :
                                {
                                    read_pos += op.Length;
                                    break;
                                }
                                case 'D': case 'N':
                                {
                                    if(op.Type == 'D' && (int)(ref_pos + op.Length) > my_vcf_entry.pos) // target vcf postion fall in this cigar region
                                        vcf_del = true;
                                    ref_pos += op.Length;
                                    break;
                                }
                                case 'H':
                                    break;
                                default:
                                    break;
                            }
                        }
	                    if(vcf_del)
	                    {
	                    	base_count['-'] ++;
	                    }
	                    else // a read could count both for ACGTacgt and +
	                    {
	                     	if(vcf_ins)
	                    	{
	                    		base_count['+'] ++;
	                    	}
	                     	char cur_alt = my_bam_alignment.QueryBases[vcf_pos];
		                    int cur_bq = my_bam_alignment.Qualities[vcf_pos];
		                    if(cur_bq >= base_quality_threshold)
		                    {
		                        qdepth++;
		                        if(my_bam_alignment.IsReverseStrand())
		                            base_count[tolower(cur_alt)] ++;
		                        else
		                            base_count[toupper(cur_alt)] ++;
		                    }
	                  	}
                    }
                    if(qdepth >= minimum_coverage_threshold)
                    {
	                    ostringstream output_buffer;
                        output_buffer << my_vcf_entry.chrom << "\t" << (my_vcf_entry.pos + 1) << "\t" << refseq_ref << "\t" << my_vcf_entry.alt << "\t" << (my_vcf_entry.ref.length() == my_vcf_entry.alt.length() ? (int)base_index[refseq_ref[0]] : -1) << "\t" << tdepth << "\t" << depth << "\t" << qdepth;
	                    for(size_t i = 0; i < base_size; i++)
                            output_buffer << "\t" << base_count[bases[i]];

	                    output_buffer << "\t" << my_vcf_entry.id << endl;
                        if(sort_output)  // save output for sorting
                        {
                            my_vcf_entry.output = output_buffer.str();
#pragma omp critical(output_counts)
                            {
                                output_sort_vec.push_back(my_vcf_entry);
                            }
                        }
                        else  // write output to file
                        {
#pragma omp critical(output_counts)
                            {
                                (*output_fs) << output_buffer.str();
                            }
                        }
                    }
                }
            }
            vcf_block.clear();
        }
        my_bam_reader.Close();
    }
    if(sort_output)
    {
        sort(output_sort_vec.begin(), output_sort_vec.end(), sort_vcf_by_pos);
        for(size_t output_index = 0; output_index < output_sort_vec.size(); output_index++)
            (*output_fs) << output_sort_vec[output_index].output;
    }
    my_vcf.close();
    if(output_gz_fs)
        output_gz_fs.close();
    if(output_txt_fs)
        output_txt_fs.close();
    cout << "[INFO] Finished processing bam file" << endl;
}



int main(int argc, const char * argv[])
{
    parseOption(argc, argv);
    GetBaseCounts();
    return 0;
}



