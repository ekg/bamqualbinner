#include <iostream>
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <map>
#include <vector>

#include "api/BamReader.h"
#include "api/BamWriter.h"
#include "api/BamAlignment.h"


using namespace std;
using namespace BamTools;


short qualityChar2ShortInt(char c) {
    return static_cast<short>(c) - 33;
}

char shortInt2QualityChar(short i) {
    return static_cast<char>(i + 33);
}

void countMismatchesAndGaps(
    BamAlignment& alignment,
    vector<CigarOp>& cigarData,
    string referenceSequence,
    int& mismatches,
    int& gaps,
    int& gapslen,
    int& softclips,
    int& mismatchQsum,
    int& softclipQsum
    ) {

    int sp = 0;
    int rp = 0;
    for (vector<CigarOp>::const_iterator c = cigarData.begin();
	 c != cigarData.end(); ++c) {
        int l = c->Length;
        char t = c->Type;
        if (t == 'M') { // match or mismatch
            for (int i = 0; i < l; ++i) {
                if (alignment.QueryBases.at(rp) != referenceSequence.at(sp)) {
                    ++mismatches;
		    mismatchQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		}
                ++sp;
                ++rp;
            }
        } else if (t == 'D') { // deletion
            ++gaps;
	    gapslen += l;
            sp += l;  // update reference sequence position
        } else if (t == 'I') { // insertion
	    ++gaps;
	    gapslen += l;
	    rp += l;  // update read position
	} else if (t == 'S') { // soft clip, clipped sequence present in the read not matching the reference
	    softclips += l;
	    for (int i = 0; i < l; ++i) {
		softclipQsum += qualityChar2ShortInt(alignment.Qualities.at(rp));
		++rp;
	    }
	} else if (t == 'H') { // hard clip on the read, clipped sequence is not present in the read
	} else if (t == 'N') { // skipped region in the reference not present in read, aka splice
	    sp += l;
	}
    }

}


inline void binQuals(BamAlignment& alignment, const vector<char>& bins) {
    for (string::iterator q = alignment.Qualities.begin(); q != alignment.Qualities.end(); ++q) {
        *q = bins.at(qualityChar2ShortInt(*q));
    }
}


void printUsage(char** argv) {
    cerr << "usage: [BAM data stream] | " << argv[0] << " [options]" << endl
         << endl
         << "author: Erik Garrison <erik.garrison@bc.edu>" << endl
         << endl
         << "Bins quality scores in alignments according to the Illumina quality binning scheme." << endl
         << "Old Quality Score range -> New quality score:" << endl
         << endl
         << "N (no call) -> N" << endl
         << "2–9 -> 6" << endl
         << "10–19 -> 15" << endl
         << "20–24 -> 22" << endl
         << "25–29 -> 27" << endl
         << "30–34 -> 33" << endl
         << "35–39 -> 37" << endl
         << "40 -> 40" << endl
         << endl
         << "arguments:" << endl
         << "    -d --debug                 Print debugging information about realignment process" << endl
         << "    -s --suppress-output       Don't output BAM on stdout" << endl;
}

int main(int argc, char** argv) {

    int c;

    bool suppress_output = false;
    bool debug = false;
    
    if (argc < 2) {
        printUsage(argv);
        exit(1);
    }

    while (true) {
        static struct option long_options[] =
        {
            {"help", no_argument, 0, 'h'},
            {"debug", no_argument, 0, 'd'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hd",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
 
        switch (c) {

            case 'd':
                debug = true;
                break;

            case 'h':
                printUsage(argv);
                exit(0);
                break;
              
            case '?':
                printUsage(argv);
                exit(1);
                break;
     
              default:
                abort();
                break;
        }
    }

    BamReader reader;
    if (!reader.Open("stdin")) {
        cerr << "could not open stdin for reading" << endl;
        exit(1);
    }

    BamWriter writer;
    if (!suppress_output && !writer.Open("stdout", reader.GetHeaderText(), reader.GetReferenceData())) {
        cerr << "could not open stdout for writing" << endl;
        exit(1);
    }

    // store the names of all the reference sequences in the BAM file
    map<int, string> referenceIDToName;
    vector<RefData> referenceSequences = reader.GetReferenceData();
    int i = 0;
    for (RefVector::iterator r = referenceSequences.begin(); r != referenceSequences.end(); ++r) {
        referenceIDToName[i] = r->RefName;
        ++i;
    }

    BamAlignment alignment;
    map<long unsigned int, vector<BamAlignment> > alignmentSortQueue;

    vector<char> bins;
    i = 0;

    for (; i < 2; ++i)
        bins.push_back(shortInt2QualityChar(0));
    for (; i < 10; ++i)
        bins.push_back(shortInt2QualityChar(6));
    for (; i < 20; ++i)
        bins.push_back(shortInt2QualityChar(15));
    for (; i < 25; ++i)
        bins.push_back(shortInt2QualityChar(22));
    for (; i < 30; ++i)
        bins.push_back(shortInt2QualityChar(27));
    for (; i < 35; ++i)
        bins.push_back(shortInt2QualityChar(33));
    for (; i < 40; ++i)
        bins.push_back(shortInt2QualityChar(37));
    for (; i < 90; ++i)
        bins.push_back(shortInt2QualityChar(40));


    while (reader.GetNextAlignment(alignment)) {
        binQuals(alignment, bins);
        if (!suppress_output)
            writer.SaveAlignment(alignment);
    }

    reader.Close();
    if (!suppress_output)
        writer.Close();

    return 0;

}
