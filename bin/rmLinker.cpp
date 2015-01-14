#include <iostream>
#include <string>
#include <vector>
#include <iostream> // input and output files
#include <fstream>  
#include <map>
#include <sstream>  // for to_string function
#include <unistd.h> // get options
#include <getopt.h> // get long options
#include <stdlib.h> // atoi
using namespace std;
/** Remove linker sequence from paired-end mates */

// Functions ######################################################################
//#################################################################################
void print_usage() 
{
	printf("================================================\n");
	printf("Remove linker sequence from paired-end reads.\n");
    printf("Usage: rmLinker -1 fastq1 -2 fastq2 -l forward_linker_sequence [OPTIONS]\n");
    printf("Options:\n");
    printf("\t-p --prefix\tPrefix of output files (two fastq, and one txt files). Default: 'result'.\n");
    printf("\t-m --mismatches\tNumber of mismatches allowed. Default: 2.\n");
    printf("\t-f --direction1\tComma separated list of allowed linker directions on the forward mate. Default: 'fwd,fwdC,rev,revC'.\n");
    printf("\t-r --direction2\tComma separated list of allowed linker directions on the reverse mate. Default: 'fwd,fwdC,rev,revC'.\n");
}
template <typename T>
string to_string(T value)
{
  //create an output string stream
  ostringstream os ;
  //throw the value into the string stream
  os << value ;
  //convert the string stream into a string and return
  return os.str() ;
}
//#################################################################################
int countMismatches(string const & read, string const & query)
{
    int localScore = 0;
    for (unsigned i = 0; i < query.size(); ++i)
        if (read[i] != query[i])
            ++localScore;
    return localScore;
}
//#################################################################################
struct alignment {
	int start;
	int end;
	string strand;
	string type;
};
//#################################################################################
alignment computeAlignment(map<string,string> & queries, string const & read, int const & numMismatches) 
{
//	string read = mate["seq"];
    int mismatches, minMismatches;
	int matchLength=0;
	int startR, endR, startQ, endQ;
	alignment output;
	output.start = read.size(); output.end=read.size(); output.strand="none"; output.type="none";
	string query, strand;

	minMismatches = numMismatches + 1; 
	matchLength   = 0;
	// Loop over queries
    map<string, string>::iterator iter;
    for( iter=queries.begin(); iter!=queries.end(); ++iter ){
		strand = iter->first; 
		query  = iter->second;
		// Iterate over read length
		for (unsigned i = 0; i < read.size() ; ++i){
			startR = i; // start of read
			endR   = min( read.size(), i + query.size() ); // end of read
			startQ = 0; // start of query
			endQ   = min<unsigned int>( query.size(), startQ +(endR-startR) ); // end of query
			mismatches = countMismatches( read.substr(startR, endR-startR ), query.substr(startQ,endQ-startQ) );

			if( (mismatches <= numMismatches) and (mismatches < minMismatches) and ((endR-startR)>matchLength) ){
				minMismatches = mismatches;
				matchLength   = endR - startR;
				//Requiere at least 15 nt to call a match
				if( matchLength >=15 ){
					output.start  = startR; output.end = endR; 
					output.strand = strand; output.type = "match";
				}else if (mismatches==0){ // If match happens in beyond the last 15nt, assume is only by chance and use it for trimming. Don't tolerate mismatches here.
					output.start  = startR; output.end = endR; 
					output.strand = strand; output.type = "trim";
	}	}	}	}
    return output;
}
// ##########################################################################################
struct readID{
	string read_idS; // sequence read id
	string read_idQ; // quality read id
};
readID getRead_ids(string const & read_id_original,string const & classification){
	readID ids;
	string prefix, suffix;
	stringstream read_id_stream(read_id_original.substr(1));
	getline(read_id_stream,prefix,' ');
	getline(read_id_stream,suffix,' ');
	ids.read_idS = "@" + prefix + "_type_" + classification + " " + suffix; 
	ids.read_idQ = "+" + prefix + "_type_" + classification + " " + suffix;
	return ids;
}
// ##########################################################################################
int printMates(map<string,string> & queries1, map<string,string> & queries2, map<string,string> & mate1, map<string,string> & mate2, int const & numMismatches, ofstream & detailsFile, ofstream & mate1File, ofstream & mate2File)
{
	// Compute alignments
	alignment align1 = computeAlignment( queries1, mate1["seq"], numMismatches);
	alignment align2 = computeAlignment( queries2, mate2["seq"], numMismatches);
	// Classify mates
	string classification = align1.strand+"_"+align2.strand;
	// Chop fastq sequences and qualities
	string seq1  = mate1["seq" ].substr( 0, align1.start );
	string qual1 = mate1["qual"].substr( 0, align1.start );
	string seq2  = mate2["seq" ].substr( 0, align2.start );
	string qual2 = mate2["qual"].substr( 0, align2.start );
	// Print details
	detailsFile << mate1["id"]+"\t"+to_string(align1.start)+"\t"+to_string(align1.end)+"\t"+align1.strand+"\t"+align1.type+"\t"+to_string(align2.start)+"\t"+to_string(align2.end)+"\t"+align2.strand+"\t"+align2.type+"\n";
	// If reads are to short, or linker are not oriented properly, don't print them.
	if ( (align1.start>=15) and (align2.start>=15)  ){
		readID id1 = getRead_ids( mate1["id"], classification);
		readID id2 = getRead_ids( mate2["id"], classification);
		mate1File << id1.read_idS +"\n"+seq1+"\n"+id1.read_idQ+"\n"+qual1+"\n";
		mate2File << id2.read_idS +"\n"+seq2+"\n"+id2.read_idQ+"\n"+qual2+"\n";
	}
	return 0;
}
// ##########################################################################################
string getComplement( string const & forward )
{
	string complement;
    for( string::size_type i = 0; i < forward.size(); ++i ) {
        if      ( (forward[i] == 'A') or (forward[i] == 'a') ) complement += 'T';
        else if ( (forward[i] == 'C') or (forward[i] == 'c') ) complement += 'G';
        else if ( (forward[i] == 'G') or (forward[i] == 'g') ) complement += 'C';
        else if ( (forward[i] == 'T') or (forward[i] == 't') ) complement += 'A';
    }	
	return complement;
}
// ##########################################################################################
string getReverse( string const & forward )
{
    string reversed="";
    for( unsigned i=0; i<forward.length(); ++i) {
        reversed = forward[i]+reversed ;
    }
    return reversed;
}
// Execution ################################################################################
//###########################################################################################
int main(int argc, char *argv[])
{
	// Get arguments ########################################################################
	string fastq1="none", fastq2="none", linker="none";
	string direction1="fwd,fwdC,rev,revC";
	string direction2="fwd,fwdC,rev,revC";
	string prefix = "result";
	string help="h";
	int mismatches=2;
    //Specifying the expected options
    static struct option long_options[] = {
        {"help",      no_argument,       0,  'h' },
        {"fastq1",    required_argument, 0,  '1' },
        {"fastq2",    required_argument, 0,  '2' },
        {"linker",    required_argument, 0,  'l' },
        {"prefix",    required_argument, 0,  'p' },
        {"mismatches",required_argument, 0,  'm' },
        {"direction1",required_argument, 0,  'f' },
        {"direction2",required_argument, 0,  'r' },
        {0,           0,                 0,  0   }
    };
    int opt = 0, long_index = 0;
    while ( (opt = getopt_long(argc, argv,"h1:2:l:p:m:d:f:r:", long_options, &long_index))!=-1 ){
		switch (opt) {
             case 'h' : print_usage();
			 	exit(EXIT_FAILURE);
             case '1' : fastq1 = optarg;
				 break;
             case '2' : fastq2 = optarg;
                 break;
             case 'l' : linker = optarg; 
                 break;
             case 'p' : prefix = optarg;
                 break;
             case 'm' : mismatches = atoi(optarg);
                 break;
             case 'f' : direction1 = optarg;
                 break;
             case 'r' : direction2 = optarg;
                 break;
             default: exit(EXIT_FAILURE);
        }
    }
	if ( (fastq1=="none") or (fastq2=="none") or (linker=="none") ) {
		printf("Missing one of the following options -1, -2, and/or -l.\n");
		print_usage();
		exit(EXIT_FAILURE);
	}
	// Open input and output files###########################################################
	string detailsFileName = prefix+"_details.txt";
	string mate1FileName   = prefix+"_mate1.fastq";
	string mate2FileName   = prefix+"_mate2.fastq";
	ifstream fastqFile1  (fastq1.c_str()        ); 
	ifstream fastqFile2  (fastq2.c_str()         );
	ofstream detailsFile (detailsFileName.c_str());
	ofstream mate1File   (mate1FileName.c_str()  ); 
	ofstream mate2File   (mate2FileName.c_str()  );
	// Define all possible query directions ##############################################
	map<string,string> sequences, queries1, queries2;
	sequences[ "fwd" ] = linker;                             // forward
	sequences[ "rev" ] = getReverse( sequences["fwd"]  );    // reverse
	sequences[ "fwdC" ] = getComplement( sequences["fwd"] ); // forward complement
	sequences[ "revC" ] = getComplement( sequences["rev"] ); // reverse complement
	stringstream dir1_stream(direction1), dir2_stream(direction2);
	string orientation;
	while( getline( dir1_stream,orientation,',') ) // add linker orientations specified by the user
		queries1[ orientation ] = sequences[ orientation ];
	while( getline( dir2_stream,orientation,',') ) // add liner orientations specified by the user
		queries2[ orientation ] = sequences[ orientation ];
	// Loop over fastq read mates ##########################################################
	cout << "Number of processed reads" <<endl; // Report progress
	int count = 0;
	cout << count << endl;
	string line1, line2;
	map<string,string> mate1, mate2;
	vector<string> keys;
	keys.push_back("id");keys.push_back("seq");keys.push_back("id2");keys.push_back("qual");
	while ( ! fastqFile1.eof() ) { 
		++count;
		if( count%1000000==0) cout << count << endl; // Report progress
		// Read four lines at a time
		for( unsigned i=0; i<4; ++i ){
			getline(fastqFile1,line1);	
			getline(fastqFile2,line2);
			mate1[keys[i]] = line1;
			mate2[keys[i]] = line2;	
		}
		printMates( queries1,queries2, mate1, mate2, mismatches,detailsFile,mate1File,mate2File);
	}
	// Close files #########################################################################	
	fastqFile1.close();
	fastqFile2.close();
	detailsFile.close();
	mate1File.close();
	mate2File.close();
    
	return 0;
}
