#include <iostream>
#include <string>
#include <vector>
#include <iostream> // input and output files
#include <fstream>  
#include <map>
#include <sstream> // for to_string function
#include <unistd.h> // get options
#include <getopt.h> // get long options
using namespace std;

////////////////////////////////////////////////////////////////////////////////////
// Functions
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

int countMismatches(string const & read, string const & query)
{
    int localScore = 0;
    for (unsigned i = 0; i < query.size(); ++i)
        if (read[i] != query[i])
            ++localScore;
    return localScore;
}
struct alignment {
	int start;
	int end;
	string strand;
	string type;
};
alignment computeAlignment(vector<string> const & queries, string const & read, int const & numMismatches) 
{
//	string read = mate["seq"];
    int mismatches, minMismatches;
	int matchLength=0;
	int startR, endR, startQ, endQ;
	alignment output;
	output.start = read.size(); output.end=read.size(); output.strand="none"; output.type="none";
	string strand;

	minMismatches = numMismatches + 1; 
	matchLength   = 0;
	// Loop over queries
	for( unsigned j=0; j < queries.size(); ++j){
		// Determine query and strand 
		string query = queries[j];
		if      (j==0) strand = "fwd";  // forward 
		else if (j==1) strand = "fwdC"; // forward complement
		else if (j==2) strand = "rve";  // reverse 
		else if (j==3) strand = "rveC"; // reverse complement 
		// Iterate over read length
		for (unsigned i = 0; i < read.size() ; ++i){
			startR = i; // start of read
			endR   = min( read.size(), i + query.size() ); // end of read
			startQ = 0; // start of query
			endQ   = min<unsigned int>( query.size(), startQ +(endR-startR) ); // end of query
			mismatches = countMismatches( read.substr(startR, endR-startR ), query.substr(startQ,endQ-startQ) );

			if( (mismatches <= numMismatches) & (mismatches < minMismatches) & ((endR-startR)>matchLength) ){
				minMismatches = mismatches;
				matchLength   = endR - startR;
				//Requiere at least 15 nt to call a match
				if( matchLength >=15 ){
					output.start  = startR; output.end = endR; 
					output.strand = strand; output.type = "match";
				}else if (mismatches==0){ // If match happens in the last 15nt, assume is only by chance and use it for trimming, in this case don't tolerate mismatches
					output.start  = startR; output.end = endR; 
					output.strand = strand; output.type = "trim";
	}	}	}	}
    return output;
}

int printMates(vector<string> const & queries, map<string,string> & mate1, map<string,string> & mate2, int const & numMismatches, ofstream & detailsFile, ofstream & mate1File, ofstream & mate2File)
{
	// Compute alignments
	alignment align1 = computeAlignment( queries, mate1["seq"], numMismatches);
	alignment align2 = computeAlignment( queries, mate2["seq"], numMismatches);
	// Classify mates
	string classification = align1.strand+"_"+align2.strand;
	// Chop fastq sequences and qualities
	string seq1  = mate1["seq" ].substr( 0, align1.start );
	string qual1 = mate1["qual"].substr( 0, align1.start );
	string seq2  = mate2["seq" ].substr( 0, align2.start );
	string qual2 = mate2["qual"].substr( 0, align2.start );
	// Print details
	detailsFile << mate1["id"]+"\t"+to_string(align1.start)+"\t"+to_string(align1.end)+"\t"+align1.strand+"\t"+align1.type+"\t"+to_string(align2.start)+"\t"+to_string(align2.end)+"\t"+align2.strand+"\t"+align2.type+"\n";
	// Print non-ambiguous fastq mates, attaching the classification to their ids.
	if ( (align1.start>=15) & (align2.start>=15) ){
		string read_id  = mate1["id"]  + "__" + classification;
		string read_id2 = mate1["id2"] + "__" + classification;
		mate1File << read_id+"\n"+seq1+"\n"+read_id2+"\n"+qual1+"\n";
		mate2File << read_id+"\n"+seq2+"\n"+read_id2+"\n"+qual2+"\n";
	}
	return 0;
}

string getComplement( string const & forward )
{
	string complement;
    for( string::size_type i = 0; i < forward.size(); ++i ) {
        if      ( (forward[i] == 'A') | (forward[i] == 'a') ) complement += 'T';
        else if ( (forward[i] == 'C') | (forward[i] == 'c') ) complement += 'G';
        else if ( (forward[i] == 'G') | (forward[i] == 'g') ) complement += 'C';
        else if ( (forward[i] == 'T') | (forward[i] == 't') ) complement += 'A';
    }	
	return complement;
}
string getReverse( string const & forward )
{
    string reversed="";
    for( unsigned i=0; i<forward.length(); ++i) {
        reversed = forward[i]+reversed ;
    }
    return reversed;
}

/////////////////////////////////////////////////////////////////////////////////////
// Execution
int main()
{
	string fastqFileName1="/home/yu68/Stitch-seq/new_data_Aug2013/ACCTRm_dupPE_stitch_seq_R1.fastq";
	string fastqFileName2="/home/yu68/Stitch-seq/new_data_Aug2013/ACCTRm_dupPE_stitch_seq_R2.fastq";
	string detailsFileName="matesDetails.txt";
	string mate1FileName="mate1.fastq";
	string mate2FileName="mate2.fastq";
	int numMismatches = 2;
	
	ifstream fastqFile1 (fastqFileName1.c_str()); 
	ifstream fastqFile2 (fastqFileName2.c_str());
	ofstream detailsFile ( detailsFileName.c_str() );
	ofstream mate1File ( mate1FileName.c_str() ); 
	ofstream mate2File ( mate2FileName.c_str() );

	// Define all possible query orientations
	string forward = "CTAGTAGCCCATGCAATGCGAGGA";
	string reverse = getReverse(forward);
	string forward_complement = getComplement(forward);
	string reverse_complement = getComplement(reverse);
	vector<string> queries;
	queries.push_back( forward );
	queries.push_back( forward_complement );
	queries.push_back( reverse );
	queries.push_back( reverse_complement );

	// Loop over fastq read mates
	int count = 0;
	string line1, line2;
	map<string,string> mate1, mate2;
	vector<string> keys;
	keys.push_back("id");keys.push_back("seq");keys.push_back("id2");keys.push_back("qual");
	while ( ! fastqFile1.eof() ) { 
		++count;
		if( count%1000000==0) cout << count << endl;
		// Read four lines at a time
		for( unsigned i=0; i<4; ++i ){
			getline(fastqFile1,line1);	
			getline(fastqFile2,line2);
			mate1[keys[i]] = line1;
			mate2[keys[i]] = line2;	
		}
		printMates( queries, mate1, mate2, numMismatches,detailsFile,mate1File,mate2File);
	}

	// Close files	
	fastqFile1.close();
	fastqFile2.close();
	detailsFile.close();
	mate1File.close();
	mate2File.close();

    return 0;
}
