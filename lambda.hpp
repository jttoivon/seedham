#include <string>
#include <vector>

double
jussi_lambda(const std::vector<std::string>& background_seqs, 
	    const std::vector<std::string>& sequences);


double
hamming_lambda(const std::vector<std::string>& data1, 
		const std::vector<std::string>& data2,
		const std::string& consensus, int generation);
