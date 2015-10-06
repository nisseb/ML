
#include <iostream>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/seq_io.h>
#include <algorithm>

#define AMINO_SIZE 20

using namespace seqan;

void print_transition_matrix(std::vector< std::vector<double> > * matrix)
{
  std::cout << "Transition matrix:" << std::endl;
  for (std::vector< std::vector<double> >::const_iterator it = matrix->begin(); 
       it != matrix->end(); ++ it)
    {
      for (std::vector<double>::const_iterator itit = it->begin(); itit != it->end(); ++ itit)
	std::cout << *itit << '\t';
      
      std::cout << std::endl;
    }
}

void print_stationary_prob(std::vector<double> * vec)
{
  std::cout << "Stationary probabilities:" << std::endl;
  for (std::vector<double>::const_iterator it = vec->begin(); it != vec->end(); ++ it)
    {
      std::cout << *it << '\t';
    }
  std::cout << std::endl;
}

int find_character_index(const char * arr, size_t size, char c)
{
  const char* end = arr + size;
  const char* match = std::find(arr, end, c);
  return (end == match)? -1 : (match-arr);
}

void load_matrix_amino(std::istream* is,
		       std::vector< std::vector<double> >* matrix,
		       std::vector<double> *vec,
		       int *num,
		       const std::string& delim = " \t")
{
    using namespace std;

    string      line;
    string      strnum;

    // clear first
    matrix->clear();

    int j = 0;
    // parse line by line
    while (getline(*is, line) && j <= AMINO_SIZE+1)
    {
      if (line.at(0) == '#')
	{
	  //std::cout << line.at(0) << '\n';
	}
      else
	{
	  if (j < AMINO_SIZE)
	    {
	      matrix->push_back(vector<double>());
	    }

	  for (string::const_iterator i = line.begin(); i != line.end(); ++ i)
	    {
	      // If i is not a delim, then append it to strnum
	      if (delim.find(*i) == string::npos)
		{
		  strnum += *i;
		  if (i + 1 != line.end()) // If it's the last char, do not continue
                    continue;
		}

	      // if strnum is still empty, it means the previous char is also a
	      // delim (several delims appear together). Ignore this char.
	      if (strnum.empty())
                continue;

	      // If we reach here, we got a number. Convert it to double.
	      double       number;

	      istringstream(strnum) >> number;
	      if (j < AMINO_SIZE)       {matrix->back().push_back(number);}
	      else if (j == AMINO_SIZE)	{vec->push_back(number);}
	      else                      {*num = (int)number;}

	      strnum.clear();
	    }
	  j++;
	}
    }
}

/*
void get_cond_probs()
{
  //parameters
  int idx_curr = 0;
  int idx_prev = 0;
  char current_char;
  char previous_char;
  double prob;
  
  for (int i = 0; i < length(seqs); i++)
    {
      //std::cout << i << '\n';
      // first value
      probs.push_back(std::vector<double>());

      current_char = seqs[i][0];
      idx_curr = find_character_index(amino_acid_sequence, 
				      sizeof(char)*AMINO_SIZE, current_char);

      prob = stationary_prob.at(idx_curr);
      probs[i].push_back(prob);

      // all other values
  
      for (int j = 1; j < length(seqs[i]); j++)
	{
	  previous_char = current_char;
	  idx_prev = idx_curr;
	  current_char = seqs[i][j];
	  idx_curr = find_character_index(amino_acid_sequence, 
					  sizeof(char)*AMINO_SIZE, current_char);
	  //std::cout << idx_prev << ", " << idx_curr << '\n';
	  if (idx_curr < 0 || idx_prev < 0)
	    {
	      prob = 1;
	    }
	  else 
	    {
	      prob = transition_matrix[idx_prev][idx_curr];
	    }
	  probs[i].push_back(prob);
	}
  
    }
}
*/


int main(int argc, char *argv[] )
{

  if (argc != 5)
    {
      std::cout << "Proper usage: ./main fasta model1 model160 model161" << '\n';
      return 0;
    }
  else
    {
      //std::string dir = "../data/";
    }
  
  
  
  // read file, load transition matrix and close file
  std::ifstream is1(argv[2]);
  std::ifstream is160(argv[3]);
  std::ifstream is161(argv[4]);

  std::vector< std::vector<double> > transition_matrix;
  std::vector< std::vector<double> > transition_matrix160;
  std::vector< std::vector<double> > transition_matrix161;
  
  std::vector<double> stationary_prob;
  std::vector<double> stationary_prob160;
  std::vector<double> stationary_prob161;

  int num_seq;
  int num_seq160;
  int num_seq161;
  
  load_matrix_amino(&is1, &transition_matrix, &stationary_prob, &num_seq);
  load_matrix_amino(&is160, &transition_matrix160, &stationary_prob160, &num_seq160);
  load_matrix_amino(&is161, &transition_matrix161, &stationary_prob161, &num_seq161);

  is1.close();
  is160.close();
  is161.close();

  //print_transition_matrix(&transition_matrix1);
  //print_stationary_prob(&stationary_prob1);
  
  StringSet<CharString> ids;
  StringSet< String<AminoAcid> > seqs;

  SeqFileIn seqFileIn(argv[1]);
  
  try
    {
      readRecords(ids, seqs, seqFileIn);
    }
  catch (Exception const & e)
    {
      std::cout << "ERROR: " << e.what() << std::endl;
      return 1;
    }

  // Character order in transition matrix and statinary probabillity
  const char * amino_acid_sequence = {"IVLFCMAGTWSYPHEQDNKR"};

  
  
  /* Likelihood computation */
  // parameters

  
  int idx_curr = 0;
  int idx_prev = 0;
  char current_char;
  char previous_char;
  std::vector< std::vector<double> > probs;
  double prob;
  

  //std::cout << "Length ids: " << length(ids) << '\n';
  //std::cout << "Length ids2: " << length(ids[0]) << '\n';
  
  
  for (int i = 0; i < length(seqs); i++)
    {
      //std::cout << i << '\n';
      // first value
      probs.push_back(std::vector<double>());

      current_char = seqs[i][0];
      idx_curr = find_character_index(amino_acid_sequence, 
				      sizeof(char)*AMINO_SIZE, current_char);

      prob = stationary_prob.at(idx_curr);
      probs[i].push_back(prob);

      // all other values
  
      for (int j = 1; j < length(seqs[i]); j++)
	{
	  previous_char = current_char;
	  idx_prev = idx_curr;
	  current_char = seqs[i][j];
	  idx_curr = find_character_index(amino_acid_sequence, 
					  sizeof(char)*AMINO_SIZE, current_char);
	  //std::cout << idx_prev << ", " << idx_curr << '\n';
	  if (idx_curr < 0 || idx_prev < 0)
	    {
	      prob = 1;
	    }
	  else 
	    {
	      prob = transition_matrix[idx_prev][idx_curr];
	    }
	  probs[i].push_back(prob);
	}
  
    }
  

  // calculate the probabillity
  
  std::vector<double> cond_prob;
  double temp_sum;
  
  for (int i = 0; i < length(probs); i++)
    {
      temp_sum = 0;
      for (int j = 0; j < length(probs[i]); j++)
	{
	  temp_sum = temp_sum + log10(probs[i][j]);
	}
      cond_prob.push_back(temp_sum);
      std::cout << temp_sum << '\n';
    }

  for (int i = 0; i < 34; i++)
    {
      //std::cout << ids[i] << '\t' << seqs[i] << '\n';
    }
  

  std::cout << '\n';
  std::cout << "End of program" << '\n';

  return 1;
}


/*
std::cout << "dim1: " << length(transition_matrix[0]) << '\n';
std::cout << "dim2: " << length(transition_matrix) << '\n';
std::cout << "num_seq: " << num_seq << '\n';
std::cout << transition_matrix[19][19] << '\n';
*/
