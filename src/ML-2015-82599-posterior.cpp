
#include <iostream>
#include <stdio.h>
#include <seqan/sequence.h>  // CharString, ...
#include <seqan/stream.h>    // to stream a CharString into cout
#include <seqan/seq_io.h>
#include <algorithm>
#include <gmp.h>


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

void get_cond_probs(StringSet< String<AminoAcid> > seqs,
		    std::vector< std::vector<double> > transition_matrix,
		    std::vector<double> stationary_prob,
		    const char * amino_acid_sequence,
		    std::vector< std::vector<double> > & probs)
{
  //parameters
  int idx_curr = 0;
  int idx_prev = 0;
  char current_char;
  char previous_char;
  double prob;
  
  for (int i = 0; i < length(seqs); i++)
    {

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

void log_cond_probs(std::vector<double> & cond_prob, 
		    std::vector< std::vector<double> > probs)
{
  
  // log all conditional probabillities
  double temp_sum;
  
  for (int i = 0; i < length(probs); i++)
    {
      temp_sum = 0;
      for (int j = 0; j < length(probs[i]); j++)
	{
	  temp_sum = temp_sum + log10(probs[i][j]);
	}
      cond_prob.push_back(temp_sum);
      //std::cout << temp_sum << '\n';
    }
}

void calc_sequence_prob()
{
  
}

void print_result_to_file(double * prob1,
			  double * prob160,
			  double * prob161,
			  StringSet<CharString> ids,
			  char *filename,
			  int length)
{
  //std::ofstream out(filename);
  FILE *fp = fopen(filename, "w");
  
  fprintf(fp, "ID\tCOG1\tCOG160\tCOG161\n");
  
  
  for (int i = 0; i < length; i++)
    {
      fprintf(fp, "SEQ %d \t", i);
      fprintf(fp, "%.4f \t %.4f \t %.4f\n", prob1[i], prob160[i], prob161[i]);
    }

  //out.close();
  fclose(fp);
}

void print_result_to_file_gmp(mpf_t * prob1,
			      mpf_t * prob160,
			      mpf_t * prob161,
			      StringSet<CharString> ids,
			      char *filename,
			      int length)
{
  std::ofstream out(filename);
  
  out << "ID" << '\t';
  out << "COG1" << '\t' << "COG160" << '\t' << "COG161" << '\n';
  
  for (int i = 0; i < length; i++)
    {
      out << ids[i] << '\t';
      out << prob1[i] << '\t' << prob160[i] << '\t' << prob161[i] << '\n';
    }

  out.close();
}

void get_cond_prob_mpf(mpf_t * cond_prob_mpf, 
		    std::vector< std::vector<double> > probs)
{
  
  mpf_t a, b;
  mpf_init(a);
  mpf_init(b);
  
  for (int i = 0; i < probs.size(); i++)
    {
  
      mpf_set_d(b,1);
      for (int j = 0; j < probs[i].size(); j++)
	{
	  mpf_set_d(a, probs[i][j]);
	  mpf_mul(b,b,a);
	}

      mpf_init(cond_prob_mpf[i]);
      mpf_set(cond_prob_mpf[i],b);
    }

  mpf_clear(a);
  mpf_clear(b);
}

void final_calc(mpf_t *cond_prob_mpf, mpf_t *cond_prob_mpf160, mpf_t *cond_prob_mpf161,
		mpf_t *final_prob, mpf_t *specific_prob, mpf_t specific_pModel, int length,
		mpf_t pModel, mpf_t pModel160, mpf_t pModel161)
{
  mpf_t temp_u, temp_l, temp_lmul, temp_lmul160, temp_lmul161;
  mpf_init(temp_u);
  mpf_init(temp_l);
  mpf_init(temp_lmul);  
  mpf_init(temp_lmul160);
  mpf_init(temp_lmul161);
  
  for (int i = 0; i < length; i++)
    {
      /* upper */
      mpf_mul(temp_u, specific_prob[i], specific_pModel);
      /* lower */
      
      // mult
      mpf_mul(temp_lmul, cond_prob_mpf[i], pModel);
      mpf_mul(temp_lmul160, cond_prob_mpf160[i], pModel160);
      mpf_mul(temp_lmul161, cond_prob_mpf161[i], pModel161);
      
      // add
      mpf_add(temp_l, temp_lmul160, temp_lmul161);
      mpf_add(temp_l, temp_l, temp_lmul);
      
      // upper / lower
      mpf_init(final_prob[i]);
      mpf_div(final_prob[i], temp_u, temp_l);
    }
  
  mpf_clear(temp_u);
  mpf_clear(temp_l);
  mpf_clear(temp_lmul);
  mpf_clear(temp_lmul160);
  mpf_clear(temp_lmul161);
}

void get_doubles(double * res, mpf_t * vals, int length)
{
  for (int i = 0; i < length; i++)
    {
      res[i] = mpf_get_d(vals[i]);
    }
}

void log_probs(double * res, int length)
{
  for (int i = 0; i < length; i++)
    {
      res[i] = log10(res[i]);
    }
}


void clear_data(mpf_t * cond_prob_data, int length)
{
  for (int i = 0; i < length; i++)
    {
      mpf_clear(cond_prob_data[i]);
    }
}




int main(int argc, char *argv[] )
{

  if (argc != 6)
    {
      std::cout << "Proper usage: ./main fasta model1 model160 model161 output_file" << '\n';
      return 0;
    }

  /* load transition matrix and sequence count */
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

  /* Read amino sequences from file */
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

  /* Character order in transition matrix and statinary probabillity */
  const char * amino_acid_sequence = {"IVLFCMAGTWSYPHEQDNKR"};

  
  /* Likelihood computation using gmp */

  // parameters
  std::vector< std::vector<double> > probs;
  std::vector< std::vector<double> > probs160;
  std::vector< std::vector<double> > probs161;
  
  get_cond_probs(seqs, transition_matrix, stationary_prob,
		 amino_acid_sequence, probs);
  get_cond_probs(seqs, transition_matrix160, stationary_prob160,
		 amino_acid_sequence, probs160);
  get_cond_probs(seqs, transition_matrix161, stationary_prob161,
		 amino_acid_sequence, probs161);
  

  // calculate the log probabillity into cond_prob
  
  std::vector<double> cond_prob;
  std::vector<double> cond_prob160;
  std::vector<double> cond_prob161;

  log_cond_probs(cond_prob, probs);
  log_cond_probs(cond_prob160, probs160);
  log_cond_probs(cond_prob161, probs161);

  // P(model)
  double temp_sum = num_seq + num_seq160 + num_seq161;
  double modP = num_seq / temp_sum;
  double modP160 = num_seq160 / temp_sum;
  double modP161 = num_seq161 / temp_sum;

  // sum P(seq|model) * P(model) for every model
  
  calc_sequence_prob();

  /* Full calculation with gmp library */
  
  // Initiate variables and allocate

  mpf_t *cond_prob_mpf;
  mpf_t *cond_prob_mpf160;
  mpf_t *cond_prob_mpf161;
  
  cond_prob_mpf = new mpf_t[probs.size()];
  cond_prob_mpf160 = new mpf_t[probs160.size()];
  cond_prob_mpf161 = new mpf_t[probs161.size()];

  mpf_t *final_prob;
  mpf_t *final_prob160;
  mpf_t *final_prob161;
  
  final_prob = new mpf_t[probs.size()];
  final_prob160 = new mpf_t[probs160.size()];
  final_prob161 = new mpf_t[probs161.size()];
  
  mpf_t pModel, pModel160, pModel161;
  mpf_init(pModel);
  mpf_init(pModel160);
  mpf_init(pModel161);

  mpf_set_d(pModel, modP);
  mpf_set_d(pModel160, modP160);
  mpf_set_d(pModel161, modP161);

  // get p(sequence|model)
  get_cond_prob_mpf(cond_prob_mpf, probs);
  get_cond_prob_mpf(cond_prob_mpf160, probs160);
  get_cond_prob_mpf(cond_prob_mpf161, probs161);
    
  // calculate p(model|sequence)
  final_calc(cond_prob_mpf, cond_prob_mpf160, cond_prob_mpf161,
	     final_prob, cond_prob_mpf, pModel, probs.size(),
	     pModel, pModel160, pModel161);
  final_calc(cond_prob_mpf, cond_prob_mpf160, cond_prob_mpf161,
	     final_prob160, cond_prob_mpf160, pModel160, probs160.size(),
	     pModel, pModel160, pModel161);
  final_calc(cond_prob_mpf, cond_prob_mpf160, cond_prob_mpf161,
	     final_prob161, cond_prob_mpf161, pModel161, probs161.size(),
	     pModel, pModel160, pModel161);

  
  double *res;
  double *res160;
  double *res161;

  res = new double[probs.size()];
  res160 = new double[probs160.size()];
  res161 = new double[probs161.size()];
  
  get_doubles(res, final_prob, probs.size());
  get_doubles(res160, final_prob160, probs160.size());
  get_doubles(res161, final_prob161, probs161.size());

  /* Test
  for (int i = 0; i < probs.size(); i++)
    {
      std::cout << i << ": " << res[i] << std::endl;
    }
  */

  log_probs(res, probs.size());
  log_probs(res160, probs160.size());
  log_probs(res161, probs161.size());

  /* Test
  for (int i = 0; i < probs.size(); i++)
    {
      std::cout << i << ": " << res[i] << std::endl;
    }
  */
  
  /* Printing */
  // gmp, not log scale
  //print_result_to_file_gmp(final_prob, final_prob160, final_prob161, ids, argv[5], probs.size());
  
  // log scale
  print_result_to_file(res, res160, res161, ids, argv[5], probs.size());
  
  /* Clean up */
  clear_data(cond_prob_mpf, probs.size());
  clear_data(cond_prob_mpf160, probs160.size());
  clear_data(cond_prob_mpf161, probs161.size()); 

  clear_data(final_prob, probs.size());
  clear_data(final_prob160, probs160.size());
  clear_data(final_prob161, probs161.size());
  
  delete [] cond_prob_mpf;
  delete [] cond_prob_mpf160;
  delete [] cond_prob_mpf161;

  delete [] final_prob;
  delete [] final_prob160;
  delete [] final_prob161;

  delete [] res;
  delete [] res160;
  delete [] res161;

  std::cout << "End of program" << '\n';
  return 1;
}

