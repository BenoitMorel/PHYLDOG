
extern "C" {
#include <pll/pll.h>
#include <pllmodules/pll.h>
#include <pllmodules/pllmod_algorithm.h>
#include <pllmodules/pll_binary.h>
#include <pllmodules/pll_msa.h>
#include <pllmodules/pll_optimize.h>
#include <pllmodules/pll_tree.h>
#include <pllmodules/pllmod_util.h>
}

#include <iostream>
#include <vector>
#include <map>
#define TOLERANCE 0.5

#define DEF_LH_EPSILON            0.1
#define OPT_LH_EPSILON            0.1
#define RAXML_PARAM_EPSILON       0.001  //0.01
#define RAXML_BFGS_FACTOR         1e7

#define RAXML_BRLEN_SMOOTHINGS    32
#define RAXML_BRLEN_DEFAULT       0.1
#define RAXML_BRLEN_MIN           1.0e-6
#define RAXML_BRLEN_MAX           100.
#define RAXML_BRLEN_TOLERANCE     1.0e-7

#define RAXML_FREERATE_MIN        0.001
#define RAXML_FREERATE_MAX        100.

#define RAXML_BRLEN_SCALER_MIN    0.01
#define RAXML_BRLEN_SCALER_MAX    100.

struct PLLData {
  pllInstanceAttr attributes;
  pllInstance *instance;
  pllAlignmentData *alignmentData;
  pllNewickTree *newick;
  partitionList *partitions;
  pllQueue *partitionInfo;
};


void optimize_libpll2(pllmod_treeinfo_t *treeinfo)
{
  // This code comes from RaxML
  double new_loglh;
  unsigned int params_to_optimize = treeinfo->params_to_optimize[0]; 

  /* optimize SUBSTITUTION RATES */
  if (params_to_optimize & PLLMOD_OPT_PARAM_SUBST_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_subst_rates_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_SUBST_RATE,
        PLLMOD_OPT_MAX_SUBST_RATE,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);
  }
  std::cout << "    subst rates " << new_loglh << std::endl;

  /* optimize BASE FREQS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREQUENCIES)
  {
    new_loglh = -1 * pllmod_algo_opt_frequencies_treeinfo(treeinfo,
        0,
        PLLMOD_OPT_MIN_FREQ,
        PLLMOD_OPT_MAX_FREQ,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);
  }
  std::cout << "    param frequ " << new_loglh << std::endl;

  std::cout << " alpha before: " << treeinfo->alphas[0] << std::endl;

  /* optimize ALPHA */
  if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_ALPHA,
        PLLMOD_OPT_MIN_ALPHA,
        PLLMOD_OPT_MAX_ALPHA,
        RAXML_PARAM_EPSILON);
  }

  std::cout << " alpha after: " << treeinfo->alphas[0] << std::endl;
  
  std::cout << "        alpha ll " << new_loglh << std::endl;
  
  if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_PINV,
        PLLMOD_OPT_MIN_PINV,
        PLLMOD_OPT_MAX_PINV,
        RAXML_PARAM_EPSILON);
  }
  std::cout << "           pinv " << new_loglh << std::endl;

  /* optimize FREE RATES and WEIGHTS */
  if (params_to_optimize & PLLMOD_OPT_PARAM_FREE_RATES)
  {
    new_loglh = -1 * pllmod_algo_opt_rates_weights_treeinfo (treeinfo,
        RAXML_FREERATE_MIN,
        RAXML_FREERATE_MAX,
        RAXML_BFGS_FACTOR,
        RAXML_PARAM_EPSILON);

    /* normalize scalers and scale the branches accordingly */
    if (treeinfo->brlen_linkage == PLLMOD_TREE_BRLEN_SCALED &&
        treeinfo->partition_count > 1)
      pllmod_treeinfo_normalize_brlen_scalers(treeinfo);

  }
  
  std::cout << "    free rates  " << new_loglh << std::endl;

  if (params_to_optimize & PLLMOD_OPT_PARAM_BRANCHES_ITERATIVE)
  {
    double brlen_smooth_factor = 0.25; // magical number from raxml
    new_loglh = -1 * pllmod_opt_optimize_branch_lengths_local_multi(treeinfo->partitions,
        treeinfo->partition_count,
        treeinfo->root,
        treeinfo->param_indices,
        treeinfo->deriv_precomp,
        treeinfo->brlen_scalers,
        RAXML_BRLEN_MIN,
        RAXML_BRLEN_MAX,
        TOLERANCE,
        brlen_smooth_factor * RAXML_BRLEN_SMOOTHINGS,
        -1,  /* radius */
        1,    /* keep_update */
        treeinfo->parallel_context,
        treeinfo->parallel_reduce_cb
        );
  }
  
  std::cout << "    bl iterativ " << new_loglh << std::endl;

  /* optimize brlen scalers, if needed */
  if (treeinfo->brlen_linkage == PLLMOD_TREE_BRLEN_SCALED &&
      treeinfo->partition_count > 1)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
        PLLMOD_OPT_PARAM_BRANCH_LEN_SCALER,
        RAXML_BRLEN_SCALER_MIN,
        RAXML_BRLEN_SCALER_MAX,
        RAXML_PARAM_EPSILON);

    /* normalize scalers and scale the branches accordingly */
    pllmod_treeinfo_normalize_brlen_scalers(treeinfo);
  }
  std::cout << "    bl linkage  " << new_loglh << std::endl;
}


void execute_pll(const char *newick, const char *part, const char *msa)
{
  PLLData plldata;
  plldata.attributes.rateHetModel     = PLL_GAMMA;
  plldata.attributes.fastScaling      = PLL_TRUE;
  plldata.attributes.saveMemory       = PLL_FALSE;
  plldata.attributes.useRecom         = PLL_FALSE;
  plldata.attributes.randomNumberSeed = 0xDEADBEEF;
  plldata.attributes.numberOfThreads  = 1;            /* This only affects the pthreads version */
  plldata.instance = pllCreateInstance (&plldata.attributes);
  
  plldata.alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, msa);
  if (!plldata.alignmentData) {
    std::cerr << "Error in pllParseAlignmentFile " << msa << std::endl;
    return;
  }
  
  
  plldata.partitionInfo = pllPartitionParse (part);

  if (!pllPartitionsValidate (plldata.partitionInfo, plldata.alignmentData)) {
    std::cerr << "Error in pllPartitionParse " << part << std::endl;
    return;
  }

  plldata.partitions = pllPartitionsCommit (plldata.partitionInfo, plldata.alignmentData);
  plldata.partitions->partitionData[0]->frequencies = 0;
  pllQueuePartitionsDestroy (&plldata.partitionInfo);
  pllAlignmentRemoveDups (plldata.alignmentData, plldata.partitions);

  plldata.newick = pllNewickParseFile (newick);
  if (!plldata.newick || !pllValidateNewick (plldata.newick)) {
    std::cerr << "Error in pllNewickParseFile " << newick << std::endl;
    return;
  }

  
  pllTreeInitTopologyNewick (plldata.instance, plldata.newick, PLL_FALSE);
  if (!pllLoadAlignment (plldata.instance, plldata.alignmentData, plldata.partitions)) {
    std::cerr << "Error in pllLoadAlignment" << std::endl;
    return;
  }
  pllNewickParseDestroy(&plldata.newick);

  pllInitModel(plldata.instance, plldata.partitions);
  
  std::cout << "PLL initial likelihood: " << plldata.instance->likelihood << std::endl;
  pllOptimizeModelParameters(plldata.instance, plldata.partitions, TOLERANCE);
  std::cout << "PLL final likelihood: " << plldata.instance->likelihood << std::endl;
}

struct pll_sequence{
  pll_sequence(char *label, char *seq, unsigned int len):
    label(label),
    seq(seq),
    len(len) {}
  char *label;
  char *seq;
  unsigned int len;
  ~pll_sequence() {
    free(label);
    free(seq);
  }
};

typedef std::vector<pll_sequence *> pll_sequences;

unsigned int * read_from_fasta(const char *fasta_file, pll_sequences &sequences)
{
  pll_fasta_t * fasta = pll_fasta_open(fasta_file, pll_map_fasta);
  char * head;
  long head_len;
  char *seq;
  long seq_len;
  long seqno;
  int length;
  while (pll_fasta_getnext(fasta, &head, &head_len, &seq, &seq_len, &seqno))
  {
    sequences.push_back(new pll_sequence(head, seq, seq_len));
    length = seq_len;
  }
  int count = sequences.size();;
  char** buffer = (char**)malloc(count * sizeof(char *));
  for (unsigned int i = 0; i < count; ++i) {
    buffer[i] = sequences[i]->seq;
  }
  unsigned int *weights = pll_compress_site_patterns(buffer, pll_map_nt, count, &length);
  for (unsigned int i = 0; i < count; ++i) {
    sequences[i]->len = length;
  }
  free(buffer);
  pll_fasta_close(fasta);
  return weights;
}

void execute_libpll2(const char *newick, const char *part, const char *msa)
{
  std::cout << "start execute libpll2 " << std::endl;
  unsigned int partitions_number = 1; 

  pll_sequences sequences;
  unsigned int *pattern_weights = read_from_fasta(msa, sequences);
  if (!pattern_weights) {
    std::cout << "Error in read_from_fasta "  << msa << std::endl;
    return;
  }
  
  pll_utree_t * utree = pll_utree_parse_newick(newick);
  if (!utree) {
    std::cout << "Error in pll_rtree_parse_newick_string " << newick << std::endl;
    return;
  }
  
  std::map<std::string, int> labelling;
  for (unsigned int i = 0 ; i < utree->tip_count; ++i)
  {
    pll_unode_t *node = utree->nodes[i];
    assert(!node->next);
    labelling[node->label] = i;
  }
  unsigned int brlen_linkage = PLLMOD_TREE_BRLEN_LINKED; 
  pll_unode_t *uroot = utree->nodes[utree->tip_count + utree->inner_count - 1] ; 
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(uroot, 
      utree->tip_count, partitions_number, brlen_linkage);
  if (!treeinfo) {
    std::cout << "Error cannot create treeinfo " << std::endl;
    return;
  }
  unsigned int attribute = PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_SITE_REPEATS;

  unsigned int numberOfStates = 4;
  const unsigned int *charmap = pll_map_nt;
  pll_partition_t *partition = pll_partition_create(utree->tip_count,
      utree->inner_count,
      numberOfStates,                // states.
      sequences[0]->len, // sites
      1,                // rate_matrices
      utree->edge_count,// prob_matrices
      4,                // categories
      utree->inner_count,// scalers
      attribute);       // attr
  std::cout << "partition create " << utree->tip_count << " " << utree->inner_count << " "
    << numberOfStates << " " << sequences[0]->len << " " << 1 << " " << utree->edge_count << 
    " " << 4 << " " << utree->inner_count << " " << attribute << std::endl;


  pll_set_pattern_weights(partition, pattern_weights);
  free(pattern_weights);
  for (unsigned int i = 0; i < sequences.size(); ++i) {
    unsigned int tip = labelling[sequences[i]->label];
    pll_set_tip_states(partition, tip, charmap, sequences[i]->seq);
    delete sequences[i];
  }

  std::cout << "Warning, hardcoded stuff here" << std::endl;
  std::cout << "THink about getting these parameters from the PLL instance" << std::endl;
  double alpha = 0.497118;
  double subst[6] = {1.143690, 3.515523, 0.721250, 0.864108, 4.319078, 1.000000};
  double gammaRates[4] = {0.136954, 0.476752, 1, 2.38629};
  double frequencies[4] = {0.168201, 0.321044, 0.340331, 0.170424};
  pll_compute_gamma_cats(alpha, 4, gammaRates, PLL_GAMMA_RATES_MEAN);
  pll_set_category_rates(partition, gammaRates);
  pll_set_frequencies(partition, 0, frequencies);
  pll_set_subst_params(partition, 0, subst);
  
  int params_to_optimize = PLLMOD_OPT_PARAM_ALL & ~PLLMOD_OPT_PARAM_FREE_RATES;
  unsigned int params_indices[4] = {0,0,0,0}; 
  pllmod_treeinfo_init_partition(treeinfo, 0, partition,
      params_to_optimize,
      PLL_GAMMA_RATES_MEAN,
      alpha, 
      params_indices,
      0); 
  // set missing branch lengths
  double length = 0.00001;
  for (unsigned int i = 0; i < utree->tip_count; ++i)
    if (!utree->nodes[i]->length)
      utree->nodes[i]->length = length;
  for (unsigned int i = utree->tip_count; i < utree->tip_count + utree->inner_count; ++i)
  {
    if (!utree->nodes[i]->length)
      utree->nodes[i]->length = length;
    if (!utree->nodes[i]->next->length)
      utree->nodes[i]->next->length = length;
    if (!utree->nodes[i]->next->next->length)
      utree->nodes[i]->next->next->length = length;
  } 

  double ll = pllmod_treeinfo_compute_loglh(treeinfo, 0);
  double newll = 0;
  std::cout << "LIBPLL-2 initial likelihood: " << ll << std::endl;
  unsigned int it = 0;
  do {
    ll = newll;
    optimize_libpll2(treeinfo);
    newll = pllmod_treeinfo_compute_loglh(treeinfo, 0);
    std::cout << "iteration " << it++ << " likelihood " << newll << std::endl;
  } while (fabs(newll - ll) > TOLERANCE);
  std::cout << "LIBPLL-2 final likelihood: " << newll << std::endl;
}

int main(int args, char ** argv)
{
  const char *newick = "data_investigation/plltree.newick";
  const char *msa = "data_investigation/pllmsa.fasta";
  const char *part = "data_investigation/partition.txt";
  execute_libpll2(newick, part,  msa);
  execute_pll(newick, part,  msa);
  std::cout << "end of the story" << std::endl;
  return 0;
}



