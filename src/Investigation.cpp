
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

#define TOLERANCE 0.5

struct PLLData {
  pllInstanceAttr attributes;
  pllInstance *instance;
  pllAlignmentData *alignmentData;
  pllNewickTree *newick;
  partitionList *partitions;
  pllQueue *partitionInfo;
};

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

  // tree
  /*
  pll_utree_t * utree = createUtreeFromBPP(bppTree);  
  currentUtree = utree;
  std::map<std::string, int> labelling;
  for (unsigned int i = 0 ; i < utree->tip_count; ++i)
  {
    pll_unode_t *node = utree->nodes[i];
    assert(!node->next);
    labelling[node->label] = i;
  }

  unsigned int brlen_linkage = PLLMOD_TREE_BRLEN_SCALED; 
  pll_unode_t *uroot = getUtreeRoot(utree); 
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(uroot, 
      utree->tip_count, partitions_number, brlen_linkage);

  // pll_attribute
  unsigned int attribute = getBestLibpllAttribute();


  // pll_partitions
  unsigned int numberOfStates = sites->getAlphabet()->getSize();
  const unsigned int *charmap = 0;
  if (numberOfStates == 4) {
    charmap = pll_map_nt;
  } else if (numberOfStates == 20) {
    charmap = pll_map_aa;
  } else {
    std::cerr << "Error: Phyldog detected " << numberOfStates << " states " << std::endl;
    std::cout << "Currently supported states number: 4 or 20" << std::endl;
    return 0;
  }
  pll_partition_t *partition = pll_partition_create(utree->tip_count,
      utree->inner_count,
      numberOfStates,                // states.
      sequences[0]->len, // sites
      1,                // rate_matrices
      utree->edge_count,// prob_matrices
      4,                // categories
      utree->edge_count,// scalers
      attribute);       // attr

  pll_set_pattern_weights(partition, pattern_weights);
  free(pattern_weights);

  // add sequences to partitions
  for (unsigned int i = 0; i < sequences.size(); ++i) {
    unsigned int tip = labelling[strictToReal[sequences[i]->label]];
    pll_set_tip_states(partition, tip, charmap, sequences[i]->seq);
    delete sequences[i];
  }

  // model
  if (!PLL_partitions || !PLL_partitions->partitionData||  ! PLL_partitions->partitionData[0]
      || !PLL_partitions->partitionData[0]->substRates
      || !PLL_partitions->partitionData[0]->frequencies
      || !PLL_partitions->partitionData[0]->substRates
     ) {
    double subst[6] = {1, 1, 1, 1, 1, 1};
    double gammaRates[4] = {0.136954, 0.476752, 1, 2.38629};
    pll_set_category_rates(partition, gammaRates);
    pll_set_frequencies(partition, 0, &(substitutionModel->getFrequencies()[0]));
    pll_set_subst_params(partition, 0, subst);

  } else { 
    pInfo *pll_part = PLL_partitions->partitionData[0];
    pll_set_category_rates(partition, pll_part->gammaRates);
    pll_set_frequencies(partition, 0, pll_part->frequencies);
    pll_set_subst_params(partition, 0, pll_part->substRates);
  }

  // treeinfo and partition
  int params_to_optimize = PLLMOD_OPT_PARAM_ALL; 
  unsigned int params_indices[4] = {0,0,0,0}; 
  pllmod_treeinfo_init_partition(treeinfo, 0, partition,
      params_to_optimize,
      PLL_GAMMA_RATES_MEAN,
      1.0, 
      params_indices,
      0); 
  needFullOptim = true;
  return treeinfo;
*/
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



