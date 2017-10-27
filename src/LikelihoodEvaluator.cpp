/*
Copyright or © or Copr. Centre National de la Recherche Scientifique
contributor : Bastien Boussau (2009-2013)

bastien.boussau@univ-lyon1.fr

This software is a bioinformatics computer program whose purpose is to
simultaneously build gene and species trees when gene families have
undergone duplications and losses. It can analyze thousands of gene
families in dozens of genomes simultaneously, and was presented in
an article in Genome Research. Trees and parameters are estimated
in the maximum likelihood framework, by maximizing theprobability
of alignments given the species tree, the gene trees and the parameters
of duplication and loss.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/



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

/*extern "C" {
#include <pll/pllInternal.h>
}*/

#include <iostream>
#include <cstdio>
#include <string>
#include <fstream>
#include <boost/graph/graph_traits.hpp>
#include <algorithm>

#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/Node.h>
#include <Bpp/Phyl/TreeTemplateTools.h>

#include "Constants.h"
#include "LikelihoodEvaluator.h"
#include "ReconciliationTools.h"


// constants taken from RAXML
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


using namespace std;
using namespace bpp;



int LikelihoodEvaluator::hackmode = 0;

void print_PLL_param(pInfo *partition)
{
  print<double>("frequencies: ", partition->frequencies, 4);
  print<double>("gamma rates: ", partition->gammaRates, 4);
  print<double>("subst : ", partition->substRates, 6);
}

class LibpllNodeProperty: public Clonable {
  public:
    LibpllNodeProperty(unsigned int id) : id(id) {}

    virtual Clonable *  clone () const {
      return new LibpllNodeProperty(id);
    }
  unsigned int id;
};

void LikelihoodEvaluator::PLL_initializePLLInstance(){
  WHEREAMI( __FILE__ , __LINE__ );
  /* Set the PLL instance attributes */
  PLL_attributes.rateHetModel     = PLL_GAMMA;
  PLL_attributes.fastScaling      = PLL_TRUE;
  PLL_attributes.saveMemory       = PLL_FALSE;
  PLL_attributes.useRecom         = PLL_FALSE;
  PLL_attributes.randomNumberSeed = 0xDEADBEEF;
  PLL_attributes.numberOfThreads  = 8;            /* This only affects the pthreads version */

  PLL_instance = pllCreateInstance (&PLL_attributes);
}


LikelihoodEvaluator::LikelihoodEvaluator(map<string, string> params):
  params(params), alternativeTree(00), initialized(false), aligmentFilesForPllWritten_(false), logLikelihood(0), pll_model_already_initialized_(false), currentTreeinfo(0)
{
  WHEREAMI( __FILE__ , __LINE__ );
  loadDataFromParams();
  tolerance_ = 0.5;
}

void LikelihoodEvaluator::loadDataFromParams(){
  WHEREAMI( __FILE__ , __LINE__ );

  rollbackRootInfo.edge = 0;
  rollbackRootInfo.son = 0;
  movesNumber = 0;
  
  // set the name of this evaluator
  istringstream tempName(ApplicationTools::getStringParameter("input.sequence.file",params,"rnd"));
  while(std::getline(tempName,name,'/'))
    ;
  tempName.str(name);
  tempName.clear();
  std::getline(tempName,name,'.');
  if(name.size() == 0)
    name = "unnamed";

  cout << "@@Instanciating a Likelihood Evaluator named " << name << endl;

  string methodString = ApplicationTools::getStringParameter("likelihood.evaluator",params,"PLL");
  method = (methodString == "PLL"? PLL:BPP);

  // loading data, imported from GeneTreeLikelihood (Bastien)

  std::vector <std::string> spNames;
  bool cont = false;
  alphabet =  getAlphabetFromOptions(params, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Alphabet"));

  sites = getSequencesFromOptions(params, alphabet, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Sequences"));

  substitutionModel = getModelFromOptions(params, alphabet, sites, cont);

  if (!cont)
    throw(Exception("Unable to load this family: Model"));


  rateDistribution = getRateDistributionFromOptions(params, substitutionModel, cont);


  //For the moment, we have not filtered sequences,
  //so we don't want to spend time computing a tree that we may discard soon.
  //Instead we just create a random tree.
  std::vector<std::string> leafNames = sites->getSequencesNames();
  tree = TreeTemplateTools::getRandomTree(leafNames, false);

  //Get the scaler variable
  scaler_ = ApplicationTools::getDoubleParameter("sequence.likelihood.scaler", params, 1.0, "", false, false);


/*
  try
  {
    bool cont = true;
    tree = getTreeFromOptions(params, alphabet, sites, substitutionModel, rateDistribution, cont);
  }
  catch (std::exception& e)
  {
    std::cout << e.what() <<"; Unable to get a proper gene tree for family <<file<< avoiding this family."<<std::endl;
    cont=false;
  }*/
}


// TODO: is clonable?
// LikelihoodEvaluator::LikelihoodEvaluator(LikelihoodEvaluator &levaluator){
//   nniLk = levaluator.nniLk->clone();
//   tree = levaluator.tree->clone();
// }

void LikelihoodEvaluator::PLL_loadAlignment(string path)
{
  WHEREAMI( __FILE__ , __LINE__ );
  /* Parse a PHYLIP/FASTA file */
  PLL_alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, path.c_str());
  if (!PLL_alignmentData)
  {
    throw Exception("PLL: Error while parsing " + path);
  }
}


void LikelihoodEvaluator::PLL_loadNewick_fromFile(string path)
{
  WHEREAMI( __FILE__ , __LINE__ );
  PLL_newick = pllNewickParseFile(path.c_str());
  if (!PLL_newick)
  {
    throw Exception("PLL: Error while parsing newick file");
  }
  if (!pllValidateNewick (PLL_newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    throw Exception("PLL: Invalid phylogenetic tree.");
  }
}

void LikelihoodEvaluator::PLL_loadNewick_fromString(string newick)
{
  WHEREAMI( __FILE__ , __LINE__ );
  PLL_newick = pllNewickParseString (newick.c_str());
  if (!PLL_newick)
  {
    throw Exception("PLL: Error while parsing newick string: " + newick);
  }
  if (!pllValidateNewick (PLL_newick))  /* check whether the valid newick tree is also a tree that can be processed with our nodeptr structure */
  {
    throw Exception("PLL: Invalid phylogenetic tree.");
  }
}


void LikelihoodEvaluator::PLL_loadPartitions(string path)
{
  WHEREAMI( __FILE__ , __LINE__ );
  /* Parse the partitions file into a partition queue structure */
  PLL_partitionInfo = pllPartitionParse (path.c_str());
  

  /* Validate the partitions */
  if (!pllPartitionsValidate (PLL_partitionInfo, PLL_alignmentData))
  {
    throw Exception("PLL: Partitions do not cover all sites.");
  }

  /* Commit the partitions and build a partitions structure */
  PLL_partitions = pllPartitionsCommit (PLL_partitionInfo, PLL_alignmentData);
  PLL_partitions->partitionData[0]->frequencies = 0;
  /* We don't need the the intermedia partition queue structure anymore */
  pllQueuePartitionsDestroy (&PLL_partitionInfo);

  /* eliminate duplicate sites from the alignment and update weights vector */
  pllAlignmentRemoveDups (PLL_alignmentData, PLL_partitions);
}

void LikelihoodEvaluator::PLL_connectTreeAndAlignment()
{
  WHEREAMI( __FILE__ , __LINE__ );
  //std::cout << "** LikelihoodEvaluator::PLL_connectTreeAndAlignment " << std::endl;
  
  pllTreeInitTopologyNewick (PLL_instance, PLL_newick, PLL_FALSE);
  WHEREAMI( __FILE__ , __LINE__ );

  // cout << "PLL: Connect the alignment and partition structure with the tree structure" << std::endl ;
  /* Connect the alignment and partition structure with the tree structure */
  if (!pllLoadAlignment (PLL_instance, PLL_alignmentData, PLL_partitions))
  {
    throw Exception("PLL: Incompatible tree/alignment combination.");
  }

  pllNewickParseDestroy(&PLL_newick);

  WHEREAMI( __FILE__ , __LINE__ );
}


void LikelihoodEvaluator::initialize_BPP_nniLk()
{
  WHEREAMI( __FILE__ , __LINE__ );
  //std::cout << "** LikelihoodEvaluator::initialize_BPP_nniLk " << std::endl;
  nniLk = new NNIHomogeneousTreeLikelihood(*tree, *sites, substitutionModel, rateDistribution, mustUnrootTrees, verbose);

  nniLk->initParameters();
  nniLk->initialize();
  if(logLikelihood == 0)
  logLikelihood = - nniLk->getValue() * scaler_;
}


struct pll_sequence{
  pll_sequence(char *label, char *seq, unsigned int len):
    label(label),
    seq(seq),
    len(len) {}
  char *label;
  char *seq;
  unsigned int len;
};

typedef std::vector<pll_sequence> pll_sequences;

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
    sequences.push_back(pll_sequence(head, seq, seq_len));
    length = seq_len;
  }
  int count = sequences.size();;
  char** buffer = (char**)malloc(count * sizeof(char *));
  for (unsigned int i = 0; i < count; ++i) {
    buffer[i] = sequences[i].seq;
  }
  unsigned int *weights = pll_compress_site_patterns(buffer, pll_map_nt, count, &length);
  for (unsigned int i = 0; i < count; ++i) {
    sequences[i].len = length;
  }
  free(buffer);
  pll_fasta_close(fasta);
  return weights;
}

pll_unode_t *LikelihoodEvaluator::get_pll_utree_root(pll_utree_t * utree)
{
  return utree->nodes[utree->tip_count + utree->inner_count - 1];
}

pll_utree_t * LikelihoodEvaluator::create_utree(bpp::TreeTemplate< bpp::Node > *bpptree)
{
  std::string newick = bpp::TreeTemplateTools::treeToParenthesis(*bpptree);
  pll_rtree_t * rtree = pll_rtree_parse_newick_string(newick.c_str());
  pll_utree_t * utree = pll_rtree_unroot(rtree);
  pll_rtree_destroy(rtree, free);
  pll_utree_reset_template_indices(get_pll_utree_root(utree), utree->tip_count);
  return utree;
}


void fill_leaves_rec(pll_unode_t *node, vector<string> &leaves)
{
  if (!node->next) { // leaf
    leaves.push_back(string(node->label));
  } else {
    fill_leaves_rec(node->next->back, leaves);
    fill_leaves_rec(node->next->next->back, leaves);
  }
}


int cbtrav(pll_unode_t *node) {
  return 1;
}

void LikelihoodEvaluator::mapUtreeToBPPTree(pll_utree_t *utree, bpp::TreeTemplate< bpp::Node > *bpptree, bool bppStrict)
{
  if (hackmode <= 1) {
    return;
  }
  //std::cout << "** LikelihoodEvaluator::mapUtreeToBPPTree" << std::endl;
  unsigned int temp = 0;
  pll_utree_traverse(currentTreeinfo->root, PLL_TREE_TRAVERSE_POSTORDER, cbtrav, utree->nodes, &temp);
  //std::cout << "  bpp tree to map " << printer.getBPPNodeString((bpptree)->getRootNode(), false, true) << std::endl;
  //std::cout << "  utree    to map " << printer.getTreeinfoString(currentTreeinfo, false, true) << std::endl;
  std::vector< bpp::Node * > nodes = bpptree->getNodes();
  std::map<vector<string>, int> bppLeavesToId;
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    vector<string> bppLeaves = TreeTemplateTools::getLeavesNames(*nodes[i]);
    for (unsigned int j = 0; j < bppLeaves.size(); ++j) {
      if (!bppStrict) {
        bppLeaves[j] = realToStrict[bppLeaves[j]];
      }
      if (!bppLeaves.size()) {
          std::cout << "invalid conv in mapfct" << std::endl;
      }
    }
    std::sort(bppLeaves.begin(), bppLeaves.end());
    bppLeavesToId[bppLeaves] = nodes[i]->getId();
    //std::cout << "bpp " << nodes[i]->getId() << " leaves : ";
    //print_v<string>(bppLeaves);
  }

  unsigned int count = 0;
  for (unsigned int i = 0; i < utree->tip_count + utree->inner_count; ++i) {
    //pll_unode_t *node = (i == utree->tip_count + utree->inner_count) ? currentTreeinfo->root : utree->nodes[i];
    pll_unode_t *node = utree->nodes[i];
    vector<string> leaves;
    fill_leaves_rec(node, leaves);
    std::sort(leaves.begin(), leaves.end());
    if (!bppLeavesToId[leaves] && node->next) {
      leaves.clear();
      fill_leaves_rec(node->next, leaves);
      std::sort(leaves.begin(), leaves.end());
    }
    if (!bppLeavesToId[leaves] && node->next) {
      leaves.clear();
      fill_leaves_rec(node->next->next, leaves);
      std::sort(leaves.begin(), leaves.end());
    }
    //std::cout << "libpll " << i << " leaves : ";
    //print_v<string>(leaves);
    LibpllNodeProperty prop(i);
    if (bppLeavesToId[leaves])
      count++;
    bpptree->getNode(bppLeavesToId[leaves])->setNodeProperty("libpll", prop);
    //std::cout << i << "libll " << node->node_index << " to bpp " << bppLeavesToId[leaves] << std::endl;
  }
  //std::cout << "sucess count " << count << std::endl;
}

void LikelihoodEvaluator::reset_libpll_tree()
{
  if (hackmode <= 1) {
    return;
  }
  destroy_treeinfo();
}

void LikelihoodEvaluator::destroy_treeinfo()
{
  //std::cout << "** LikelihoodEvaluator::destroy_treeinfo" << std::endl;
  if (currentTreeinfo)
  {
    for (unsigned int i = 0; i < currentTreeinfo->partition_count; ++i)
    {
      if (currentTreeinfo->partitions[i])
        pll_partition_destroy(currentTreeinfo->partitions[i]);
    }

    pll_utree_graph_destroy(currentTreeinfo->root, NULL);
    pllmod_treeinfo_destroy(currentTreeinfo);
  }
  currentTreeinfo = 0;

}

pllmod_treeinfo_t * LikelihoodEvaluator::build_treeinfo(bool alternative)
{
  //std::cout << "** LikelihoodEvaluator::build_treeinfo" << std::endl;


  WHEREAMI( __FILE__ , __LINE__ );
  // partitions descriptors
  unsigned int partitions_number = 1; // todobenoit handle partitions!!!

  // sequences 
  std::string fileName = fileNamePrefix + "alignment.fasta"; 
  const char* fasta_file = fileName.c_str();
  pll_sequences sequences;
  unsigned int *pattern_weights = read_from_fasta(fasta_file, sequences);

  // tree
  bpp::TreeTemplate<Node> *treeToBuild = alternative ? alternativeTree : tree;
  pll_utree_t * utree = create_utree(treeToBuild);  
  currentUtree = utree;
  std::map<std::string, int> labelling;
  for (unsigned int i = 0 ; i < utree->tip_count; ++i)
  {
    pll_unode_t *node = utree->nodes[i];
    assert(!node->next);
    labelling[node->label] = i;
  }
 
  unsigned int brlen_linkage = PLLMOD_TREE_BRLEN_SCALED; // todobenoit is it the same model as raxml?
  pll_unode_t *uroot = get_pll_utree_root(utree); //todobenoit why does choice of the root matter?
  pllmod_treeinfo_t * treeinfo = pllmod_treeinfo_create(uroot, 
      utree->tip_count, partitions_number, brlen_linkage);

  // pll_attribute
  unsigned int attribute = PLL_ATTRIB_ARCH_AVX | PLL_ATTRIB_SITE_REPEATS;

  // pll_partitions
  pll_partition_t *partition = pll_partition_create(utree->tip_count,
      utree->inner_count,
      4,                // states.
      sequences[0].len, // sites
      1,                // rate_matrices
      utree->edge_count,// prob_matrices
      4,                // categories
      utree->edge_count,// scalers
      attribute);       // attr
  
  pll_set_pattern_weights(partition, pattern_weights);
  // add sequences to partitions
  const unsigned int *charmap = pll_map_nt; // todobenoit do not hardcode that
  for (unsigned int i = 0; i < sequences.size(); ++i) {
    unsigned int tip = labelling[strictToReal[sequences[i].label]];
    pll_set_tip_states(partition, tip, charmap, sequences[i].seq);
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
    std::cout << pll_part->frequencies << " " << pll_part->gammaRates <<" " << pll_part->substRates << std::endl;
    pll_set_category_rates(partition, pll_part->gammaRates);
    pll_set_frequencies(partition, 0, pll_part->frequencies);
    pll_set_subst_params(partition, 0, pll_part->substRates);
  }

  

  // treeinfo and partition
  int params_to_optimize = PLLMOD_OPT_PARAM_ALL; // todobenoit see what we should optimize
  unsigned int params_indices[4] = {0,0,0,0}; // todobenoit do not hardcode
  pllmod_treeinfo_init_partition(treeinfo, 0, partition,
        params_to_optimize,
        PLL_GAMMA_RATES_MEAN, // todobenoit: to check
        1.0, // todobenoit what is this alpha
        params_indices,
        0); // todobenoit check that we don't need it
  needFullOptim = true;
  return treeinfo;
}

void LikelihoodEvaluator::optimize_treeinfo(pllmod_treeinfo_t *treeinfo)
{
  double previousLogl = get_likelihood_treeinfo(treeinfo); 
  double newLogl = previousLogl;
  unsigned int iterations = 0;
  std::cout << "before ll = " << get_likelihood_treeinfo(currentTreeinfo) << std::endl;
  do {
    previousLogl = newLogl;
    optimize_treeinfo_iter(treeinfo);
    newLogl = get_likelihood_treeinfo(treeinfo);
    iterations++;
  } while (false ); //newLogl - previousLogl > tolerance_);
  std::cout << "Optimization iterations " << iterations << std::endl;
  logLikelihood = get_likelihood_treeinfo(currentTreeinfo);
  std::cout << "after: ll = " << logLikelihood<< std::endl;
}


double LikelihoodEvaluator::optimize_treeinfo_iter(pllmod_treeinfo_t *treeinfo)
{
  double new_loglh;
  unsigned int params_to_optimize = treeinfo->params_to_optimize[0]; // todobenoit read it from treeinfo
  
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

  /* optimize ALPHA */
  if (params_to_optimize & PLLMOD_OPT_PARAM_ALPHA)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
                                                      PLLMOD_OPT_PARAM_ALPHA,
                                                      PLLMOD_OPT_MIN_ALPHA,
                                                      PLLMOD_OPT_MAX_ALPHA,
                                                      RAXML_PARAM_EPSILON);
  }

  if (params_to_optimize & PLLMOD_OPT_PARAM_PINV)
  {
    new_loglh = -1 * pllmod_algo_opt_onedim_treeinfo(treeinfo,
                                                      PLLMOD_OPT_PARAM_PINV,
                                                      PLLMOD_OPT_MIN_PINV,
                                                      PLLMOD_OPT_MAX_PINV,
                                                      RAXML_PARAM_EPSILON);
  }

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
                                                                    tolerance_,
                                                                    brlen_smooth_factor * RAXML_BRLEN_SMOOTHINGS,
                                                                    -1,  /* radius */
                                                                    1,    /* keep_update */
                                                                    treeinfo->parallel_context,
                                                                    treeinfo->parallel_reduce_cb
                                                                    );
  }

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

}


double LikelihoodEvaluator::get_likelihood_treeinfo(pllmod_treeinfo_t *treeinfo)
{
  return pllmod_treeinfo_compute_loglh(treeinfo, 0);
}

void LikelihoodEvaluator::initialize_PLL()
{
  WHEREAMI( __FILE__ , __LINE__ );

  //std::cout << "** LikelihoodEvaluator::initialize_PLL" << std::endl;
  loadStrictNamesFromAlignment_forPLL();
  writeAlignmentFilesForPLL();
  alpha_ = 1.0;
  PLL_initializePLLInstance();
  PLL_loadAlignment(fileNamePrefix + "alignment.fasta");
  PLL_loadPartitions(fileNamePrefix + "partition.txt");
  if(logLikelihood == 0) {
    logLikelihood = PLL_evaluate(&tree) * scaler_;
  }  
}


void LikelihoodEvaluator::setTree(TreeTemplate<Node> * newTree)
{
  WHEREAMI( __FILE__ , __LINE__ );
  //std::cout << "LikelihoodEvaluator::setTree" << std::endl;
  //std::cout << newTree << " " << printer.getBPPNodeString(newTree->getRootNode(), true, true) << std::endl;
  if(!isInitialized()){
    if(tree)
      delete tree;
    tree = newTree->clone();
  }else
    throw Exception("This evaluator has already be initialized. It is forbidden to modify it now.");
}
  
pll_unode_t *LikelihoodEvaluator::getLibpllNode(bpp::Node *node)
{
  unsigned int libpllid =  dynamic_cast<LibpllNodeProperty *>
    (node->getNodeProperty("libpll"))->id;
  return currentUtree->nodes[libpllid];
}

bool areEquals(pll_unode_t *n1, pll_unode_t *n2) 
{
  if (n2->next)
    return n1== n2 || n1 == n2->next || n1 == n2->next->next;
  else
    return n1 == n2;
}

// returns the branch between n1 and n2 with the same index as n2
pll_unode_t *get_branch(pll_unode_t *n1, pll_unode_t *n2) {
 

  if (!n2->next) {
    return areEquals(n2->back, n1) ? n2 : 0;
  }
  
  if (areEquals(n2->back, n1))   {
    return n2;
  }  else if (areEquals(n2->next->back, n1)) {
    return n2->next;
  }  else if (areEquals(n2->next->next->back, n1)) {
    return n2->next->next;
  }  else {
    return 0;
  }
}

void LikelihoodEvaluator::applyNNIRoot(bpp::Node *bppParent,
  bpp::Node *bppGrandParent,
  bpp::Node *bppSon, bpp::Node *bppUncle)
{
  std::cout << "** LikelihoodEvaluator::applyNNIRoot" << std::endl;
  pll_unode_t *parent, *son, *uncle;
  try {
    parent = getLibpllNode(bppParent);
    son = getLibpllNode(bppSon);
    uncle = getLibpllNode(bppUncle);
  }
  catch (Exception e) {
    std::cout << "Exception ! " << e.what()<< std::endl;
    return;
  }
  pll_unode_t *edge = get_branch(uncle, parent);
  if (!edge) {
    std::cout << "LikelihoodEvaluator::applyNNIRoot: Impossible to find the good NNI move" << std::endl;
    return;
  }
  if (areEquals(edge->next->back, son)) {
    son = edge->next->back;
  } else if (areEquals(edge->next->next->back, son)) {
    son = edge->next->next->back;
  } else {
    std::cout << "Impossible to find son" << std::endl;
  }
  double t1 = edge->length;
  double t2 = son->length;
  rollbackRootInfo.t1 = t1;
  rollbackRootInfo.t2 = t2;
  rollbackRootInfo.edge = edge;
  rollbackRootInfo.son = son;
  pllmod_utree_set_length(son, t1/2.0 + t2);
  pllmod_utree_set_length(edge, t1/2.0);
}

  
void LikelihoodEvaluator::rebuildTreeinfoFromTree()
{
  std::cout << "**LikelihoodEvaluator::rebuildTreeinfoFromTree" << std::endl;
  destroy_treeinfo();
  currentTreeinfo = build_treeinfo(false);
}
  
void LikelihoodEvaluator::applyNNI(bpp::Node *bppParent, 
    bpp::Node *bppGrandParent,
    bpp::Node *bppSon, bpp::Node *bppUncle,
    bpp::Node *bppRoot)
{
  if (hackmode <= 1) {
    return;
  }
  //std::cout << "** LikelihoodEvaluator::applyNNI" << std::endl;
  movesNumber++;
  //std::cout << "moves: " << movesNumber << std::endl;
  if (bppParent->getId() == bppRoot->getId() ||
      bppGrandParent->getId() == bppRoot->getId()) {
    if (bppGrandParent->getId() == bppRoot->getId()) {
      applyNNIRoot(bppParent, bppGrandParent, bppSon,bppUncle);
      return;
    }
  }
  pll_unode_t *parent, *grandParent, *son, *uncle;
  try {
    parent = getLibpllNode(bppParent);
    grandParent = getLibpllNode(bppGrandParent);
    son = getLibpllNode(bppSon);
    uncle = getLibpllNode(bppUncle);
  }
  
  catch (Exception e) {
    std::cout << "Exception ! " << e.what()<< std::endl;
    return;
  }
  pll_unode_t *edge = get_branch(grandParent, parent);
  
  if (!edge) { 
    std::cout << "LikelihoodEvaluator::applyNNI: Impossible to find the good NNI move" << std::endl;
    return;
  }
 
  bool sonNext = areEquals(edge->next->back, son);
  bool uncleNext = areEquals(edge->back->next->back, uncle);

  unsigned int move = (sonNext == uncleNext) ?
    PLL_UTREE_MOVE_NNI_LEFT : PLL_UTREE_MOVE_NNI_RIGHT;
    
  
   
  if (!pllmod_utree_nni(edge, move, &rollbackInfo)) {
    std::cout << "failed applying nni : " << pll_errmsg << std::endl;
  }
  
}

void LikelihoodEvaluator::rollbackLastMove()
{
  if (hackmode <= 1) {
    return;
  }
  //std::cout << "** LikelihoodEvaluator::rollbackLastMove " << movesNumber << std::endl;
  movesNumber--;
  if (rollbackRootInfo.edge) {
    pllmod_utree_set_length(rollbackRootInfo.edge, rollbackRootInfo.t1);
    pllmod_utree_set_length(rollbackRootInfo.son, rollbackRootInfo.t2);
    rollbackRootInfo.edge = 0;
    return;
  }
  if (rollbackInfo.NNI.edge == 0) {
    std::cout << "no move to rollback" << std::endl;
    return;
  }
  if (PLL_SUCCESS != pllmod_tree_rollback(&rollbackInfo)) {
    std::cout << "rollback FAILED" << std::endl;
  }
  rollbackInfo.NNI.edge = 0;

}

void LikelihoodEvaluator::utreeRealToStrict(pllmod_treeinfo_t *treeinfo)
{
  pll_utree_t * utree = pll_utree_wraptree(treeinfo->root, treeinfo->tip_count);
  for (unsigned int i = 0; i < utree->tip_count; ++i) {
    pll_unode_t *tip = utree->nodes[i];
    std::string newlabel = realToStrict[tip->label];
    free(tip->label);
    unsigned int len = newlabel.size();
    tip->label = (char *)calloc(len + 1, sizeof(char));
    memcpy(tip->label, newlabel.c_str(), len);
  }
  // todobenoit free utree
}



void reroot(bpp::TreeTemplate<bpp::Node>* tree, 
    const set<string> &leaves1,
    const set<string> &leaves2)
{
    // the plyogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (tree)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());

      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        tree->newOutGroup(*currSon);

    }

    // if not, we will try all the internal branches as potential roots
    if(!tree->isRooted())
    {
      vector<Node*> outgroupCandidates = tree->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());
      
        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        tree->newOutGroup(*currCandidate);
      }
    }


    if(!tree->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    }


}

void saveRoot(bpp::TreeTemplate<bpp::Node>* tree, set<string> &leaves1, set<string> &leaves2)
{
  Node* root = tree->getRootNode();
  Node* son1 = root->getSon(0);
  Node* son2 = root->getSon(1);
  vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
  vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
  leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
  leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
  tree->unroot();
}

void updateTreeToEvaluate(bpp::TreeTemplate<Node> **treeToEvaluate, pllmod_treeinfo_t *treeinfo, BenoitPrinter &printer) {
  std::string newStr = printer.getTreeinfoString(treeinfo, true, false);
  stringstream outputNewickString;
  outputNewickString.str(newStr);
  Newick outputNewick;
  delete *treeToEvaluate;
  *treeToEvaluate = outputNewick.read(outputNewickString);

}


double LikelihoodEvaluator::libpll_evaluate_iterative(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)

{
  //std::cout << "LikelihoodEvaluator::libpll_evaluate_iterative" << std::endl;
  /*
  std::cout << "first with real PLL ";
  bpp::TreeTemplate<bpp::Node>* temp = (*treeToEvaluate)->clone();
  std::cout << realPLL_evaluate(&temp) << std::endl;
  delete temp;
  */

  if (!currentTreeinfo) {
    double res = libpll_evaluate_fromscratch(treeToEvaluate);
    return res;
  }
  
  bool wasRooted = ((*treeToEvaluate)->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted){
    saveRoot((*treeToEvaluate), leaves1, leaves2);
  }

  double result_ll = get_likelihood_treeinfo(currentTreeinfo);
  //std::cout << "ll it =" << result_ll << std::endl;
  if (needFullOptim) {
    optimize_treeinfo(currentTreeinfo); 
    needFullOptim = false;
  }
  result_ll = get_likelihood_treeinfo(currentTreeinfo);
  //std::cout << "libpll ll = " << result_ll << std::endl;
 
   updateTreeToEvaluate(treeToEvaluate, currentTreeinfo, printer);

  if (wasRooted) {
    reroot(*treeToEvaluate, leaves1, leaves2);
  }
  mapUtreeToBPPTree(currentUtree, *treeToEvaluate, true);
  //std::cout << "  output tree : "
  //  << printer.getBPPNodeString((*treeToEvaluate)->getRootNode(), false, true)
  //  << std::endl;
  return result_ll;
}

double LikelihoodEvaluator::libpll_evaluate_fromscratch(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)
{
  WHEREAMI( __FILE__ , __LINE__ );
  std::cout << "LikelihoodEvaluator::libpll_evaluate_fromscratch" << std::endl;
  /*
  std::cout << "first with real PLL " << std::endl;
  bpp::TreeTemplate<bpp::Node>* temp = (*treeToEvaluate)->clone();
  std::cout << realPLL_evaluate(&temp) << std::endl;
  delete temp;
  */
  if (currentTreeinfo) {
    destroy_treeinfo();
  }
  TreeTemplate<Node>* inputBPPTree = (*treeToEvaluate)->clone();
  Newick inputNewick;
  stringstream inputNewickString;
  inputNewick.write(**treeToEvaluate,inputNewickString);
  // save root to find it back
  bool wasRooted = (inputBPPTree->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted){
    saveRoot(inputBPPTree, leaves1, leaves2);
  }
 
  currentTreeinfo = build_treeinfo(alternativeTree != 0);
  double result_ll = get_likelihood_treeinfo(currentTreeinfo);
  std::cout << "libpll ll = " << result_ll << std::endl;
  optimize_treeinfo(currentTreeinfo);
  needFullOptim = false;
  result_ll = get_likelihood_treeinfo(currentTreeinfo);
  std::cout << "libpll ll = " << result_ll << std::endl;

  std::string newStr = printer.getTreeinfoString(currentTreeinfo, true, false, false);
  stringstream outputNewickString;
  outputNewickString.str(newStr);
  Newick outputNewick;
  delete *treeToEvaluate;
  *treeToEvaluate = outputNewick.read(outputNewickString);
  //std::cout << "PLL tree: " << printer.getTreeinfoString(currentTreeinfo, true, false, false, REAL_TO_STRICT) << std::endl;
  if (wasRooted) {
    reroot(*treeToEvaluate, leaves1, leaves2);
  }
  //std::cout << "output bpp tree " << printer.getBPPNodeString((*treeToEvaluate)->getRootNode(), false, false) << std::endl;
  mapUtreeToBPPTree(currentUtree, *treeToEvaluate, true);
  delete inputBPPTree;
 

  return result_ll;
}

double LikelihoodEvaluator::PLL_evaluate(TreeTemplate<Node>** treeToEvaluate)
{
  WHEREAMI( __FILE__ , __LINE__ );
  static unsigned int count = 0;
  //std::cout << "** LikelihoodEvaluator::PLL_evaluate" << count++ << std::endl;
  if (hackmode == 1) {
    return libpll_evaluate_fromscratch(treeToEvaluate);
  } else if (hackmode == 2) {
    return libpll_evaluate_iterative(treeToEvaluate);
  } else {
    return realPLL_evaluate(treeToEvaluate);
  }
}
 
double LikelihoodEvaluator::realPLL_evaluate(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)
{

  std::cout << "Real PLL evaluate " << std::endl;
  //TODO debug remove
  Newick debugTree;
  stringstream debugSS;
  debugTree.write(**treeToEvaluate,debugSS);
//  cout << "tree to evaluate:\n" << debugSS.str() << endl;

  // preparing the tree
  TreeTemplate<Node>* treeForPLL = (*treeToEvaluate)->clone();
  
  
  Newick newickForLibpll;
  stringstream newickStingForLibpll;
  newickForLibpll.write(**treeToEvaluate,newickStingForLibpll); 
  
  convertTreeToStrict(treeForPLL);

  // getting the root
  bool wasRooted = (treeForPLL->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted)
  {
    Node* root = treeForPLL->getRootNode();
    Node* son1 = root->getSon(0);
    Node* son2 = root->getSon(1);
    vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
    vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
    leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
    leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
    treeForPLL->unroot();
  }


  Newick newickForPll;
  stringstream newickStingForPll;
  newickForPll.write(*treeForPLL,newickStingForPll);
  PLL_loadNewick_fromString(newickStingForPll.str());
  // delete treeForPLL; //todobenoit (temporary moved to the end)

  // processing by PLL
  PLL_connectTreeAndAlignment();
  // pllSetFixedAlpha(alpha_, 0, PLL_partitions, PLL_instance);
  // pllSetFixedBaseFrequencies(baseFreq_, 4, 0, PLL_partitions, PLL_instance);
  // pllSetFixedSubstitutionMatrix(subsMatrix_, 6, 0, PLL_partitions, PLL_instance);
  WHEREAMI( __FILE__ , __LINE__ );
  if(!pll_model_already_initialized_){
    std::cout << "init model" << std::endl;
    pllInitModel(PLL_instance, PLL_partitions);
    pll_model_already_initialized_ = true;
  }
  else
    pllEvaluateLikelihood (PLL_instance, PLL_partitions, PLL_instance->start, PLL_TRUE, PLL_FALSE);
  WHEREAMI( __FILE__ , __LINE__ );

  //std::cout << "Old partition before opt:" << std::endl;
  //printOldPll(PLL_partitions);

 // pllOptimizeBranchLengths (PLL_instance, PLL_partitions, 64);
 // pllOptimizeModelParameters(PLL_instance, PLL_partitions, 0.1);
 
/**
 *  hackmode == 0:old PLL
 *  hackmode == 1: new libpll2
 * */

  double result_ll = 0.0;
  char *newickStr = 0;
    std::cout << "PLL ll before opt = " << PLL_instance->likelihood << std::endl;
    pllOptimizeModelParameters(PLL_instance, PLL_partitions, tolerance_);
    std::cout << "PLL ll after  opt = " << PLL_instance->likelihood << std::endl;
    result_ll = PLL_instance->likelihood;
    pllTreeToNewick(PLL_instance->tree_string, PLL_instance, PLL_partitions, PLL_instance->start->back, true, true, 0, 0, 0, true, 0,0);
    newickStingForPll.str(PLL_instance->tree_string);
    //std::cout << "PLL tree: " << newickStingForPll.str() << std::endl;
  
  // getting the new tree with new branch lengths

  
  //debug
  //cout << "DEBUG returned tree from PLL, LikelihoodEvaluator l367 \n" << newickStingForPll.str() << endl;

  delete *treeToEvaluate;
  *treeToEvaluate = newickForPll.read(newickStingForPll);
  

  // getting the likelihood and then deleting PLL_instance
  double PLL_instance_likelihood = result_ll * scaler_;

  //cout << "DEBUG PLL loglk LikelihoodEvaluator l375 \n" << PLL_instance_likelihood << endl;


  //re-rooting if needed
  if(wasRooted)
  {
    // the plyogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (*treeToEvaluate)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());

      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currSon);

    }

    // if not, we will try all the internal branches as potential roots
    if(!(*treeToEvaluate)->isRooted())
    {
      vector<Node*> outgroupCandidates = (*treeToEvaluate)->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());



        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currCandidate);
      }
    }


    if(!(*treeToEvaluate)->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    }

  }

  delete treeForPLL; //todobenoit (temporary moved to the end)

  //std::cout << "DEBUG Before restoring: "<< TreeTemplateTools::treeToParenthesis(**treeToEvaluate) <<std::endl;
  restoreTreeFromStrict(*treeToEvaluate);
  //std::cout << "DEBUG AFTER restoring: "<< TreeTemplateTools::treeToParenthesis(**treeToEvaluate) <<std::endl;

  //TODO debug remove
 /* stringstream debugSS2;
  debugTree.write(**treeToEvaluate,debugSS2);
  cout << "Final tree for BPP" << debugSS2.str() << endl;
  */
  free(newickStr);
  return(PLL_instance_likelihood);
}

double LikelihoodEvaluator::BPP_evaluate(TreeTemplate<Node>** treeToEvaluate)
{
  WHEREAMI( __FILE__ , __LINE__ );

  // preparing the tree
  TreeTemplate<Node>* treeForBPP = (*treeToEvaluate)->clone();

  // getting the root
  bool wasRooted = (treeForBPP->isRooted() ? true : false);
  set<string> leaves1, leaves2;
  if(wasRooted)
  {
    Node* root = treeForBPP->getRootNode();
    Node* son1 = root->getSon(0);
    Node* son2 = root->getSon(1);
    vector<string> leaves1Vector = TreeTemplateTools::getLeavesNames(*son1);
    vector<string> leaves2Vector = TreeTemplateTools::getLeavesNames(*son2);
    leaves1.insert(leaves1Vector.begin(),leaves1Vector.end());
    leaves2.insert(leaves2Vector.begin(),leaves2Vector.end());
    treeForBPP->unroot();
  }


  NNIHomogeneousTreeLikelihood * drlk = 0;
  drlk  = new NNIHomogeneousTreeLikelihood (**treeToEvaluate,
                                              *(nniLk->getData()),
                                              nniLk->getSubstitutionModel(),
                                              nniLk->getRateDistribution(),
                                              true, false);

  drlk->initialize();
  auto_ptr<BackupListener> backupListener;
  int tlEvalMax = 100;
  OutputStream* messageHandler = 0 ;
  OptimizationTools::optimizeBranchLengthsParameters(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*> (drlk), drlk->getParameters(), backupListener.get(), tolerance_, tlEvalMax, messageHandler, messageHandler, 0);

  delete *treeToEvaluate;
  *treeToEvaluate = static_cast< TreeTemplate<Node>* > (drlk->getTree().clone() );

 //re-rooting if needed
  if(wasRooted)
  {
    // the phylogenetic root is between the current topological
    // root and one of the three sons
    // so, one of the three sons of a root is the outgroup
    vector<Node*> rootSons = (*treeToEvaluate)->getRootNode()->getSons();
    for(vector<Node*>::iterator currSon = rootSons.begin(); currSon != rootSons.end(); currSon++)
    {
      vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currSon);
      set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());

      if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currSon);

    }

    // if not, we will try all the internal branches as potential roots
    if(!(*treeToEvaluate)->isRooted())
    {
      vector<Node*> outgroupCandidates = (*treeToEvaluate)->getNodes();
      for(vector<Node*>::iterator currCandidate = outgroupCandidates.begin(); currCandidate != outgroupCandidates.end(); currCandidate++)
      {
        vector<string> currLeavesVector = TreeTemplateTools::getLeavesNames(**currCandidate);
        set<string> currLeavesSet(currLeavesVector.begin(),currLeavesVector.end());

        if((currLeavesSet.size() == leaves1.size() && currLeavesSet == leaves1) || (currLeavesSet.size() == leaves2.size() && currLeavesSet == leaves2))
        (*treeToEvaluate)->newOutGroup(*currCandidate);
      }
    }


    if(!(*treeToEvaluate)->isRooted())
    {
      cout << "Unable to re-root the tree, I will give up, sorry." << endl;
      cout << "l1.s = " << leaves1.size() << ". l2.s = " << leaves2.size() << endl;
      throw Exception("Unable to re-root the tree.");
    }

  }

  return -drlk->getValue() * scaler_;


}





// CLONABLE?
// LikelihoodEvaluator* LikelihoodEvaluator::clone()
// {
//   return(new LikelihoodEvaluator(this));
// }


void LikelihoodEvaluator::unload()
{
  if(!initialized)
    return;
  initialized = false;
  if(method == PLL){
    //initialized = true;
    //return;
    if(pll_model_already_initialized_){
      pllAlignmentDataDestroy(PLL_alignmentData);
      pllPartitionsDestroy(PLL_instance, &PLL_partitions);
      pllDestroyInstance(PLL_instance);
      pll_model_already_initialized_ = false;
    }
  }
  else
  {
    delete nniLk;
    if(nniLkAlternative)
      delete nniLkAlternative;
  }
}

LikelihoodEvaluator::~LikelihoodEvaluator()
{
  unload();
  delete tree;
  if(alternativeTree)
    delete alternativeTree;
  if(method=PLL){
    remove(string(fileNamePrefix + "alignment.fasta").c_str());
    remove(string(fileNamePrefix + "partition.txt").c_str());
  }
  WHEREAMI( __FILE__ , __LINE__ );

}

void LikelihoodEvaluator::initialize()
{
  WHEREAMI( __FILE__ , __LINE__ );
  //std::cout << "LikelihoodEvaluator::initialize" << std::endl;
  if(initialized)
    return;
  //checking the alignment and the tree contain the same number of sequences
  if(sites->getNumberOfSequences() != tree->getNumberOfLeaves()){
    ostringstream errorMessage;
    errorMessage << "\nNumber of sequences (here: "<< sites->getNumberOfSequences() << "sequences) must match to number of leaves in the tree (here: "<< tree->getNumberOfLeaves() << "leaves). I give up.";
    throw Exception(errorMessage.str());
  }

  // ### common requirements for initialization
  if(method == PLL)
    initialize_PLL();
  else
    initialize_BPP_nniLk();

  //
  initialized = true;
}


void LikelihoodEvaluator::setAlternativeTree(TreeTemplate< Node >* newAlternative)
{
  WHEREAMI( __FILE__ , __LINE__ );
  //std::cout << "** LikelihoodEvaluator::setAlternativeTree" << std::endl;
  if(!initialized)
    initialize();
  if(alternativeTree != 00)
    delete alternativeTree;
  alternativeTree = newAlternative->clone();

  if(method == PLL){
    alternativeLogLikelihood = PLL_evaluate( &alternativeTree ) * scaler_;
  }
  else
  {
  /*  if ( nniLkAlternative)
      delete nniLkAlternative;*/
    alternativeLogLikelihood = BPP_evaluate( &alternativeTree ) * scaler_;
   /* nniLkAlternative =  new NNIHomogeneousTreeLikelihood (*alternativeTree, *sites, substitutionModel, rateDistribution, mustUnrootTrees, verbose);
    alternativeLogLikelihood = nniLkAlternative->getLogLikelihood();*/
  }
}

void LikelihoodEvaluator::acceptAlternativeTree()
{
  WHEREAMI( __FILE__ , __LINE__ );
  delete tree;
  WHEREAMI( __FILE__ , __LINE__ );
  tree = alternativeTree;
  alternativeTree = 00;
  logLikelihood = alternativeLogLikelihood;
}





TreeTemplate< Node >* LikelihoodEvaluator::getAlternativeTree()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return alternativeTree;
}

bool LikelihoodEvaluator::isInitialized()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return initialized;
}

Alphabet* LikelihoodEvaluator::getAlphabet()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return alphabet;
}

double LikelihoodEvaluator::getAlternativeLogLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return alternativeLogLikelihood;
}

double LikelihoodEvaluator::getLogLikelihood()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return logLikelihood;
}

DiscreteDistribution* LikelihoodEvaluator::getRateDistribution()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return rateDistribution;
}

VectorSiteContainer* LikelihoodEvaluator::getSites()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return sites;
}

TreeTemplate< Node >* LikelihoodEvaluator::getTree()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return tree;
}

SubstitutionModel* LikelihoodEvaluator::getSubstitutionModel()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return substitutionModel;
}


void LikelihoodEvaluator::loadStrictNamesFromAlignment_forPLL()
{
  WHEREAMI( __FILE__ , __LINE__ );
  vector<string> seqNames = sites->getSequencesNames();
  string currName;
  ostringstream currStrictName;
  for(unsigned int currIdx = 0; currIdx != seqNames.size(); currIdx++)
  {
    currName = seqNames.at(currIdx);
    if(currName.at(currName.size()-1) == '\r')
      currName = currName.substr(0,(currName.size()-1));

//     cout << "Curr sequence name = ..." << currName << "..." << endl;


//     // determining current postfix with letters instead of integers.
//     // Thanks PLL to be so strict!...
//     // (it is not a clean base conversion)
//     ostringstream currPrefix;
//     unsigned int currIntegerPrefix = currIdx+1;
//     while(currIntegerPrefix != 0)
//     {
//       unsigned int unit = currIntegerPrefix % 10;
//       currPrefix << (char)('A' + unit);
//       currIntegerPrefix = currIntegerPrefix / 10;
//     }
//     cout << "curr prefix = " << currPrefix.str() << endl;

    currStrictName.clear();
    currStrictName.str("");
    currStrictName << "seq" << currIdx;
    realToStrict[currName] = currStrictName.str();
    strictToReal[currStrictName.str()] = currName;

  }
  printer = BenoitPrinter(realToStrict, strictToReal);
}

void LikelihoodEvaluator::convertTreeToStrict(TreeTemplate< Node >* targetTree)
{
  WHEREAMI( __FILE__ , __LINE__ );

  vector<Node*> leaves = targetTree->getLeaves();
  for(vector<Node*>::iterator currLeaf = leaves.begin(); currLeaf != leaves.end(); currLeaf++)
  {
    string currLeafName = (*currLeaf)->getName();
    if(currLeafName.at(currLeafName.size()-1) == '\r')
    {
      currLeafName = currLeafName.substr(0,(currLeafName.size()-1));
      (*currLeaf)->setName(currLeafName);
    }
    map<string,string>::iterator found = realToStrict.find(currLeafName);
    if(found == realToStrict.end())
    {
      cout << "Unable to find sequence named ++" << currLeafName << "++ in the alignment." << endl;
    }
    else
    {
      (*currLeaf)->setName(found->second);
    }
  }
}

void LikelihoodEvaluator::restoreTreeFromStrict(TreeTemplate< Node >* targetTree)
{
  WHEREAMI( __FILE__ , __LINE__ );
  vector<Node*> leaves = targetTree->getLeaves();
  for(vector<Node*>::iterator currLeaf = leaves.begin(); currLeaf != leaves.end(); currLeaf++){
    (*currLeaf)->setName(strictToReal[(*currLeaf)->getName()]);
  }
}


void LikelihoodEvaluator::writeAlignmentFilesForPLL()
{
  if(aligmentFilesForPllWritten_)
    return;
  WHEREAMI( __FILE__ , __LINE__ );
  fileNamePrefix = "tmpPLL_" + name + "_" ;
  ofstream alignementFile(string(fileNamePrefix + "alignment.fasta").c_str(), ofstream::out);

  //preparing the file for the alignment
  BasicSequence currSequence(sites->getAlphabet());
  //DEBUG
  cout << "Writing an alignment for PLL with " << sites->getNumberOfSequences() << " sequences. File: " << string(fileNamePrefix + "alignment.fasta") << endl;

  for(unsigned int currSeqIndex = 0; currSeqIndex != sites->getNumberOfSequences(); currSeqIndex++)
  {
    currSequence = sites->getSequence(currSeqIndex);
    string currSequenceName = currSequence.getName();
    if(currSequenceName.at(currSequenceName.size()-1) == '\r')
      currSequenceName = currSequenceName.substr(0,(currSequenceName.size()-1));
    string currSequenceStr = currSequence.toString();
    if(alphabet->getSize() != 4)
        replace(currSequenceStr.begin(), currSequenceStr.end(), '*', 'X');
    alignementFile << ">" << realToStrict[currSequenceName] << "\n" << currSequenceStr << "\n";
  }
  alignementFile.close();
  ofstream partitionFile(string(fileNamePrefix + "partition.txt").c_str(), ofstream::out);
if (alphabet->getSize() == 4) {
    if (substitutionModel->getName()!="GTR") {
     std::cout << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes GTR for nucleotide models." <<std::endl;
     cout.flush();
     std::cerr << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes GTR for nucleotide models." <<std::endl;
     cerr.flush();
     //MPI::COMM_WORLD.Abort(1); //SHOULD BE CORRECTED 13062017
     exit(-1);
    }
    if (ApplicationTools::getBooleanParameter("codon.partition", params, false, "") == true ) {
      partitionFile << "DNA, p1=1-" << sites->getNumberOfSites() << "/3\n";
      partitionFile << "DNA, p2=2-" << sites->getNumberOfSites() << "/3\n";
      partitionFile << "DNA, p3=3-" << sites->getNumberOfSites() << "/3\n";
    }
     else
      partitionFile << "DNA, p1=1-" << sites->getNumberOfSites() << "\n";
}
else if (alphabet->getSize() == 20) {
    if (substitutionModel->getName().substr(0,4)=="LG08") {
        partitionFile << "LG, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else if (substitutionModel->getName().substr(0,5)=="WAG01") {
        partitionFile << "WAG, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else if (substitutionModel->getName().substr(0,5)=="JTT92") {
        partitionFile << "JTT, p1=1-" << sites->getNumberOfSites() << "\n";
    }
    else {
     std::cout << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes LG08, WAG01, JTT92 for protein models." <<std::endl;
     cout.flush();
     std::cerr << "Error: model unrecognized for optimization with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only recognizes LG08, WAG01, JTT92 for protein models." <<std::endl;
     cerr.flush();
     //MPI::COMM_WORLD.Abort(1); //SHOULD BE CORRECTED 13062017  
     exit(-1);
    }
}
else {
  std::cout << "Error: alphabet incompatible with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only works with DNA/RNA or Protein alphabets." <<std::endl;
  cout.flush();
  std::cerr << "Error: alphabet incompatible with PLL. Maybe you want to use BPP by setting the option: likelihood.evaluator=BPP. PLL only works with DNA/RNA or Protein alphabets." <<std::endl;
  cerr.flush();
  //MPI::COMM_WORLD.Abort(1); //SHOULD BE CORRECTED 13062017  
  exit(-1);
}
  partitionFile.close();
  aligmentFilesForPllWritten_ = true;
}

LikelihoodEvaluator::LikelihoodEvaluator(LikelihoodEvaluator const &leval):
params(leval.params), initialized(false), PLL_instance(00), PLL_alignmentData(00), PLL_newick(00), PLL_partitions(00), PLL_partitionInfo(00), tree(00), alternativeTree(00), nniLk(00), nniLkAlternative(00), substitutionModel(00), rateDistribution(00), sites(00), alphabet(00), aligmentFilesForPllWritten_(false), logLikelihood(0), pll_model_already_initialized_(false)
{
  WHEREAMI( __FILE__ , __LINE__ );

  loadDataFromParams();
  tree = leval.tree;
  sites = leval.sites;

  if(leval.initialized)
    initialize();
}

LikelihoodEvaluator* LikelihoodEvaluator::clone()
{
  WHEREAMI( __FILE__ , __LINE__ );
  return new LikelihoodEvaluator(*this);
}

LikelihoodEvaluator::LikelihoodEvaluator(const Tree* tree, const SiteContainer* alignment, SubstitutionModel* model, DiscreteDistribution* rateDistribution, std::map<std::string, std::string> par, bool mustUnrootTrees, bool verbose):
initialized(false), PLL_instance(00), PLL_alignmentData(00), PLL_newick(00), PLL_partitions(00), PLL_partitionInfo(00), tree(00), alternativeTree(00), nniLk(00), nniLkAlternative(00), substitutionModel(00), rateDistribution(00), sites(00), alphabet(00), params(par), aligmentFilesForPllWritten_(false), logLikelihood(0), pll_model_already_initialized_(false)
{
  WHEREAMI( __FILE__ , __LINE__ );
  rollbackRootInfo.edge = 0;
  rollbackRootInfo.son = 0;
  movesNumber = 0;
  this->tree = dynamic_cast<TreeTemplate<Node> *>(tree->clone());
  this->substitutionModel = model->clone();
  this->rateDistribution = rateDistribution->clone();
  this->sites = dynamic_cast<VectorSiteContainer*>(alignment->clone());

  string methodString = ApplicationTools::getStringParameter("likelihood.evaluator",params,"PLL");
  this->method = (methodString == "PLL"? PLL:BPP);

  initialize();

}
