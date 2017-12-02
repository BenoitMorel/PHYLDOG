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


/**
 * I am duplicating the same function on purpose for profiling purpose
 * */
double getTreeinfoLikelihoodFull(pllmod_treeinfo_t *treeinfo, bool update_matrix)
{
  return pllmod_treeinfo_compute_loglh_flex(treeinfo, 0, update_matrix);
}

double getTreeinfoLikelihoodIncr(pllmod_treeinfo_t *treeinfo, bool update_matrix)
{
  return pllmod_treeinfo_compute_loglh_flex(treeinfo, 1, update_matrix);
}

double getTreeinfoLikelihood(pllmod_treeinfo_t *treeinfo, bool incremental = false, bool update_matrix = true)
{
  if (incremental) {
    return getTreeinfoLikelihoodIncr(treeinfo, update_matrix);
  } else {
    return getTreeinfoLikelihoodFull(treeinfo, update_matrix);
  }
}

double getTreeinfoLikelihoodNoRecompute(pllmod_treeinfo_t *treeinfo) 
{
  return treeinfo->partition_loglh[0];
}

/*
 * @brief Rollback for an identity move (no move)
 */
class IdentityRollback: public LikelihoodEvaluator::Rollback {
  public:
    IdentityRollback(double likelihood):
      LikelihoodEvaluator::Rollback(likelihood) {}

    virtual bool applyRollback() {
      return true;      
    }

  protected:
    virtual void localOptimizeNoCheck() {}
};

/*
 * @brief Rollback for a classical SPR move
 */
class SPRRollback: public LikelihoodEvaluator::Rollback {
  public:
    SPRRollback(pllmod_treeinfo_t *treeinfo,
        const pll_tree_rollback_t &rb,
        double likelihood):
      LikelihoodEvaluator::Rollback(likelihood),
      rb(rb),
      treeinfo(treeinfo) {}

    virtual bool applyRollback() {
      return PLL_SUCCESS == pllmod_tree_rollback(&rb);      
    }

  protected:
    /*
     * Optimize the three branches around the regrafted edge
     */
    virtual void localOptimizeNoCheck() {
      getTreeinfoLikelihood(treeinfo);
      pll_unode_t *edge = (pll_unode_t *)rb.SPR.prune_edge;
      pllmod_treeinfo_set_root(treeinfo, edge);
      getTreeinfoLikelihood(treeinfo, true);
      unsigned int params_indices[4] = {0,0,0,0}; 
      pllmod_opt_optimize_branch_lengths_local(
        treeinfo->partitions[0],
        treeinfo->root,
        params_indices,
        RAXML_BRLEN_MIN,
        RAXML_BRLEN_MAX,
        RAXML_BRLEN_TOLERANCE,
        RAXML_BRLEN_SMOOTHINGS,
        1,
        true);
    }

  private:
    pll_tree_rollback_t rb;
    pllmod_treeinfo_t *treeinfo;
};

/*
 * Special rollback used for SPR moves that don't change the
 * unrooted tree likelihood. In HYBRID mode, we set some branch 
 * lengths to 0.1 to be consistent with PLL/BPP. This requires
 * a rollback even if the tree didn't change
 */
class SPRIdentityHybridRollback: public LikelihoodEvaluator::Rollback {
  public:
    SPRIdentityHybridRollback(double previousLikelihood, pll_unode_t *nbf, pll_unode_t *ob) :
      Rollback(previousLikelihood),
      nbf(nbf), ob(ob), 
      nbf1(nbf->back->length),
      nbf2(nbf->back->next->length),
      nbf3(nbf->back->next->next->length),
      ob1(ob->back->length) {}

    virtual bool applyRollback() {
      pllmod_utree_set_length(nbf->back, nbf1);
      pllmod_utree_set_length(nbf->back->next, nbf2);
      pllmod_utree_set_length(nbf->back->next->next, nbf3);
      pllmod_utree_set_length(ob, ob1);
      return true;
    }
    
  protected:
    virtual void localOptimizeNoCheck() {
    }
  private:
    pll_unode_t *nbf; // newBrotherFather
    pll_unode_t *ob; // oldBrother
    double nbf1; // branch length (BL)
    double nbf2; // BL
    double nbf3; // BL
    double ob1;  // BL
};

/*
 * Rollback for a NNI move
 */
class NNIRollback: public LikelihoodEvaluator::Rollback {
  public:
    NNIRollback(pllmod_treeinfo_t *treeinfo, const pll_tree_rollback_t &rb, double likelihood): 
      LikelihoodEvaluator::Rollback(likelihood),
      rb(rb),
      treeinfo(treeinfo) {}

    virtual bool applyRollback() {return PLL_SUCCESS == pllmod_tree_rollback(&rb);}

  protected:
    /*
     * Optimizes the 5 branches around the NNI move
     */
    virtual void localOptimizeNoCheck() {
      unsigned int params_indices[4] = {0,0,0,0}; 
      pll_unode_t *edge = (pll_unode_t *)rb.NNI.edge;
      std::vector<pll_unode_t *> nodesToOptimize;
      nodesToOptimize.push_back(edge);
      if (edge->next) {
        nodesToOptimize.push_back(edge->next);
        nodesToOptimize.push_back(edge->next->next);
      }
      if (edge->back && edge->back->next) {
        nodesToOptimize.push_back(edge->back->next);
        nodesToOptimize.push_back(edge->back->next->next);
      }
      // might be incremental
      for (unsigned int i = 0; i < nodesToOptimize.size(); ++i) {
        pllmod_treeinfo_set_root(treeinfo, nodesToOptimize[i]);
        getTreeinfoLikelihood(treeinfo, i);
        
        pllmod_opt_optimize_branch_lengths_local(
          treeinfo->partitions[0],
          treeinfo->root,
          params_indices,
          RAXML_BRLEN_MIN,
          RAXML_BRLEN_MAX,
          RAXML_BRLEN_TOLERANCE,
          RAXML_BRLEN_SMOOTHINGS,
          0,
          true);
      }
    }
  private:
    pll_tree_rollback_t rb;
    pllmod_treeinfo_t *treeinfo;
};

/*
 * Rollback for a special case of NNI moves: when the root of the rooted tree
 * is implied in the NNI move, the unrooted tree topology is not changed but
 * stuff still happen on branch lengths...
 * This is similar to SPRIdentityHybridRollback
 */
class NNIRootRollback: public LikelihoodEvaluator::Rollback {
  public:
    NNIRootRollback(pll_unode_t *edge, pll_unode_t *son, double t1, double t2, double likelihood):
      LikelihoodEvaluator::Rollback(likelihood),
      edge(edge), son(son), t1(t1), t2(t2) {}

    virtual bool applyRollback() {
      pllmod_utree_set_length(edge, t1);
      pllmod_utree_set_length(son, t2);
      return true;
    }

  protected:
    virtual void localOptimizeNoCheck() {}
  private:
    pll_unode_t *edge;
    pll_unode_t *son;
    double t1;
    double t2;
};

/*
 * Node property used to map a libpll node to a bpp node
 */
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
  PLL_partitions = 0;

  destroyRollbacks();
  
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
  std::cout << " METHOD " << methodString << std::endl;
  if (methodString == "PLL") {
    this->method = PLL;
  } else if (methodString == "LIBPLL2") {
    this->method = LIBPLL2;
  } else if (methodString == "HYBRID") {
    this->method = HYBRID;
  } else {
    this->method = BPP;
  }

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
  WHEREAMI( __FILE__ , __LINE__ );

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

pll_unode_t *LikelihoodEvaluator::getUtreeRoot(pll_utree_t * utree)
{
  return utree->nodes[utree->tip_count + utree->inner_count - 1];
}

pll_utree_t * LikelihoodEvaluator::createUtreeFromBPP(const bpp::TreeTemplate< bpp::Node > *bpptree)
{
  std::string newick = bpp::TreeTemplateTools::treeToParenthesis(*bpptree);
  pll_rtree_t * rtree = pll_rtree_parse_newick_string(newick.c_str());
  pll_utree_t * utree = pll_rtree_unroot(rtree);
  pll_rtree_destroy(rtree, free);
  pll_utree_reset_template_indices(getUtreeRoot(utree), utree->tip_count);
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

void LikelihoodEvaluator::mapUtreeToBPPTree(pll_utree_t *utree, bpp::TreeTemplate< bpp::Node > *bpptree)
{
  //std::cout << "LikelihoodEvaluator::mapUtreeToBPPTree" << std::endl;
  if (method != LIBPLL2 && method != HYBRID) {
    return;
  }
  unsigned int temp = 0;
  pll_utree_traverse(currentTreeinfo->root, PLL_TREE_TRAVERSE_POSTORDER, cbtrav, utree->nodes, &temp);
  std::vector< bpp::Node * > nodes = bpptree->getNodes();
  std::map<vector<string>, int> bppLeavesToId;
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    vector<string> bppLeaves = TreeTemplateTools::getLeavesNames(*nodes[i]);
    std::sort(bppLeaves.begin(), bppLeaves.end());
    bppLeavesToId[bppLeaves] = nodes[i]->getId();
    //std::cout << "bpp " << nodes[i]->getId() << " leaves : ";
    //print_v<string>(bppLeaves);
  }
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
    bpptree->getNode(bppLeavesToId[leaves])->setNodeProperty("libpll", prop);
    //std::cout << i << "libll " << node->node_index << " to bpp " << bppLeavesToId[leaves] << std::endl;
  }
}

void LikelihoodEvaluator::destroyTreeinfo()
{
  if (currentTreeinfo)
  {
    for (unsigned int i = 0; i < currentTreeinfo->partition_count; ++i)
    {
      if (currentTreeinfo->partitions[i])
        pll_partition_destroy(currentTreeinfo->partitions[i]);
    }

    pll_utree_graph_destroy(currentTreeinfo->root, NULL);
    pllmod_treeinfo_destroy(currentTreeinfo);
    currentTreeinfo = 0;
  }
  if (currentUtree) {
    free(currentUtree->nodes);
    free(currentUtree);
    //pll_utree_destroy(currentUtree, free);
    currentUtree = 0;
  }
}

unsigned int getBestLibpllAttribute() {
  pll_hardware_probe();
  unsigned int arch = PLL_ATTRIB_ARCH_CPU;
  std::cout << "Libpll-2 will use ";
  if (pll_hardware.avx_present) {
    arch = PLL_ATTRIB_ARCH_AVX;
    std::cout << "avx";
  } else if (pll_hardware.sse_present) {
    arch = PLL_ATTRIB_ARCH_SSE;
    std::cout << "sse";
  } else {
    std::cout << "no";
  }
  std::cout << " vectorization" << std::endl;
  return PLL_ATTRIB_SITE_REPEATS | arch;
}

pllmod_treeinfo_t * LikelihoodEvaluator::buildTreeinfo(const bpp::TreeTemplate<Node> *bppTree)
{
  WHEREAMI( __FILE__ , __LINE__ );
  unsigned int partitions_number = 1; 

  // sequences 
  std::string fileName = fileNamePrefix + "alignment.fasta"; 
  const char* fasta_file = fileName.c_str();
  pll_sequences sequences;
  unsigned int *pattern_weights = read_from_fasta(fasta_file, sequences);

  // tree
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
}

void LikelihoodEvaluator::fullOptimizeTreeinfo(pllmod_treeinfo_t *treeinfo)
{
  if (method != LIBPLL2) {
    return;
  }

  double previousLogl = getTreeinfoLikelihood(treeinfo); 



  double newLogl = previousLogl;
  std::cout << "LikelihoodEvaluator::fullOptimizeTreeinfo before: ll = " 
    << getTreeinfoLikelihood(treeinfo) << std::endl;
  std::cout << "Genes in the tree: " << treeinfo->tip_count << std::endl;
  unsigned int it = 0;
  unsigned int max_it = 50;
  unsigned int tolerance = 0.5;
  do {
    previousLogl = newLogl;
    fullOptimizeTreeinfoIter(treeinfo);
    newLogl = getTreeinfoLikelihood(treeinfo);
    ++it;
  } while (it < max_it && newLogl - previousLogl > tolerance);
  std::cout << "LikelihoodEvaluator::fullOptimizeTreeinfo after: ll = " << newLogl <<
    " (" << it << " iterations)" << std::endl;
}
  
  
double LikelihoodEvaluator::optimizeTreeinfoLocal(pllmod_treeinfo_t *treeinfo)
{
  if (rollbacks_.size()) {
    rollbacks_.top()->localOptimize();
  }
  return getTreeinfoLikelihoodNoRecompute(currentTreeinfo);
}


double LikelihoodEvaluator::fullOptimizeTreeinfoIter(pllmod_treeinfo_t *treeinfo)
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


void LikelihoodEvaluator::initialize_LIBPLL2() 
{
  loadStrictNamesFromAlignment_forPLL();
  writeAlignmentFilesForPLL();
  if(logLikelihood == 0) {
    logLikelihood = evaluate(&tree) * scaler_;
  }   
}

void LikelihoodEvaluator::initialize_PLL()
{
  WHEREAMI( __FILE__ , __LINE__ );
  loadStrictNamesFromAlignment_forPLL();
  writeAlignmentFilesForPLL();
  alpha_ = 1.0;
  PLL_initializePLLInstance();
  PLL_loadAlignment(fileNamePrefix + "alignment.fasta");
  PLL_loadPartitions(fileNamePrefix + "partition.txt");
  
  if(logLikelihood == 0) {
    logLikelihood = evaluate(&tree) * scaler_;
  }  
}


void LikelihoodEvaluator::setTree(TreeTemplate<Node> * newTree)
{
  WHEREAMI( __FILE__ , __LINE__ );

  destroyRollbacks();
  if(!isInitialized()){
    if(tree)
      delete tree;
    tree = newTree->clone();
  }else
    throw Exception("This evaluator has already be initialized. It is forbidden to modify it now.");
}
  
pll_unode_t *LikelihoodEvaluator::getLibpllNode(bpp::Node *node)
{
  Clonable *prop = 0;
  try {
    prop = node->getNodeProperty("libpll");
  } catch (Exception e) {
    std::cout << "Error: cannot find libpll property for node " << node->getId() << std::endl;
    return 0;
  }
  if (!prop) {
    std::cout << "Error: cannot find libpll property for node " << node->getId() << std::endl;
    return 0;
  }
  unsigned int libpllid =  dynamic_cast<LibpllNodeProperty *>(prop)->id;
  return currentUtree->nodes[libpllid];
}

bool areEquals(pll_unode_t *n1, pll_unode_t *n2) 
{
  if (n2->next)
    return n1== n2 || n1 == n2->next || n1 == n2->next->next;
  else
    return n1 == n2;
}

bool sprYeldsSameTree(pll_unode_t *n1, pll_unode_t *n2)
{
  return (n2 == n1) || (n2 == n1->next) || (n2 == n1->next->next)
    || (n2->back == n1) || (n2->back == n1->next) || (n2->back == n1->next->next);
}

// returns the branch between n1 and n2 with the same index as n2
pll_unode_t *getBranchFromLibpll(pll_unode_t *n1, pll_unode_t *n2) {
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

// returns the branch between n1 and n2 with the same index as n2
pll_unode_t *LikelihoodEvaluator::getBranch(bpp::Node *n1, bpp::Node *n2) {
  return getBranchFromLibpll(getLibpllNode(n1), getLibpllNode(n2));
}

void LikelihoodEvaluator::applyNNI(bpp::TreeTemplate<bpp::Node> &tree,
    bpp::Node *bppSon)
{
  Node * bppParent = bppSon->getFather();
  Node * bppGrandParent = bppParent->getFather();
  Node * bppUncle = bppGrandParent->getSon(0);
  if (bppUncle == bppParent) {
    bppUncle = bppGrandParent->getSon(1);
  }
  Node * bppRoot = tree.getRootNode();
  return applyNNIIntern(bppParent, bppGrandParent, bppSon, bppUncle, bppRoot);
}
void LikelihoodEvaluator::applyNNIRoot(bpp::Node *bppParent,
  bpp::Node *bppGrandParent,
  bpp::Node *bppSon, bpp::Node *bppUncle)
{
  //std::cout << "** LikelihoodEvaluator::applyNNIRoot" << std::endl;
  pll_unode_t *parent, *son, *uncle;
  try {
    parent = getLibpllNode(bppParent);
    son = getLibpllNode(bppSon);
    uncle = getLibpllNode(bppUncle);
  }
  catch (Exception e) {
    std::cout << "Error in getBranch: " << e.what()<< std::endl;
    return;
  }
  double previousLikelihood = getTreeinfoLikelihoodNoRecompute(currentTreeinfo);
  pll_unode_t *edge = getBranchFromLibpll(uncle, parent);
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
  pushRollback(new NNIRootRollback(edge, son, t1, t2, previousLikelihood));
  pllmod_utree_set_length(son, t1/2.0 + t2);
  pllmod_utree_set_length(edge, t1/2.0);
}

  
void LikelihoodEvaluator::applyNNIIntern(bpp::Node *bppParent, 
    bpp::Node *bppGrandParent,
    bpp::Node *bppSon, bpp::Node *bppUncle,
    bpp::Node *bppRoot)
{
  if (method != LIBPLL2 && method != HYBRID) {
    return;
  }
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
  } catch (Exception e) {
    std::cout << "Error in applyNNI: " << e.what()<< std::endl;
    return;
  }
  double previousLikelihood = getTreeinfoLikelihoodNoRecompute(currentTreeinfo);
  pll_unode_t *edge = getBranchFromLibpll(grandParent, parent);
  if (!edge) { 
    std::cout << "LikelihoodEvaluator::applyNNI: Impossible to find the good NNI move" << std::endl;
    return;
  }
  bool sonNext = areEquals(edge->next->back, son);
  bool uncleNext = areEquals(edge->back->next->back, uncle);
  unsigned int move = (sonNext == uncleNext) ?
    PLL_UTREE_MOVE_NNI_LEFT : PLL_UTREE_MOVE_NNI_RIGHT;
  pll_tree_rollback_t rollbackInfo;
  if (!pllmod_utree_nni(edge, move, &rollbackInfo)) {
    std::cout << "failed applying nni : " << pll_errmsg << std::endl;
  }
  pushRollback(new NNIRollback(currentTreeinfo, rollbackInfo, previousLikelihood));
}

bpp::Node *getBrother(bpp::Node *node) 
{
  if (node->getFather()->getSon(0) != node) {
    return node->getFather()->getSon(0);
  } else {
    return node->getFather()->getSon(1);
  }
}

bpp::Node *getFatherOrBrotherIfRoot(bpp::Node *node)
{
  if (node->getFather()->hasFather()) {
    return node->getFather();
  } else {
    return getBrother(node);
  }
}

void LikelihoodEvaluator::applySPR(bpp::Node *bppToCut,
      bpp::Node *bppNewBrother)
{
  //std::cout << "LikelihoodEvaluator::applySPR " << std::endl;
  if (method != LIBPLL2 && method != HYBRID) {
    return;
  }
  double previousLikelihood = getTreeinfoLikelihoodNoRecompute(currentTreeinfo);
  //std::cout << "previousLikelihood " << previousLikelihood << std::endl;
  // get bpp nodes 
  bpp::Node *bppOldFather = getFatherOrBrotherIfRoot(bppToCut);
  bpp::Node *bppOldBrother = getBrother(bppToCut);
  
  // this one is the father of the new brother BEFORE spr
  bpp::Node *bppNewBrotherFather = getFatherOrBrotherIfRoot(bppNewBrother);

  // get libpll nodes
  pll_unode_t *toCut = getBranch(bppToCut, bppOldFather );
  pll_unode_t *newBrotherFather = getBranch(bppNewBrother, bppNewBrotherFather);
  pll_unode_t *oldBrother = getBranch(bppOldFather, bppOldBrother); 

  if (!toCut || !newBrotherFather) {
    std::cout << "Error, null pll_unode_t in applySPR " << toCut << " " << newBrotherFather << std::endl;
    pushRollback(new IdentityRollback(previousLikelihood));
    return;
  }

  if (sprYeldsSameTree(toCut, newBrotherFather)) { // will yeld same tree
    if (method == HYBRID) {
      if (!oldBrother) {
        std::cout << "Warning, null oldBrother. Next hybrid mismatch might be ok" << std::endl;
        pushRollback(new IdentityRollback(previousLikelihood));
      } else {
        pushRollback(new SPRIdentityHybridRollback(previousLikelihood,  
                          newBrotherFather, oldBrother));
      }
    } else {
      pushRollback(new IdentityRollback(previousLikelihood));
    }
  } else {
    pll_tree_rollback_t rb;
    if (pllmod_utree_spr(toCut, newBrotherFather, &rb)!= PLL_SUCCESS) {
      std::cout << "Error: cannot apply SPR" << std::endl;
      std::cout << pll_errmsg << std::endl;
      pushRollback(new IdentityRollback(previousLikelihood));
      return;
    }
    pushRollback(new SPRRollback(currentTreeinfo, rb, previousLikelihood));
  }
  if (method == HYBRID && oldBrother) { // be consistent with the bpp and pll trees
    // newBrotherFather->back is the father of toCut after SPR
    bool atRoot = !bppNewBrother->getFather()->hasFather();
    double bl = 0.1;
    if (atRoot) {
      bl += bppNewBrotherFather->getDistanceToFather();
    }
    pllmod_utree_set_length(newBrotherFather->back, bl);
    pllmod_utree_set_length(newBrotherFather->back->next, 0.1);
    pllmod_utree_set_length(newBrotherFather->back->next->next, 0.1);
    pllmod_utree_set_length(oldBrother, 0.1);
  }
  //std::cout << "After applySPR: " <<  printer.getTreeinfoString(currentTreeinfo, false, false) << std::endl;
}

void LikelihoodEvaluator::destroyRollbacks()
{
  while (rollbacks_.size()) {
    delete rollbacks_.top();
    rollbacks_.pop();
  }
}

void LikelihoodEvaluator::pushRollback(Rollback *rollback)
{
  rollbacks_.push(rollback);
}

void LikelihoodEvaluator::rollbackAllMoves()
{
  if (!rollbacks_.size()) {
    return;
  }
  while (rollbacks_.size()) {
    if (!rollbackLastMove()) {
      return;
    }
  }
}

bool LikelihoodEvaluator::rollbackLastMove()
{
  static int maxSize = 0;
  if (rollbacks_.size() > maxSize) {
    maxSize = rollbacks_.size();
    std::cout << "new rb max size " << maxSize << std::endl;
  }
  //std::cout << "LikelihoodEvaluator::rollbackLastMove" << std::endl;
  if (method != LIBPLL2 && method != HYBRID) {
    return false;
  }
  if (!rollbacks_.size()) {
    std::cout << "Error: no rollback info" << std::endl;
    return false;
  }
  if (!rollbacks_.top()->applyRollback()) {
    std::cout << "An error occured while trying to rollback libpll NNI move" << std::endl;
    return false;
  }
  double previousLikelihood = rollbacks_.top()->getPreviousLikelihood();
  delete rollbacks_.top();
  rollbacks_.pop();
  double newLL = getTreeinfoLikelihood(currentTreeinfo);
  if (fabs(newLL - previousLikelihood) > 0.001) {
    std::cout << "Error: different likelihood after rollback " << std::endl;
    std::cout << "  before: " << previousLikelihood << std::endl;
    std::cout << "  after : " << newLL << std::endl;
    return false;
  }
  return true;
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


double LikelihoodEvaluator::libpllEvaluateIterative(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)

{
  //std::cout << "LikelihoodEvaluator::libpllEvaluateIterative" << std::endl;

  if (!currentTreeinfo) {
    return libpllEvaluateFromScratch(treeToEvaluate);
  }
  
  //std::cout << "TreeInfo" <<  printer.getTreeinfoString(currentTreeinfo, true, false, false) << std::endl;
  //std::cout << "Tree BPP" <<  printer.getBPPNodeString((*treeToEvaluate)->getRootNode(), true, false) << std::endl;
 
  if (needFullOptim && method != HYBRID) {
    fullOptimizeTreeinfo(currentTreeinfo); 
    needFullOptim = false;
  } else if (method != HYBRID) {
    optimizeTreeinfoLocal(currentTreeinfo);
  }
  double result_ll = getTreeinfoLikelihood(currentTreeinfo);
  updateTreeToEvaluate(treeToEvaluate, currentTreeinfo, printer);
  mapUtreeToBPPTree(currentUtree, *treeToEvaluate);
  return result_ll;
}

double LikelihoodEvaluator::libpllEvaluateFromScratch(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)
{
  std::cout << "LikelihoodEvaluator::libpllEvaluateFromScratch" << std::endl;
  WHEREAMI( __FILE__ , __LINE__ );
  if (currentTreeinfo) {
    destroyTreeinfo();
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
 
  currentTreeinfo = buildTreeinfo(alternativeTree ? alternativeTree : tree);
  double result_ll = getTreeinfoLikelihood(currentTreeinfo, false);
  fullOptimizeTreeinfo(currentTreeinfo);
  needFullOptim = false;
  result_ll = getTreeinfoLikelihood(currentTreeinfo);

  std::string newStr = printer.getTreeinfoString(currentTreeinfo, true, false, false);
  stringstream outputNewickString;
  outputNewickString.str(newStr);
  Newick outputNewick;
  delete *treeToEvaluate;
  *treeToEvaluate = outputNewick.read(outputNewickString);
  if (wasRooted) {
    reroot(*treeToEvaluate, leaves1, leaves2);
  }
  mapUtreeToBPPTree(currentUtree, *treeToEvaluate);
  delete inputBPPTree;
 

  return result_ll;
}

double LikelihoodEvaluator::HYBRID_evaluate(TreeTemplate<Node>** treeToEvaluate)
{
  TreeTemplate<Node> *pllTree = (*treeToEvaluate)->clone();
  double pllRes = PLL_evaluate(&pllTree);
  delete pllTree;
  double libpllRes = LIBPLL2_evaluate(treeToEvaluate);
  if (fabs(pllRes - libpllRes) > 0.0001) {
    std::cout << "error: different results with PLL and libpll-2" << std::endl;
  }
  return libpllRes;

}

double LikelihoodEvaluator::LIBPLL2_evaluate(TreeTemplate<Node>** treeToEvaluate)
{
  if (!currentTreeinfo) {
    return libpllEvaluateFromScratch(treeToEvaluate);
  } else {
    return libpllEvaluateIterative(treeToEvaluate);
  }
}

 
double LikelihoodEvaluator::PLL_evaluate(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)
{
  WHEREAMI( __FILE__ , __LINE__ );

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

  // processing by PLL
  PLL_connectTreeAndAlignment();
  // pllSetFixedAlpha(alpha_, 0, PLL_partitions, PLL_instance);
  // pllSetFixedBaseFrequencies(baseFreq_, 4, 0, PLL_partitions, PLL_instance);
  // pllSetFixedSubstitutionMatrix(subsMatrix_, 6, 0, PLL_partitions, PLL_instance);
  WHEREAMI( __FILE__ , __LINE__ );
  if(!pll_model_already_initialized_){
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
 
  double result_ll = 0.0;
  char *newickStr = 0;
  
  result_ll = PLL_instance->likelihood;
  if (method != HYBRID) {
    pllOptimizeModelParameters(PLL_instance, PLL_partitions, tolerance_);
    result_ll = PLL_instance->likelihood;
  }
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

  delete treeForPLL;

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
  std::cout << "LikelihoodEvaluator::BPP_evaluate " << std::endl;

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





void LikelihoodEvaluator::unload()
{
  
 
  if(!initialized)
    return;
  initialized = false;
  
  if(method == PLL || method == HYBRID){
    if(pll_model_already_initialized_){
      pllAlignmentDataDestroy(PLL_alignmentData);
      pllPartitionsDestroy(PLL_instance, &PLL_partitions);
      pllDestroyInstance(PLL_instance);
      pll_model_already_initialized_ = false;
    }
  }
  
  if (method == BPP)
  {
    delete nniLk;
    if(nniLkAlternative)
      delete nniLkAlternative;
  }
}

LikelihoodEvaluator::~LikelihoodEvaluator()
{
  WHEREAMI( __FILE__ , __LINE__ );
  unload();
  destroyTreeinfo(); 
  while (rollbacks_.size()) {
    delete rollbacks_.top();
    rollbacks_.pop();
  }
  delete tree;
  if(alternativeTree)
    delete alternativeTree;
  if(method != BPP){
    remove(string(fileNamePrefix + "alignment.fasta").c_str());
    remove(string(fileNamePrefix + "partition.txt").c_str());
  }
  delete sites;
  delete alphabet;
  
  //delete substitutionModel;
  //delete rateDistribution;
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
  if(method == LIBPLL2) {
    initialize_LIBPLL2();
  }  else if (method == BPP) {
    initialize_BPP_nniLk();
  } else {
    initialize_PLL();
  } 

  //
  initialized = true;
}

double LikelihoodEvaluator::evaluate(bpp::TreeTemplate<bpp::Node>** treeToEvaluate)
{
  if (method == BPP) {
    return BPP_evaluate(treeToEvaluate);
  } else if (method == PLL) {
    return PLL_evaluate(treeToEvaluate);
  } else if (method == LIBPLL2) {
    return LIBPLL2_evaluate(treeToEvaluate);
  } else if (method == HYBRID) {
    return HYBRID_evaluate(treeToEvaluate);
  }
}

void LikelihoodEvaluator::setAlternativeTree(TreeTemplate< Node >* newAlternative)
{
  WHEREAMI( __FILE__ , __LINE__ );
  if(!initialized)
    initialize();
  if(alternativeTree != 00)
    delete alternativeTree;
  alternativeTree = newAlternative->clone();
  alternativeLogLikelihood = evaluate( &alternativeTree ) * scaler_;
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
  if(aligmentFilesForPllWritten_)
    return;
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
  this->tree = dynamic_cast<TreeTemplate<Node> *>(tree->clone());
  this->substitutionModel = model->clone();
  this->rateDistribution = rateDistribution->clone();
  this->sites = dynamic_cast<VectorSiteContainer*>(alignment->clone());

  string methodString = ApplicationTools::getStringParameter("likelihood.evaluator",params,"PLL");
  std::cout << " METHOD " << methodString << std::endl;
  if (methodString == "PLL") {
    this->method = PLL;
  } else if (methodString == "LIBPLL2") {
    this->method = LIBPLL2;
  } else if (methodString == "HYBRID") {
    this->method = HYBRID;
  } else {
    this->method = BPP;
  }
  initialize();

}
