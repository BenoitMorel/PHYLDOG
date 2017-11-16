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

/* This file contains various functions useful for reconciliations, such as reconciliation computation, printing of trees with integer indexes, search of a root with reconciliation...*/

#ifndef _RECONCILIATIONTOOLS_H_
#define _RECONCILIATIONTOOLS_H_

#include <queue>
#include <set>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/Nhx.h>


#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/NumConstants.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

// From Utils:
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Clonable.h>
#include <Bpp/Numeric/Number.h>
//#include <Bpp/Clonable.h>

#include <Bpp/BppString.h>
#include <Bpp/BppVector.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/TextTools.h>

#include "Constants.h"
#include "GeneTreeAlgorithms.h"



using namespace bpp;



void writeReconciledGeneTree ( map<string, string >& params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string>& seqSp, bool temporary ) ;
void writeReconciledGeneTreeToFile ( map<string, string >& params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string>& seqSp, string& reconcTree );
void assignArbitraryBranchLengths ( TreeTemplate<Node> & tree );
void reNumber ( TreeTemplate<Node> & tree, Node * noeud, int & index );
void reNumber ( TreeTemplate<Node> & tree );
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree );
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities );
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities, std::vector <double> & coalBl );
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector <double> & coalBl );
std::map <int, std::vector <int> > breadthFirstreNumberAndResetProperties ( TreeTemplate<Node> & tree );

void printVector ( std::vector<int> & v );
void resetLossesAndDuplications ( TreeTemplate<Node> & tree, /*std::vector <int> &lossNumbers, */std::vector <double> &lossProbabilities, /*std::vector <int> &duplicationNumbers, */std::vector <double> &duplicationProbabilities );
void resetLossesDuplicationsSpeciations ( TreeTemplate<Node> & tree, std::vector <int> &lossNumbers, std::vector <double> &lossProbabilities, std::vector <int> &duplicationNumbers, std::vector <double> &duplicationProbabilities, std::vector <int> &branchNumbers );
void resetLossesDuplicationsSpeciationsForGivenNodes ( TreeTemplate<Node> & tree, std::vector <int> & lossNumbers, std::vector <double> & lossProbabilities, std::vector <int> & duplicationNumbers, std::vector <double> & duplicationProbabilities, std::vector <int> & branchNumbers, std::vector <int> nodesToUpdate, std::map <int, std::vector<int> > & geneNodeIdToLosses, std::map <int, int > & geneNodeIdToDuplications, std::map <int, std::vector<int> > & geneNodeIdToSpeciations );
void resetVector ( std::vector<unsigned int> & v );
void resetVector ( std::vector<int> & v );
void resetVector ( std::vector<double> & v );
void resetSpeciesIds ( TreeTemplate<Node> & tree );
void resetSpeciesIdsForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds );
void resetSpeciesIdsForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate );
//void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::vector<int> &branchNumbers, std::map <int,int> & geneNodeIdToDuplications, std::map <int, std::vector <int> > & geneNodeIdToLosses, std::map <int, std::vector <int> > & geneNodeIdToSpeciations) ;
//void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::map <int,int> &geneNodeIdToDuplications, std::map <int, std::vector <int> > &geneNodeIdToLosses, std::map <int, std::vector <int> > &geneNodeIdToSpeciations) ;
void computeDuplicationAndLossProbabilities ( const int i, const int j, const int k, double & lossProbability, double & duplicationProbability );
void computeDuplicationAndLossProbabilitiesForAllBranches ( const std::vector <int> &numOGenes, const std::vector <int> &num1Genes, const std::vector <int> &num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities );
void computeAverageDuplicationAndLossProbabilitiesForAllBranches ( const std::vector <int> &numOGenes, const std::vector <int> &num1Genes, const std::vector <int> &num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities );
double computeBranchProbability ( const double &duplicationProbability, const double &lossProbability, const int numberOfLineages );
double computeLogBranchProbability ( const double &duplicationProbability, const double &lossProbability, const int numberOfLineages );
double computeBranchProbabilityAtRoot ( const double &duplicationProbability, double &lossProbability, const int numberOfLineages );
double computeLogBranchProbabilityAtRoot ( const double &duplicationProbability, const double &lossProbability, const int numberOfLineages );
void computeScenarioScore ( TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::vector<int> &branchNumbers, std::map <int, std::vector <int> > &geneNodeIdToSpeciations, const std::vector<double> &duplicationProbabilities, const std::vector <double> &lossProbabilities, std::vector <int> &num0lineages, std::vector <int> &num1lineages, std::vector <int> &num2lineages );
std::string nodeToParenthesisWithIntNodeValues ( const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName ) throw ( NodeNotFoundException );
std::string treeToParenthesisWithIntNodeValues ( const Tree & tree, bool bootstrap, const std::string & propertyName );
std::string nodeToParenthesisWithDoubleNodeValues ( const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName ) throw ( NodeNotFoundException );
std::string treeToParenthesisWithDoubleNodeValues ( const Tree & tree, bool bootstrap, const std::string & propertyName );
void setLossesAndDuplications ( TreeTemplate<Node> & tree,
                                std::vector <int> &lossNumbers,
                                std::vector <int> &duplicationNumbers );
void assignNumLineagesOnSpeciesTree ( TreeTemplate<Node> & tree,
                                      std::vector <int> &num0Lineages,
                                      std::vector <int> &num1Lineages,
                                      std::vector <int> &num2Lineages );
double makeReconciliationAtGivenRoot ( TreeTemplate<Node> * tree,
                                       TreeTemplate<Node> * geneTree,
                                       const std::map<std::string, std::string > seqSp,
                                       std::vector< double> lossProbabilities,
                                       std::vector < double> duplicationProbabilities,
                                       int MLindex );
/*double findMLReconciliation (TreeTemplate<Node> * spTree,
                             TreeTemplate<Node> * geneTreeSafe,
                             std::map<std::string, std::string > seqSp,
                             std::vector<int> & lossNumbers,
                             std::vector< double> lossProbabilities,
                             std::vector< int> & duplicationNumbers,
                             std::vector < double> duplicationProbabilities,
                             int & MLindex,
                             std::vector<int> &branchNumbers,
                             int speciesIdLimitForRootPosition,
                             std::vector <int> &num0lineages,
                             std::vector <int> &num1lineages,
                             std::vector <int> &num2lineages,
                             std::set <int> &nodesToTryInNNISearch);*/
int assignSpeciesIdToLeaf ( Node * node,
                            const std::map<std::string,
                            std::string > & seqSp,
                            const std::map<std::string, int > & spID );
struct ReconciliationCache;

void recoverLosses(Node *& node, int & a, const int & b, int & olda, const int & a0,
                   const TreeTemplate<Node> & tree,
                   double & likelihoodCell,
                   const std::vector< double> & lossRates,
                   const std::vector< double> & duplicationRates,
                   ReconciliationCache &cache);

void recoverLossesWithDuplication ( const Node * nodeA,
                                    const int &a,
                                    const int &olda,
                                    const TreeTemplate<Node> & tree,
                                    double & likelihoodCell,
                                    const std::vector< double> & lossRates,
                                    const std::vector< double> & duplicationRates );
double computeConditionalLikelihoodAndAssignSpId ( TreeTemplate<Node> & tree,
        std::vector <Node *> sons,
        double & rootLikelihood,
        double & son0Likelihood,
        double & son1Likelihood,
        const std::vector< double> & lossRates,
        const std::vector< double> & duplicationRates,
        int & rootSpId,
        const int & son0SpId,
        const int & son1SpId,
        int & rootDupData,
        int & son0DupData,
        int & son1DupData,
        bool atRoot,
        ReconciliationCache &cache);
double computeSubtreeLikelihoodPostorder ( TreeTemplate<Node> & spTree,
        TreeTemplate<Node> & geneTree,
        Node * node,
        const std::map<std::string, std::string > & seqSp,
        const std::map<std::string, int > & spID,
        std::vector <std::vector<double> > & likelihoodData,
        const std::vector< double> & lossRates,
        const std::vector < double> & duplicationRates,
        std::vector <std::vector<int> > & speciesIDs,
        std::vector <std::vector<int> > & dupData );
double computeSubtreeLikelihoodPostorderIter ( TreeTemplate<Node> & spTree,
        TreeTemplate<Node> & geneTree,
        Node * node,
        const std::map<std::string, std::string > & seqSp,
        const std::map<std::string, int > & spID,
        std::vector <std::vector<double> > & likelihoodData,
        const std::vector< double> & lossRates,
        const std::vector < double> & duplicationRates,
        std::vector <std::vector<int> > & speciesIDs,
        std::vector <std::vector<int> > & dupData,
        ReconciliationCache &cache);
void computeRootingLikelihood ( TreeTemplate<Node> & spTree,
                                Node * node,
                                std::vector <std::vector<double> > & likelihoodData,
                                const std::vector< double> & lossRates,
                                const std::vector < double> & duplicationRates,
                                std::vector <std::vector<int> > & speciesIDs,
                                std::vector <std::vector<int> > & dupData,
                                int sonNumber,
                                std::map <double, Node*> & LksToNodes,
                                ReconciliationCache &cache);
void computeSubtreeLikelihoodPreorder ( TreeTemplate<Node> & spTree,
                                        TreeTemplate<Node> & geneTree,
                                        Node * node,
                                        const std::map<std::string, std::string > & seqSp,
                                        const std::map<std::string, int > & spID,
                                        std::vector <std::vector<double> > & likelihoodData,
                                        const std::vector< double> & lossRates,
                                        const std::vector < double> & duplicationRates,
                                        std::vector <std::vector<int> > & speciesIDs,
                                        std::vector <std::vector<int> > & dupData,
                                        int sonNumber,
                                        std::map <double, Node*> & LksToNodes );
void computeSubtreeLikelihoodPreorderIter ( TreeTemplate<Node> & spTree,
                                        TreeTemplate<Node> & geneTree,
                                        Node * node,
                                        const std::map<std::string, std::string > & seqSp,
                                        const std::map<std::string, int > & spID,
                                        std::vector <std::vector<double> > & likelihoodData,
                                        const std::vector< double> & lossRates,
                                        const std::vector < double> & duplicationRates,
                                        std::vector <std::vector<int> > & speciesIDs,
                                        std::vector <std::vector<int> > & dupData,
                                        int sonNumber,
                                        std::map <double, Node*> & LksToNodes, 
                                        ReconciliationCache &cache);
void recoverLossesAndLineages ( Node *& node, int & a, const int & b, int & olda,
                                int & a0,
                                const TreeTemplate<Node> & tree,
                                int & dupData,
                                std::vector<int> &num0lineages,
                                std::vector<int> &num1lineages );
void recoverLossesAndLineagesWithDuplication ( const Node *& nodeA,
        const int &a,
        const int &olda,
        const TreeTemplate<Node> & tree,
        std::vector <int> &num0lineages );
void computeNumbersOfLineagesInASubtree ( TreeTemplate<Node> & tree,
        std::vector <Node *> sons,
        int & rootSpId,
        const int & son0SpId,
        const int & son1SpId,
        int & rootDupData,
        int & son0DupData,
        int & son1DupData,
        bool atRoot,
        std::vector <int> &num0lineages,
        std::vector <int> &num1lineages,
        std::vector <int> &num2lineages,
        std::set <int> &branchesWithDuplications );
void computeNumbersOfLineagesFromRoot ( TreeTemplate<Node> * spTree,
                                        TreeTemplate<Node> * geneTree,
                                        Node * node,
                                        const std::map<std::string, std::string > seqSp,
                                        const std::map<std::string, int > spID,
                                        std::vector <int> &num0lineages,
                                        std::vector <int> &num1lineages,
                                        std::vector <int> &num2lineages,
                                        std::vector <std::vector<int> > & speciesIDs,
                                        std::vector <std::vector<int> > & dupData,
                                        std::set <int> & branchesWithDuplications );
void computeNumbersOfLineagesFromRootIter ( TreeTemplate<Node> * spTree,
                                        TreeTemplate<Node> * geneTree,
                                        Node * node,
                                        const std::map<std::string, std::string > seqSp,
                                        const std::map<std::string, int > spID,
                                        std::vector <int> &num0lineages,
                                        std::vector <int> &num1lineages,
                                        std::vector <int> &num2lineages,
                                        std::vector <std::vector<int> > & speciesIDs,
                                        std::vector <std::vector<int> > & dupData,
                                        std::set <int> & branchesWithDuplications );
double findMLReconciliationDR ( TreeTemplate<Node> * spTree,
                                TreeTemplate<Node> * geneTree,
                                const std::map<std::string, std::string > seqSp,
                                const std::map<std::string, int > spID,
                                std::vector< double> lossRates,
                                std::vector < double> duplicationRates,
                                int & MLindex,
                                std::vector <int> &num0lineages,
                                std::vector <int> &num1lineages,
                                std::vector <int> &num2lineages,
                                std::set <int> &nodesToTryInNNISearch,
                                const bool fillTables = true );
double computeAverageLossProportionOnCompletelySequencedLineages ( const std::vector <int> & num0lineages,
        const std::vector <int> & num1lineages,
        const std::vector <int> & num2lineages,
        const std::map <std::string, int> & genomeMissing,
        const TreeTemplate<Node> & tree );
/*void alterLineageCountsWithCoverages(std::vector <int> & num0lineages,
                                     std::vector <int> & num1lineages,
                                     std::vector <int> & num2lineages,
                                     std::map <std::string, int> & genomeMissing,
                                     TreeTemplate<Node> & tree);*/
void alterLineageCountsWithCoverages ( std::vector <int> & num0lineages,
                                       std::vector <int> & num1lineages,
                                       std::vector <int> & num2lineages,
                                       const std::map <std::string, int> & genomeMissing,
                                       const TreeTemplate<Node> & tree,
                                       bool average );
/*void alterLineageCountsWithCoveragesAverage(std::vector <int> & num0lineages,
                                            std::vector <int> & num1lineages,
                                            std::vector <int> & num2lineages,
                                            std::map <std::string, int> & genomeMissing,
                                            TreeTemplate<Node> & tree);*/
void alterLineageCountsWithCoveragesInitially ( std::vector <int> & num0lineages,
        std::vector <int> & num1lineages,
        std::vector <int> & num2lineages,
        std::map <std::string, int> & genomeMissing,
        TreeTemplate<Node> & tree );
void resetLineageCounts ( std::vector <int> & num0lineages,
                          std::vector <int> & num1lineages,
                          std::vector <int> & num2lineages );
void deleteSubtreeProperties ( Node &node );
void deleteTreeProperties ( TreeTemplate<Node> & tree );
std::map <std::string, int> computeSpeciesNamesToIdsMap ( TreeTemplate<Node> & tree );
void computeDuplicationAndLossRatesForTheSpeciesTree ( std::string &branchProbaOptimization,
        std::vector <int> & num0Lineages,
        std::vector <int> & num1Lineages,
        std::vector <int> & num2Lineages,
        std::vector<double> & lossExpectedNumbers,
        std::vector<double> & duplicationExpectedNumbers,
        std::map <std::string, int> & genomeMissing,
        TreeTemplate<Node> & tree );
void computeDuplicationAndLossRatesForTheSpeciesTreeInitially ( std::string &branchProbaOptimization,
        std::vector <int> & num0Lineages,
        std::vector <int> & num1Lineages,
        std::vector <int> & num2Lineages,
        std::vector<double> & lossProbabilities,
        std::vector<double> & duplicationProbabilities,
        const std::map <std::string, int> & genomeMissing,
        const TreeTemplate<Node> & tree );


/**
  * Cleanly remove a leaf from a tree, to get a propre tree
  */

void removeLeaf ( TreeTemplate<Node> & tree, std::string toRemove );
void cleanVectorOfOptions ( std::vector<std::string> & listOptions, bool sizeConstraint );
bool sortMaxFunction ( std::pair <std::string, double> i, std::pair <std::string, double> j );
bool sortMinFunction ( std::pair <std::vector<std::string>, double> i, std::pair <std::vector<std::string>, double> j );
void generateListOfOptionsPerClient ( std::vector <std::string> listOptions,
                                      int size,
                                      std::vector <std::vector<std::string> > &listOfOptionsPerClient,
                                      std::vector <unsigned int> &numbersOfGenesPerClient );
std::string removeComments (
    const std::string & s,
    const std::string & begin,
    const std::string & end );
void annotateGeneTreeWithDuplicationEvents ( TreeTemplate<Node> & spTree,
        TreeTemplate<Node> & geneTree,
        Node * node,
        const std::map<std::string, std::string > seqSp,
        const std::map<std::string, int > spID );
void annotateGeneTreeWithScoredDuplicationEvents ( TreeTemplate<Node> & spTree,
        TreeTemplate<Node> & geneTree,
        Node * node,
        const std::map<std::string, std::string > seqSp,
        const std::map<std::string, int > spID );

/*void editDuplicationNodesMuffato2(TreeTemplate<Node> & spTree,
                 TreeTemplate<Node> & geneTree,
                 Node * node,
                 double editionThreshold) ;
void editDuplicationNodesMuffato(TreeTemplate<Node> & spTree,
                 TreeTemplate<Node> & geneTree,
                 Node * node,
                 double editionThreshold) ;
*/

//To sort in descending order
bool cmp ( int a, int b );

//To sort in ascending order
bool anticmp ( int a, int b ) ;



bpp::Alphabet* getAlphabetFromOptions ( std::map <std::string, std::string>  params, bool& cont);

//reads an alignment from a file and parses it according to options
bpp::VectorSiteContainer *  getSequencesFromOptions ( std::map <std::string, std::string>  params, bpp::Alphabet* alphabet, bool& cont ) ;

//constructs a substitution model according to options
bpp::SubstitutionModel*   getModelFromOptions ( std::map <std::string, std::string> params, bpp::Alphabet *alphabet, bpp::VectorSiteContainer *sites, bool& cont ) ;

//constructs a rate heterogeneity distribution according to options
bpp::DiscreteDistribution* getRateDistributionFromOptions ( std::map <std::string,std::string> params, SubstitutionModel* model, bool& cont ) ;

//reads a tree from a file or builds it according to options
bpp::TreeTemplate<bpp::Node> * getTreeFromOptions ( std::map <std::string,std::string> params, bpp::Alphabet *alphabet, bpp::VectorSiteContainer * sites, bpp::SubstitutionModel* model, bpp::DiscreteDistribution* rDist, bool& cont ) ;

//Gets the map between sequence names and species names from options
void getCorrespondanceSequenceSpeciesFromOptions ( std::map< std::string, std::string > params, bool& cont, std::map <std::string, std::string> &seqSp, std::map<std::string, std::deque<std::string> > &spSeq );
//Remove sequences from the alignment, e.g. if they come from species we are not interested in
void removeUselessSequencesFromAlignment ( const TreeTemplate<Node>* spTree, bpp::VectorSiteContainer * sites, bool& cont, std::map<std::string, std::deque<std::string> > &spSeq, std::string file);

//Remove sequences from both the alignment and gene tree in case they have a very long branch
void qualityControlGeneTree ( bpp::TreeTemplate<bpp::Node>* geneTree, bpp::VectorSiteContainer * sites, bool& cont, std::string file);

// Utilitary function to recover a loss after a duplication.
int recoverLossClosestToDuplication(TreeTemplate<Node> * spTree, int a, int aold);

//Set node property with the node ID
void setNDProperty(TreeTemplate<Node> * tree);

//Recovers the events of duplication and loss for a given gene tree wrt a species tree.
vector < std::string > recoverDuplicationsAndLosses(
                    TreeTemplate<Node> * cTree,
                    TreeTemplate<Node> * spTree,
                    const string &familyName);

//Wrapper that annotates a gene tree then calls recoverDuplicationsAndLosses and finally writes it to an events file.
void outputNumbersOfEventsPerFamilyPerSpecies( map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, bool temporary );
void outputNumbersOfEventsToFile( map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, std::string& evFile ) ;

// Recovers the events of duplication for a given gene tree wrt a species tree,
// so that we can tell which genes are orthologous and which genes are paralogous.
vector < std::string > recoverOrthologsAndParalogs(
                    TreeTemplate<Node> * cTree,
                    TreeTemplate<Node> * spTree,
                    const string &familyName);

// Outputs to file orthologous and paralogous genes
void outputOrthologousAndParalogousGenes(map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, bool temporary ) ;
void outputOrthologousAndParalogousGenesToFile(map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, std::string& orFile ) ;


#endif  //_RECONCILIATIONTOOLS_H_
