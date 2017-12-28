




#ifndef _FASTRECONCILIATIONTOOLS_H_
#define _FASTRECONCILIATIONTOOLS_H_

#include <map>
#include <vector>
#include <set>

#include<Bpp/Phyl/Node.h>
#include<Bpp/Phyl/TreeTemplate.h>

using namespace bpp;

struct ReconciliationCache;

class FastReconciliationTools {

  public:
    double findMLReconciliationDR (TreeTemplate<Node> * speciesTree,
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


    static double computeLogBranchProbability ( const double & duplicationProbability, 
        const double & lossProbability, 
        const int numberOfLineages);

  private:

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
        std::vector <std::vector<int> > & dupData ,
        ReconciliationCache &cache);

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

    void computeNumbersOfLineagesFromRoot ( TreeTemplate<Node> * spTree,
        TreeTemplate<Node> * geneTree,
        Node * node,
        const std::map<std::string, std::string > &seqSp,
        const std::map<std::string, int > &spID,
        std::vector <int> &num0lineages,
        std::vector <int> &num1lineages,
        std::vector <int> &num2lineages,
        std::vector <std::vector<int> > & speciesIDs,
        std::vector <std::vector<int> > & dupData,
        std::set <int> & branchesWithDuplications );

    void computeNumbersOfLineagesFromRootIter ( TreeTemplate<Node> * spTree,
        TreeTemplate<Node> * geneTree,
        Node * node,
        const std::map<std::string, std::string > &seqSp,
        const std::map<std::string, int > &spID,
        std::vector <int> &num0lineages,
        std::vector <int> &num1lineages,
        std::vector <int> &num2lineages,
        std::vector <std::vector<int> > & speciesIDs,
        std::vector <std::vector<int> > & dupData,
        std::set <int> & branchesWithDuplications );

    int assignSpeciesIdToLeaf ( Node * node,  
        const std::map<std::string, std::string > & seqSp,
        const std::map<std::string, int > & spID );


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


    void recoverLossesAndLineages ( Node *& node, int & a, const int & b, int & olda, int & a0,
        const TreeTemplate<Node> & tree,
        int & dupData, std::vector<int> &num0lineages, std::vector<int> &num1lineages );

    void recoverLossesAndLineagesWithDuplication ( const Node * nodeA,
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

    static double computeBranchProbability ( const double & duplicationProbability, const double & lossProbability, const int numberOfLineages);

    // todobenoit : template and use stl better
    static void resetVector(std::vector<unsigned int> & v);
    static void resetVector(std::vector<int> & v);
    static void resetVector(std::vector<double> & v);



};


#endif
