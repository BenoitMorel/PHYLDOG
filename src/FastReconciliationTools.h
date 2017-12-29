




#ifndef _FASTRECONCILIATIONTOOLS_H_
#define _FASTRECONCILIATIONTOOLS_H_

#include <map>
#include <vector>
#include <set>

#include<Bpp/Phyl/Node.h>
#include<Bpp/Phyl/TreeTemplate.h>

using namespace bpp;


class FastReconciliationTools {

  public:
    FastReconciliationTools(TreeTemplate<Node> * speciesTree,
        TreeTemplate<Node> * geneTree,
        const std::map<std::string, std::string > seqSp,
        const std::map<std::string, int > spID,
        std::vector< double> lossRates,
        std::vector < double> duplicationRates,
        std::vector <int> &num0lineages,
        std::vector <int> &num1lineages,
        std::vector <int> &num2lineages,
        std::set <int> &nodesToTryInNNISearch,
        const bool fillTables = true );

    ~FastReconciliationTools();
    
    double findMLReconciliationDR(int & MLindex);

  double computeLogBranchProbability (int branch, int numberOfLineages );

  private:

    void computeSubtreeLikelihoodPreorder (Node * node,
        int sonNumber);

    void computeSubtreeLikelihoodPreorderIter (Node * node,
        int sonNumber);

    void computeRootingLikelihood (Node * node,
        int sonNumber);

    double computeSubtreeLikelihoodPostorder (Node * node);

    double computeSubtreeLikelihoodPostorderIter (Node *node);

    double computeConditionalLikelihoodAndAssignSpId (double & rootLikelihood,
        double & son0Likelihood,
        double & son1Likelihood,
        int & rootSpId,
        int son0SpId,
        int son1SpId,
        int & rootDupData,
        int & son0DupData,
        int & son1DupData,
        bool atRoot);

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

    void recoverLosses(Node *& node, 
        int &a, 
        int b, 
        int &olda, 
        int a0,
        double &likelihoodCell);

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
   

    bool isDescendant(Node *father, int descendantId);

    /* precompute stuff */
    void initialize();
private:
    TreeTemplate<Node> &_speciesTree;
    TreeTemplate<Node> &_geneTree;
    std::vector<Node *> _speciesNodes;
    std::map<std::string, std::string > _seqSp;
    std::map<std::string, int > _spID;
    std::vector<double> _lossRates;
    std::vector<double> _duplicationRates;
    std::vector <std::vector<int> > _speciesIDs;
    std::vector <std::vector<int> > _dupData;
    //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
    std::map <double, Node*> _LksToNodes;
    std::vector <std::vector<double> > _likelihoodData;
    bool _fillTables;
    
    // cache
    std::vector<int> _speciesIdsPreorder;
    std::vector<int> _speciesIdsLastSon;
    std::vector<std::vector< double > > _logBranchProbabilities;

    std::vector<int> &_num0lineages;
    std::vector<int> &_num1lineages;
    std::vector<int> &_num2lineages;
    std::set<int> &_nodesToTryInNNISearch;
    

};


#endif
