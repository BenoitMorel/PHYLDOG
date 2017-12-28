#include "FastReconciliationTools.h"    
#include "Constants.h"

using namespace std;

struct ReconciliationCache {

  ReconciliationCache(unsigned int spMaxId):
    nodesIds(spMaxId + 1), logBranchProbabilities(3) 
  {  
    logBranchProbabilities[0] = std::vector<double>(spMaxId + 1, 0.0);
    logBranchProbabilities[1] = std::vector<double>(spMaxId + 1, 0.0);
    logBranchProbabilities[2] = std::vector<double>(spMaxId + 1, 0.0);
  }

  bool contains(Node *node, int element) {
    unsigned int id = node->getId();
    if (!nodesIds[id].size()) {
      TreeTemplateTools::getNodesId(*node, nodesIds[id]);
    }
    return VectorTools::contains(nodesIds[id], element);
  }

  double computeLogBranchProbabilityCached(unsigned int branch, 
      unsigned int numberOfLineages, 
      const std::vector<double> &duplicationsProbabilities, 
      const std::vector<double> &lossProbabilities)
  {
    double res = logBranchProbabilities[numberOfLineages][branch];
    if (res == 0.0) {
      res = FastReconciliationTools::computeLogBranchProbability(duplicationsProbabilities[branch], lossProbabilities[branch], numberOfLineages);
      logBranchProbabilities[numberOfLineages][branch] = res;
    }
    return res;
  }

  std::vector<std::vector< int > > nodesIds;
  std::vector<std::vector< double > > logBranchProbabilities;
};


double FastReconciliationTools::findMLReconciliationDR (TreeTemplate<Node> * spTree,
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
        const bool fillTables)
{
  assert(speciesTree->isRooted());
  assert(geneTreeTree->isRooted());
  
  std::vector <double> nodeData ( 3, 0.0 );
  std::vector <std::vector<double> > likelihoodData ( geneTree->getNumberOfNodes(), nodeData );
  std::vector <int> nodeSpId ( 3, 0 );
  std::vector <std::vector<int> > speciesIDs ( geneTree->getNumberOfNodes(), nodeSpId );
  std::vector <std::vector<int> > dupData = speciesIDs;
  //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
  std::map <double, Node*> LksToNodes;

  Node * geneRoot = geneTree->getRootNode();
  double initialLikelihood = computeSubtreeLikelihoodPostorder ( *spTree, *geneTree,
      geneRoot, seqSp, spID,
      likelihoodData, lossRates,
      duplicationRates, speciesIDs, dupData );

  std::vector <Node *> sons = geneRoot->getSons();

  if ( sons.size() !=2 ) {
    std::cerr <<"Error: "<<sons.size() << "sons at the root!"<<std::endl;
  }
  LksToNodes[initialLikelihood] = sons[0];
  //We fill the likelihood and species ID data for the root node.
  //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
  //with other nodes.
  likelihoodData[geneRoot->getId()][1] = likelihoodData[geneRoot->getSon ( 1 )->getId()][0];
  likelihoodData[geneRoot->getId()][2] = likelihoodData[geneRoot->getSon ( 0 )->getId()][0];
  speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon ( 1 )->getId()][0];
  speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon ( 0 )->getId()][0];
  dupData[geneRoot->getId()][1] = dupData[geneRoot->getSon ( 1 )->getId()][0];
  dupData[geneRoot->getId()][2] = dupData[geneRoot->getSon ( 0 )->getId()][0];

  for ( unsigned int i = 0; i< sons.size(); i++ ) {
    for ( unsigned int j =0; j<sons[i]->getNumberOfSons(); j++ ) {
      computeSubtreeLikelihoodPreorder ( *spTree, *geneTree,
          sons[i], seqSp, spID,
          likelihoodData,
          lossRates, duplicationRates,
          speciesIDs, dupData, j, LksToNodes );
    }
  }
  vector<Node*> nodes = geneTree->getNodes();
  for ( unsigned int i = 0 ; i < nodes.size() ; i++ ) {
    if ( nodes[i]->hasNodeProperty ( "outgroupNode" ) ) {
      nodes[i]->deleteNodeProperty ( "outgroupNode" );
      break;
    }
  }

  LksToNodes.rbegin()->second->setNodeProperty ( "outgroupNode", BppString ( "here" ) );

  if ( fillTables ) {
    //Now the best root has been found. I can thus run a function with this best root to fill all the needed tables. This additional tree traversal could be avoided.
    //To this end, the needed tables should be filled by the postfix and prefix traversals. This has not been done yet.
    //resetting
    speciesIDs = std::vector<std::vector<int> > ( geneTree->getNumberOfNodes(), nodeSpId );
    dupData = speciesIDs;
    // Getting a well-rooted tree
    TreeTemplate<Node > * tree = geneTree->clone();
    tree->newOutGroup ( LksToNodes.rbegin()->second->getId() );
    nodesToTryInNNISearch.clear();
    //Resetting numLineages std::vectors
    resetVector ( num0lineages );
    resetVector ( num1lineages );
    resetVector ( num2lineages );
    computeNumbersOfLineagesFromRoot ( spTree, tree,
        tree->getRootNode(),
        seqSp, spID,
        num0lineages, num1lineages,
        num2lineages, speciesIDs,
        dupData, nodesToTryInNNISearch );
    delete tree;
  }
  
  //We return the best likelihood
  MLindex = LksToNodes.rbegin()->second->getId();
  return LksToNodes.rbegin()->first;
}


/*****************************************************************************
 * This function computes the lower conditional likelihood of a subtree and
 * assigns its summit node a species ID. Notations are influenced by
 * Zmasek and Eddy algorithm (2001).
 *
 ****************************************************************************/
double FastReconciliationTools::computeConditionalLikelihoodAndAssignSpId ( TreeTemplate<Node> & tree,
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
    ReconciliationCache &cache)
{
  if ( rootLikelihood == 0.0 ) {
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;

    Node * temp0 = tree.getNode ( son0SpId );
    Node * temp1 = tree.getNode ( son1SpId );


    double oldLL = rootLikelihood;
    while ( a!=b ) { //There have been losses !
      if ( a>b ) {
        recoverLosses( temp0, a, b, olda, a0, tree, rootLikelihood, lossRates, duplicationRates, cache);
      }
      else {
        recoverLosses( temp1, b, a, oldb, b0, tree, rootLikelihood, lossRates, duplicationRates, cache);
      }
    }
    rootSpId = a;
    if ( ( a==a0 ) || ( b==b0 ) ) { //There has been a duplication !
      if ( ( a==a0 ) && ( b==b0 ) ) {
        rootDupData += son0DupData+son1DupData;
        rootLikelihood-= ( computeLogBranchProbability ( duplicationRates[a0], lossRates[a0], son0DupData ) +
            computeLogBranchProbability ( duplicationRates[b0], lossRates[b0], son1DupData ) );
      }//there has been no loss, here
      else if ( b==b0 ) { //The loss has occured before a0
        rootDupData += son1DupData+1;
        rootLikelihood-=computeLogBranchProbability ( duplicationRates[b0], lossRates[b0], son1DupData );
        recoverLossesWithDuplication ( temp0, a, olda, tree, rootLikelihood, lossRates, duplicationRates );
      }
      else { //The loss has occured before b0
        rootDupData += son0DupData+1;
        rootLikelihood-=computeLogBranchProbability ( duplicationRates[a0], lossRates[a0], son0DupData );
        recoverLossesWithDuplication ( temp1, b, oldb, tree, rootLikelihood, lossRates, duplicationRates );
      }
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
      }
      else {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
      }
    }
    else { //there was no duplication
      rootDupData = 1;
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
      }
      else {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
      }
    }
    //Setting the lower conditional likelihood for the node of interest.
    rootLikelihood += son0Likelihood + son1Likelihood;
  }
  return ( rootLikelihood );
}




/*****************************************************************************
 * This function performs a postorder tree traversal in order to find
 * likelihoods for rootings.
 * When followed by the preorder tree traversal function,
 * likelihoods for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 *
 ****************************************************************************/
double FastReconciliationTools::computeSubtreeLikelihoodPostorder ( TreeTemplate<Node> & spTree,
    TreeTemplate<Node> & geneTree,
    Node * node,
    const std::map<std::string, std::string > & seqSp,
    const std::map<std::string, int > & spID,
    std::vector <std::vector<double> > & likelihoodData,
    const std::vector< double> & lossRates,
    const std::vector < double> & duplicationRates,
    std::vector <std::vector<int> > & speciesIDs,
    std::vector <std::vector<int> > & dupData )
{
  ReconciliationCache cache(TreeTools::getMaxId(spTree, spTree.getRootId()));
  return computeSubtreeLikelihoodPostorderIter(spTree, geneTree, node, seqSp, spID,
      likelihoodData, lossRates, duplicationRates, speciesIDs, dupData, cache);
}

double FastReconciliationTools::computeSubtreeLikelihoodPostorderIter ( TreeTemplate<Node> & spTree,
    TreeTemplate<Node> & geneTree,
    Node * node,
    const std::map<std::string, std::string > & seqSp,
    const std::map<std::string, int > & spID,
    std::vector <std::vector<double> > & likelihoodData,
    const std::vector< double> & lossRates,
    const std::vector < double> & duplicationRates,
    std::vector <std::vector<int> > & speciesIDs,
    std::vector <std::vector<int> > & dupData ,
    ReconciliationCache &cache)
{
  int id=node->getId();
  if ( node->isLeaf() ) {
    if ( likelihoodData[id][0]==0.0 ) {
      speciesIDs[id][0]=speciesIDs[id][1]=speciesIDs[id][2]=assignSpeciesIdToLeaf ( node, seqSp, spID );
      likelihoodData[id][0]=likelihoodData[id][1]=likelihoodData[id][2]=computeLogBranchProbability ( duplicationRates[speciesIDs[id][0]], lossRates[speciesIDs[id][0]], 1 );
      dupData[id][0] = dupData[id][1] = dupData[id][2] = 1;
    }
    return ( likelihoodData[id][0] );
  }
  else {
    std::vector <Node *> sons = node->getSons();
    for ( unsigned int i = 0; i< sons.size(); i++ ) {
      computeSubtreeLikelihoodPostorderIter ( spTree, geneTree,
          sons[i], seqSp,
          spID, likelihoodData,
          lossRates, duplicationRates,
          speciesIDs, dupData , cache);

    }
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0 = 0;
    unsigned int directionSon1 = 0;

    computeConditionalLikelihoodAndAssignSpId ( spTree, sons,
        likelihoodData[id][0],
        likelihoodData[idSon0][directionSon0],
        likelihoodData[idSon1][directionSon1],
        lossRates, duplicationRates,
        speciesIDs[id][0],
        speciesIDs[idSon0][directionSon0],
        speciesIDs[idSon1][directionSon1],
        dupData[id][0],
        dupData[idSon0][directionSon0],
        dupData[idSon1][directionSon1],
        TreeTemplateTools::isRoot ( *node ),
        cache);
    return ( likelihoodData[id][0] );
  }
}


/*****************************************************************************
 * This function computes the likelihood of a rooting.
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 *
 ****************************************************************************/
void FastReconciliationTools::computeRootingLikelihood ( TreeTemplate<Node> & spTree,
    Node * node,
    std::vector <std::vector<double> > & likelihoodData,
    const std::vector< double> & lossRates,
    const std::vector < double> & duplicationRates,
    std::vector <std::vector<int> > & speciesIDs,
    std::vector <std::vector<int> > & dupData,
    int sonNumber,
    std::map <double, Node*> & LksToNodes,
    ReconciliationCache &cache)
{
  int geneNodeId = node->getId();

  int directionForFather;
  std::vector <Node*> nodes;
  nodes.push_back ( node->getFather() );
  if ( sonNumber==0 ) { //If sonNumber==0, the subtree we're interested in is composed of son 1 and father of node.
    nodes.push_back ( node->getSon ( 1 ) );
  }
  else { //If sonNumber==1, the subtree we're interested in is composed of son 0 and father of node.
    nodes.push_back ( node->getSon ( 0 ) );
  }
  if ( node->getFather()->getSon ( 0 ) ==node ) {
    directionForFather = 1; //node #1 is son 0, except at the root
  }
  else {
    directionForFather = 2; //node #2 is son 1, except at the root
  }

  int idNode0, idNode1;
  idNode0 = nodes[0]->getId();
  idNode1 = nodes[1]->getId();
  unsigned int directionNode0, directionNode1;
  directionNode0 = directionForFather;
  directionNode1 = 0;

  computeConditionalLikelihoodAndAssignSpId ( spTree, nodes,
      likelihoodData[geneNodeId][sonNumber+1],
      likelihoodData[idNode0][directionNode0],
      likelihoodData[idNode1][directionNode1],
      lossRates, duplicationRates,
      speciesIDs[geneNodeId][sonNumber+1],
      speciesIDs[idNode0][directionNode0],
      speciesIDs[idNode1][directionNode1],
      dupData[geneNodeId][sonNumber+1],
      dupData[idNode0][directionNode0],
      dupData[idNode1][directionNode1], false,
      cache);
  //Now we have the conditional likelihood of the upper subtree,
  //as well as the conditional likelihood of the lower subtree (which we already had)
  //We can thus compute the total likelihood of the rooting.

  std::vector <Node*> sons;
  sons.push_back ( node );

  sons.push_back ( node->getSon ( sonNumber ) );
  int idSon0 = geneNodeId;
  int idSon1 = sons[1]->getId();
  unsigned int directionSon0, directionSon1;
  directionSon0 = sonNumber+1;
  directionSon1 = 0;

  double rootLikelihood = 0.0;
  int rootSpId;
  int rootDupData = 0;

  computeConditionalLikelihoodAndAssignSpId ( spTree, sons,
      rootLikelihood,
      likelihoodData[idSon0][directionSon0],
      likelihoodData[idSon1][directionSon1],
      lossRates, duplicationRates,
      rootSpId, speciesIDs[idSon0][directionSon0],
      speciesIDs[idSon1][directionSon1],
      rootDupData, dupData[idSon0][directionSon0],
      dupData[idSon1][directionSon1], true,
      cache);
  while ( LksToNodes.find ( rootLikelihood ) !=LksToNodes.end() ) {
    // std::cout <<"changing rootLikelihood !!!!!!!!!!!!!!!!!!!"<<std::endl;
    rootLikelihood+=SMALLPROBA;
  }
  LksToNodes[rootLikelihood]=node->getSon ( sonNumber );
}



/*****************************************************************************
 * This function performs a preorder tree traversal in order to find likelihoods for rootings.
 * When used after the postorder tree traversal function, likelihoods for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 *
 ****************************************************************************/
void FastReconciliationTools::computeSubtreeLikelihoodPreorder ( TreeTemplate<Node> & spTree,
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
    std::map <double, Node*> & LksToNodes )

{
  ReconciliationCache cache(TreeTools::getMaxId(spTree, spTree.getRootId()));
  computeSubtreeLikelihoodPreorderIter(spTree, geneTree,node,seqSp,spID,likelihoodData,lossRates,
      duplicationRates,speciesIDs,dupData,sonNumber,LksToNodes, cache);
}

void FastReconciliationTools::computeSubtreeLikelihoodPreorderIter ( TreeTemplate<Node> & spTree,
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
    std::map <double, Node*> & LksToNodes ,
    ReconciliationCache &cache)
{
  computeRootingLikelihood ( spTree, node,
      likelihoodData, lossRates,
      duplicationRates, speciesIDs,
      dupData, sonNumber, LksToNodes, cache);
  if ( node->isLeaf() ) {
    return;
  }
  Node * son;
  if ( sonNumber==1 ) {
    son= node->getSon ( 1 );
  }
  else {
    son= node->getSon ( 0 );
  }
  {
    for ( unsigned int j =0; j<son->getNumberOfSons(); j++ ) {
      computeSubtreeLikelihoodPreorderIter ( spTree, geneTree,
          son, seqSp, spID,
          likelihoodData,
          lossRates, duplicationRates,
          speciesIDs, dupData, j, LksToNodes, cache);
    }
  }
}

void FastReconciliationTools::resetVector ( std::vector<unsigned int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0;
  }
}

void FastReconciliationTools::resetVector ( std::vector<int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0;
  }
}

void FastReconciliationTools::resetVector ( std::vector<double> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0.0;
  }
}

/**************************************************************************
 * Computes the probability of a given number of lineages numberOfLineages at the end of a branch given that there was only one at the beginning of the branch.
 **************************************************************************/
double FastReconciliationTools::computeBranchProbability ( const double & duplicationProbability, const double & lossProbability, const int numberOfLineages ) {
  double elm = exp ( duplicationProbability - lossProbability );
  double A = elm -1;
  double B = duplicationProbability * elm - lossProbability;
  double C = lossProbability * elm - duplicationProbability;
  double res =0;
  if ( numberOfLineages == 0 ) {
    res = lossProbability * A / B;
  }
  else {
    res = pow ( duplicationProbability, numberOfLineages ) * lossProbability * pow ( A,1+numberOfLineages ) * pow ( B, - ( 1+numberOfLineages ) ) - pow ( duplicationProbability, numberOfLineages - 1 ) *  pow ( A,numberOfLineages - 1 ) *  pow ( B, -numberOfLineages ) * C ;
  }
  if ( ( res == 0 ) || std::isnan ( res ) ) {
    res = NumConstants::VERY_TINY() ;
  }
  return ( res );
}


/**************************************************************************
 * Computes the log of the probability.
 **************************************************************************/
double FastReconciliationTools::computeLogBranchProbability ( const double & duplicationProbability, const double & lossProbability, const int numberOfLineages ) {
  return ( log ( computeBranchProbability ( duplicationProbability, lossProbability, numberOfLineages ) ) );
}



void FastReconciliationTools::computeNumbersOfLineagesFromRoot ( TreeTemplate<Node> * spTree,
    TreeTemplate<Node> * geneTree,
    Node * node,
    const std::map<std::string, std::string > &seqSp,
    const std::map<std::string, int > &spID,
    std::vector <int> &num0lineages,
    std::vector <int> &num1lineages,
    std::vector <int> &num2lineages,
    std::vector <std::vector<int> > & speciesIDs,
    std::vector <std::vector<int> > & dupData,
    std::set <int> & branchesWithDuplications )
{
  computeNumbersOfLineagesFromRootIter(spTree,geneTree,node,seqSp,spID,num0lineages,
      num1lineages,num2lineages,speciesIDs,dupData,branchesWithDuplications);
}

void FastReconciliationTools::computeNumbersOfLineagesFromRootIter ( TreeTemplate<Node> * spTree,
    TreeTemplate<Node> * geneTree,
    Node * node,
    const std::map<std::string, std::string > &seqSp,
    const std::map<std::string, int > &spID,
    std::vector <int> &num0lineages,
    std::vector <int> &num1lineages,
    std::vector <int> &num2lineages,
    std::vector <std::vector<int> > & speciesIDs,
    std::vector <std::vector<int> > & dupData,
    std::set <int> & branchesWithDuplications )
{
  int id=node->getId();
  if ( node->isLeaf() ) {
    speciesIDs[id][0]=speciesIDs[id][1]=speciesIDs[id][2]=assignSpeciesIdToLeaf ( node, seqSp, spID );
    num1lineages[speciesIDs[id][0]]+=1;
    dupData[id][0] = dupData[id][1] = dupData[id][2] = 1;
    return;
  }
  else {
    std::vector <Node *> sons = node->getSons();
    {
      for ( unsigned int i = 0; i< sons.size(); i++ ) {
        computeNumbersOfLineagesFromRootIter ( spTree, geneTree, sons[i],
            seqSp, spID, num0lineages,
            num1lineages, num2lineages,
            speciesIDs, dupData,
            branchesWithDuplications );
      }
    }
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0 = 0;
    unsigned int directionSon1 = 0;

    computeNumbersOfLineagesInASubtree ( *spTree, sons,
        speciesIDs[id][0],
        speciesIDs[idSon0][directionSon0],
        speciesIDs[idSon1][directionSon1],
        dupData[id][0], dupData[idSon0][directionSon0],
        dupData[idSon1][directionSon1],
        TreeTemplateTools::isRoot ( *node ),
        num0lineages, num1lineages,
        num2lineages, branchesWithDuplications );
  }
}

/*****************************************************************************
 * This function returns the speciesID assigned to a leaf.
 *
 ****************************************************************************/
int FastReconciliationTools::assignSpeciesIdToLeaf ( Node * node,  
    const std::map<std::string, std::string > & seqSp,
    const std::map<std::string, int > & spID )
{
  std::map<std::string, std::string >::const_iterator seqtosp;
  seqtosp=seqSp.find ( node->getName() );
  if ( seqtosp!=seqSp.end() ) {
    std::map<std::string, int >::const_iterator sptoid;
    sptoid = spID.find ( seqtosp->second );
    if ( sptoid!=spID.end() ) {
      return ( sptoid->second );
    }
    else {
      std::cerr <<"Error in assignSpeciesIdToLeaf: "<< seqtosp->second <<" not found in std::map spID for sequence "<< seqtosp->first<<std::endl;
      exit ( -1 );
    }
  }
  else {
    std::cerr <<"Error in assignSpeciesIdToLeaf: "<< node->getName() <<" not found in std::map seqSp"<<std::endl;
    exit ( -1 );
  }
}

/*****************************************************************************
 * This function recovers gene losses by comparing a subtree in a gene tree to
 * a species tree (tree).
 *
 ****************************************************************************/
void FastReconciliationTools::recoverLosses(Node *& node, int & a, const int & b, int & olda, const int & a0,
    const TreeTemplate<Node> & tree,
    double & likelihoodCell,
    const std::vector< double> & lossRates,
    const std::vector< double> & duplicationRates,
    ReconciliationCache &cache)
{
  olda = a;
  Node* nodeA;
  if (node->hasFather()) {
    nodeA = node->getFather();
  }
  else {
    std::cout <<"Problem in recoverLosses, nodeA has no father"<<std::endl;
  }
  a = nodeA->getId();
  node = nodeA;
  int lostNodeId = -1;
  if (   (nodeA->getSon(0)->getId()==olda)
      && (!(cache.contains(nodeA->getSon(1), b) || nodeA->getSon(1)->getId() == b))
      && (b!=a))

  {
    lostNodeId=nodeA->getSon(1)->getId();
    likelihoodCell += cache.computeLogBranchProbabilityCached(lostNodeId, 0, duplicationRates, lossRates);
  }
  else  if ((nodeA->getSon(1)->getId()==olda)
      && (!(cache.contains(nodeA->getSon(0), b) || (nodeA->getSon(0)->getId() == b)))
      && (b!=a))
  {
    lostNodeId=nodeA->getSon(0)->getId();
    likelihoodCell += cache.computeLogBranchProbabilityCached(lostNodeId, 0, duplicationRates, lossRates);
  }

  if ((olda!=a0)) {
    likelihoodCell += cache.computeLogBranchProbabilityCached(olda, 1, duplicationRates, lossRates);
  }
}

/*****************************************************************************
 * This function recovers gene losses by comparing a subtree in a gene tree to
 * a species tree (tree), when a duplication has affected the subtree.
 *
 ****************************************************************************/

void FastReconciliationTools::recoverLossesWithDuplication ( const Node * nodeA,
    const int &a,
    const int &olda,
    const TreeTemplate<Node> & tree,
    double & likelihoodCell,
    const std::vector< double> & lossRates,
    const std::vector< double> & duplicationRates )
{
  //The loss has occured before a0
  const Node * nodeOldA ;
  const Node * lostNode;
  if ( nodeA->getSon ( 0 )->getId() == olda ) {
    nodeOldA = nodeA->getSon ( 0 );
    lostNode=nodeOldA->getFather()->getSon ( 1 );
  }
  else {
    nodeOldA = nodeA->getSon ( 1 );
    lostNode=nodeOldA->getFather()->getSon ( 0 );
  }
  //We need to place the loss event in the right lineage
  likelihoodCell += computeLogBranchProbability ( duplicationRates[lostNode->getId()], lossRates[lostNode->getId()], 0 );
}



void FastReconciliationTools::computeNumbersOfLineagesInASubtree ( TreeTemplate<Node> & tree,
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
    std::set <int> &branchesWithDuplications )
{
  int a, a0, olda;
  int b, b0, oldb;
  a = a0 = olda = son0SpId;
  b = b0 = oldb = son1SpId;

  Node * temp0 = tree.getNode ( son0SpId );
  Node * temp1 = tree.getNode ( son1SpId );


  while ( a!=b ) { //There have been losses !
    if ( a>b ) {
      recoverLossesAndLineages ( temp0, a, b, olda, a0, tree, son0DupData, num0lineages, num1lineages );
    }
    else {
      recoverLossesAndLineages ( temp1, b, a, oldb, b0, tree, son1DupData, num0lineages, num1lineages );
    }
  }
  rootSpId = a;
  if ( ( a==a0 ) || ( b==b0 ) ) { //There has been a duplication !
    // std::cout <<"INSERTING"<<std::endl;
    branchesWithDuplications.insert ( sons[0]->getFather()->getId() );
    branchesWithDuplications.insert ( sons[0]->getId() ); //We also include the two son nodes themselves
    branchesWithDuplications.insert ( sons[1]->getId() );
    branchesWithDuplications.insert ( sons[0]->getId() );
    if ( sons[1]->getNumberOfSons () > 0 ) { //and the grandsons, if present.
      branchesWithDuplications.insert ( sons[1]->getSon ( 0 )->getId() );
      branchesWithDuplications.insert ( sons[1]->getSon ( 1 )->getId() );
    }
    if ( sons[0]->getNumberOfSons () > 0 ) {
      branchesWithDuplications.insert ( sons[0]->getSon ( 0 )->getId() );
      branchesWithDuplications.insert ( sons[0]->getSon ( 1 )->getId() );
    }

    if ( ( a==a0 ) && ( b==b0 ) ) {
      rootDupData += son0DupData+son1DupData;
      /* rootLikelihood-=(computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData) +
       *             computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData));*/
      if ( son0DupData==1 ) {
        num1lineages[a0]=num1lineages[a0]-1;
      }
      else if ( son0DupData>=2 ) { //All branches with 2 or more lineages
        num2lineages[a0]=num2lineages[a0]-1;
      }
      if ( son1DupData==1 ) {
        num1lineages[b0]=num1lineages[b0]-1;
      }
      else if ( son1DupData>=2 ) { //All branches with 2 or more lineages
        num2lineages[b0]=num2lineages[b0]-1;
      }

    }//there has been no loss, here
    else if ( b==b0 ) { //The loss has occured before a0
      rootDupData += son1DupData+1;
      //rootLikelihood-=computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
      if ( son1DupData==1 ) {
        num1lineages[b0]=num1lineages[b0]-1;
      }
      else if ( son1DupData>=2 ) { //All branches with 2 or more lineages
        num2lineages[b0]=num2lineages[b0]-1;
      }
      recoverLossesAndLineagesWithDuplication ( temp0, a, olda, tree, num0lineages );
    }
    else { //The loss has occured before b0
      rootDupData += son0DupData+1;
      //rootLikelihood-=computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
      if ( son0DupData==1 ) {
        num1lineages[a0]=num1lineages[a0]-1;
      }
      else if ( son0DupData>=2 ) { //All branches with 2 or more lineages
        num2lineages[a0]=num2lineages[a0]-1;
      }
      recoverLossesAndLineagesWithDuplication ( temp1, b, oldb, tree, num0lineages );
    }
    //Counting the duplication(s) event(s)
    if ( rootDupData==1 ) {
      num1lineages[rootSpId]+=1;
    }
    else if ( rootDupData>=2 ) { //All branches with 2 or more lineages
      num2lineages[rootSpId]+=1;
    }
    //  }
  }
  else { //there was no duplication
    rootDupData = 1;
    num1lineages[rootSpId]+=1;
  }
  return;
}

/*****************************************************************************
 * This set of functions aims at filling the num*lineages std::vectors.
 * It performs a post-order tree traversal using the root previously found by
 * double-recursive tree traversal.
 * Meanwhile, it also lists nodes where there has been a duplication in
 * std::vector<int> branchesWithDuplications.
 * This std::vector can be useful for making NNIs only around nodes showing duplications.
 ****************************************************************************/
void FastReconciliationTools::recoverLossesAndLineages ( Node *& node, int & a, const int & b, int & olda, int & a0,
    const TreeTemplate<Node> & tree,
    int & dupData, std::vector<int> &num0lineages, std::vector<int> &num1lineages )
{
  olda=a;
  Node* nodeA;
  if ( node->hasFather() ) {
    nodeA = node->getFather();
  }
  else {
    std::cout <<"Problem in recoverLossesAndLineages, nodeA has no father"<<std::endl;
  }
  a = nodeA->getId();
  node = nodeA;
  std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId ( * ( nodeA->getSon ( 0 ) ) );
  nodesIds0.push_back ( nodeA->getSon ( 0 )->getId() );
  std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId ( * ( nodeA->getSon ( 1 ) ) );
  nodesIds1.push_back ( nodeA->getSon ( 1 )->getId() );
  int lostNodeId = -1;

  if ( ( nodeA->getSon ( 0 )->getId() ==olda ) && ( ! ( VectorTools::contains ( nodesIds1, b ) ) ) && ( b!=a ) ) {

    lostNodeId=nodeA->getSon ( 1 )->getId();
  }
  else  if ( ( nodeA->getSon ( 1 )->getId() ==olda ) && ( ! ( VectorTools::contains ( nodesIds0, b ) ) ) && ( b!=a ) ) {

    lostNodeId=nodeA->getSon ( 0 )->getId();
  }
  if ( lostNodeId!= -1 ) {
    //likelihoodCell += computeLogBranchProbability(duplicationRates[lostNodeId], lossRates[lostNodeId], 0);
    num0lineages[lostNodeId]+=1;
    //  std::cout << "A loss on branch "<<lostNodeId<<std::endl;
  }
  /*  if ((dupData>0)&&(olda==a0)) {
   *     likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], dupData);
   *     std::cout <<dupData<<" B genes on branch "<<olda<<std::endl;
   *     // dupData = 0; //Resetting the dupdata value
   }
   else {*/
  if ( ( olda!=a0 ) ) {
    //likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
    num1lineages[olda]+=1;
    // std::cout <<"1 I genes on branch "<<olda<<std::endl;
  }
  // std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;

  return;
}




/****************************************************************************/





void FastReconciliationTools::recoverLossesAndLineagesWithDuplication ( const Node * nodeA,
    const int &a,
    const int &olda,
    const TreeTemplate<Node> & tree,
    std::vector <int> &num0lineages )
{
  //The loss has occured before a0
  //  const Node * nodeA = tree.getNode(a);
  const Node * nodeOldA ;
  const Node * lostNode;
  if ( nodeA->getSon ( 0 )->getId() == olda ) {
    nodeOldA = nodeA->getSon ( 0 );
    lostNode=nodeOldA->getFather()->getSon ( 1 );
  }
  else {
    nodeOldA = nodeA->getSon ( 1 );
    lostNode=nodeOldA->getFather()->getSon ( 0 );
  }
  //We need to place the loss event in the right lineage
  //likelihoodCell += computeLogBranchProbability(duplicationRates[lostNode->getId()], lossRates[lostNode->getId()], 0);
  num0lineages[lostNode->getId()]+=1;
  return;
}





/****************************************************************************/


/****************************************************************************/



