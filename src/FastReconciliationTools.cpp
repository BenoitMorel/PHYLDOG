#include "FastReconciliationTools.h"    
#include "Constants.h"

using namespace std;

void computeMaxId(Node *node, int &maxId) 
{
  maxId = std::max(maxId, node->getId());
  for (unsigned int i = 0; i < node->getNumberOfSons(); ++i) {
    computeMaxId(node->getSon(i), maxId);
  }
}

void fillPreorderRec(Node *node, int &currentId, 
    std::vector<int> &preorderIds, 
    std::vector<int> &lastSons,
    unsigned int tab = 0) {
  preorderIds[node->getId()] = currentId++;
  if (node->isLeaf())
    return;
  for (unsigned int i = 0; i < node->getNumberOfSons(); ++i) {
    fillPreorderRec(node->getSon(i), currentId, preorderIds, lastSons, tab + 1);
  }
  lastSons[node->getId()] = currentId;
}

bool FastReconciliationTools::isDescendant(Node *father, int descendantId)
{
  if (father->getId() == descendantId)
    return true;
  if (father->isLeaf()) 
    return false;
  if (!father->hasFather()) 
    return true;
  unsigned int preorderFather = _speciesIdsPreorder[father->getId()]; 
  unsigned int preorderDescendant = _speciesIdsPreorder[descendantId]; 
  unsigned int lastSon = _speciesIdsLastSon[father->getId()];
  return preorderFather < preorderDescendant && preorderDescendant < lastSon;
}

void FastReconciliationTools::initialize()
{
  _maxSpeciesId = 1;
  computeMaxId(_speciesTree.getRootNode(), _maxSpeciesId);
  _speciesNodes = std::vector<Node *>(_maxSpeciesId + 1, 0);
  std::vector<Node *> nodes = _speciesTree.getNodes();
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    _speciesNodes[nodes[i]->getId()] = nodes[i];
  }
  _speciesIdsPreorder = std::vector<int>(_maxSpeciesId, 0);
  _speciesIdsLastSon = std::vector<int>(_maxSpeciesId, 0);
  int currentId = 1;
  fillPreorderRec(_speciesTree.getRootNode(), currentId, _speciesIdsPreorder, _speciesIdsLastSon);
  _logBranchProbabilities.resize(3);
  _logBranchProbabilities[0] = std::vector<double>(_maxSpeciesId + 1, 0.0);
  _logBranchProbabilities[1] = std::vector<double>(_maxSpeciesId + 1, 0.0);
  _logBranchProbabilities[2] = std::vector<double>(_maxSpeciesId + 1, 0.0);
}


FastReconciliationTools::FastReconciliationTools(TreeTemplate<Node> * spTree,
        TreeTemplate<Node> * geneTree,
        const std::map<std::string, std::string > seqSp,
        const std::map<std::string, int > spID,
        std::vector< double> lossRates,
        std::vector < double> duplicationRates,
        std::vector <int> &num0lineages,
        std::vector <int> &num1lineages,
        std::vector <int> &num2lineages,
        std::set <int> &nodesToTryInNNISearch,
        bool fillTables):
  _speciesTree(*spTree),
  _geneTree(*geneTree),
  _seqSp(seqSp),
  _spID(spID),
  _lossRates(lossRates),
  _duplicationRates(duplicationRates),
  _fillTables(fillTables),
  _num0lineages(num0lineages),
  _num1lineages(num1lineages),
  _num2lineages(num2lineages),
  _nodesToTryInNNISearch(nodesToTryInNNISearch)
{
  assert(speciesTree->isRooted());
  assert(geneTreeTree->isRooted());
 
  initialize();
  std::vector <double> nodeData(3, 0.0);
  _likelihoodData = std::vector <std::vector<double> >(_geneTree.getNumberOfNodes(), nodeData);
  std::vector <int> nodeSpId(3, 0);
  _speciesIDs = std::vector <std::vector<int> >(_geneTree.getNumberOfNodes(), nodeSpId);
  _dupData = _speciesIDs;
}
    
FastReconciliationTools::~FastReconciliationTools()
{
}

double FastReconciliationTools::findMLReconciliationDR(int &MLindex) {
  Node * geneRoot = _geneTree.getRootNode();
  double initialLikelihood = computeSubtreeLikelihoodPostorder (geneRoot);

  std::vector <Node *> sons = geneRoot->getSons();

  if ( sons.size() !=2 ) {
    std::cerr <<"Error: "<<sons.size() << "sons at the root!"<<std::endl;
  }
  _LksToNodes[initialLikelihood] = sons[0];
  //We fill the likelihood and species ID data for the root node.
  //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
  //with other nodes.
  _likelihoodData[geneRoot->getId()][1] = _likelihoodData[geneRoot->getSon ( 1 )->getId()][0];
  _likelihoodData[geneRoot->getId()][2] = _likelihoodData[geneRoot->getSon ( 0 )->getId()][0];
  _speciesIDs[geneRoot->getId()][1] = _speciesIDs[geneRoot->getSon ( 1 )->getId()][0];
  _speciesIDs[geneRoot->getId()][2] = _speciesIDs[geneRoot->getSon ( 0 )->getId()][0];
  _dupData[geneRoot->getId()][1] = _dupData[geneRoot->getSon ( 1 )->getId()][0];
  _dupData[geneRoot->getId()][2] = _dupData[geneRoot->getSon ( 0 )->getId()][0];

  for ( unsigned int i = 0; i< sons.size(); i++ ) {
    for ( unsigned int j =0; j<sons[i]->getNumberOfSons(); j++ ) {
      computeSubtreeLikelihoodPreorder ( sons[i], j);
    }
  }
  vector<Node*> nodes = _geneTree.getNodes();
  for ( unsigned int i = 0 ; i < nodes.size() ; i++ ) {
    if ( nodes[i]->hasNodeProperty ( "outgroupNode" ) ) {
      nodes[i]->deleteNodeProperty ( "outgroupNode" );
      break;
    }
  }

  _LksToNodes.rbegin()->second->setNodeProperty ( "outgroupNode", BppString ( "here" ) );

  if ( _fillTables ) {

    //Now the best root has been found. I can thus run a function with this best root to fill all the needed tables. This additional tree traversal could be avoided.
    //To this end, the needed tables should be filled by the postfix and prefix traversals. This has not been done yet.
    //resetting
    std::vector <int> nodeSpId(3, 0);
    _speciesIDs = std::vector<std::vector<int> > ( _geneTree.getNumberOfNodes(), nodeSpId );
    _dupData = _speciesIDs;
    // Getting a well-rooted tree
    TreeTemplate<Node > * tree = _geneTree.clone();
    tree->newOutGroup ( _LksToNodes.rbegin()->second->getId() );
    _nodesToTryInNNISearch.clear();
    //Resetting numLineages std::vectors
    resetVector ( _num0lineages );
    resetVector ( _num1lineages );
    resetVector ( _num2lineages );
    computeNumbersOfLineagesFromRoot ( &_speciesTree, tree,
        tree->getRootNode(),
        _seqSp, _spID,
        _num0lineages, _num1lineages,
        _num2lineages, _speciesIDs,
        _dupData, _nodesToTryInNNISearch );
    delete tree;
  }
  
  //We return the best likelihood
  MLindex = _LksToNodes.rbegin()->second->getId();
  return _LksToNodes.rbegin()->first;
}


/*****************************************************************************
 * This function computes the lower conditional likelihood of a subtree and
 * assigns its summit node a species ID. Notations are influenced by
 * Zmasek and Eddy algorithm (2001).
 *
 ****************************************************************************/
double FastReconciliationTools::computeConditionalLikelihoodAndAssignSpId ( 
    double & rootLikelihood,
    double son0Likelihood,
    double son1Likelihood,
    int & rootSpId,
    int son0SpId,
    int son1SpId,
    int & rootDupData,
    int son0DupData,
    int son1DupData,
    bool atRoot)
{
  if ( rootLikelihood == 0.0 ) {
    int key = std::max(son0SpId, son1SpId) * 4242 +
      std::min(son0SpId, son1SpId);
    Data data = _assignMap[key];
    if (!atRoot && data.rootLikelihood != 0.0) {
      rootLikelihood = data.rootLikelihood + son0Likelihood + son1Likelihood;
      rootDupData = data.rootDupData;
      rootSpId = data.rootSpId;
      return rootLikelihood;;
    }
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
    Node * temp0 = _speciesNodes[son0SpId];
    Node * temp1 = _speciesNodes[son1SpId];

    double oldLL = rootLikelihood;
    while ( a!=b ) { //There have been losses !
      if ( a>b ) {
        recoverLosses(temp0, a, b, olda, a0, rootLikelihood);
      }
      else {
        recoverLosses(temp1, b, a, oldb, b0, rootLikelihood);
      }
    }
    rootSpId = a;
    if ( ( a==a0 ) || ( b==b0 ) ) { //There has been a duplication !
      if ( ( a==a0 ) && ( b==b0 ) ) {
        rootDupData += son0DupData+son1DupData;
        rootLikelihood-= ( computeLogBranchProbability (a0, son0DupData ) +
            computeLogBranchProbability (b0, son1DupData ) );
      }//there has been no loss, here
      else if ( b==b0 ) { //The loss has occured before a0
        rootDupData += son1DupData+1;
        rootLikelihood-=computeLogBranchProbability (b0, son1DupData );
        recoverLossesWithDuplication ( temp0, a, olda, _speciesTree, rootLikelihood, _lossRates, _duplicationRates );
      }
      else { //The loss has occured before b0
        rootDupData += son0DupData+1;
        rootLikelihood-=computeLogBranchProbability (a0, son0DupData );
        recoverLossesWithDuplication ( temp1, b, oldb, _speciesTree, rootLikelihood, _lossRates, _duplicationRates );
      }
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbability (a, rootDupData );
      }
      else {
        rootLikelihood += computeLogBranchProbability (a, rootDupData );
      }
    }
    else { //there was no duplication
      rootDupData = 1;
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbability (a, rootDupData );
      }
      else {
        rootLikelihood += computeLogBranchProbability (a, rootDupData );
      }
    }
    //Setting the lower conditional likelihood for the node of interest.
    data.rootLikelihood = rootLikelihood;
    data.rootDupData = rootDupData;
    data.rootSpId = rootSpId;
    _assignMap[key] =  data;
    rootLikelihood = data.rootLikelihood + son0Likelihood + son1Likelihood;
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
double FastReconciliationTools::computeSubtreeLikelihoodPostorder (Node *node) 
{
  return computeSubtreeLikelihoodPostorderIter(node);
}

double FastReconciliationTools::computeSubtreeLikelihoodPostorderIter (Node *node)
{
  int id=node->getId();
  if ( node->isLeaf() ) {
    if ( _likelihoodData[id][0]==0.0 ) {
      _speciesIDs[id][0] = _speciesIDs[id][1] = _speciesIDs[id][2]= assignSpeciesIdToLeaf ( node, _seqSp, _spID );
      _likelihoodData[id][0]= _likelihoodData[id][1]= _likelihoodData[id][2] =
        computeLogBranchProbability ( _speciesIDs[id][0], 1 );
      _dupData[id][0] = _dupData[id][1] = _dupData[id][2] = 1;
    }
    return _likelihoodData[id][0];
  }
  else {
    std::vector <Node *> sons = node->getSons();
    for ( unsigned int i = 0; i< sons.size(); i++ ) {
      computeSubtreeLikelihoodPostorderIter (sons[i]);
    }
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0 = 0;
    unsigned int directionSon1 = 0;

    computeConditionalLikelihoodAndAssignSpId (
        _likelihoodData[id][0],
        _likelihoodData[idSon0][directionSon0],
        _likelihoodData[idSon1][directionSon1],
        _speciesIDs[id][0],
        _speciesIDs[idSon0][directionSon0],
        _speciesIDs[idSon1][directionSon1],
        _dupData[id][0],
        _dupData[idSon0][directionSon0],
        _dupData[idSon1][directionSon1],
        TreeTemplateTools::isRoot ( *node ));
    return _likelihoodData[id][0];
  }
}


/*****************************************************************************
 * This function computes the likelihood of a rooting.
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 *
 ****************************************************************************/
void FastReconciliationTools::computeRootingLikelihood(Node * node,
    int sonNumber)
{
  int geneNodeId = node->getId();

  int directionForFather;
  Node *node0 = node->getFather();
  Node *node1 = 0;
  if ( sonNumber==0 ) { //If sonNumber==0, the subtree we're interested in is composed of son 1 and father of node.
    node1 = node->getSon(1);
  }
  else { //If sonNumber==1, the subtree we're interested in is composed of son 0 and father of node.
    node1 = node->getSon(0);
  }
  if ( node->getFather()->getSon ( 0 ) == node ) {
    directionForFather = 1; //node #1 is son 0, except at the root
  }
  else {
    directionForFather = 2; //node #2 is son 1, except at the root
  }

  int idNode0, idNode1;
  idNode0 = node0->getId();
  idNode1 = node1->getId();
  unsigned int directionNode0, directionNode1;
  directionNode0 = directionForFather;
  directionNode1 = 0;

  computeConditionalLikelihoodAndAssignSpId(
      _likelihoodData[geneNodeId][sonNumber+1],
      _likelihoodData[idNode0][directionNode0],
      _likelihoodData[idNode1][directionNode1],
      _speciesIDs[geneNodeId][sonNumber+1],
      _speciesIDs[idNode0][directionNode0],
      _speciesIDs[idNode1][directionNode1],
      _dupData[geneNodeId][sonNumber+1],
      _dupData[idNode0][directionNode0],
      _dupData[idNode1][directionNode1], false);
  //Now we have the conditional likelihood of the upper subtree,
  //as well as the conditional likelihood of the lower subtree (which we already had)
  //We can thus compute the total likelihood of the rooting.

  node0 = node;
  node1 = node->getSon(sonNumber);
  int idSon0 = geneNodeId;
  int idSon1 = node1->getId();
  unsigned int directionSon0, directionSon1;
  directionSon0 = sonNumber+1;
  directionSon1 = 0;

  double rootLikelihood = 0.0;
  int rootSpId;
  int rootDupData = 0;

  computeConditionalLikelihoodAndAssignSpId (
      rootLikelihood,
      _likelihoodData[idSon0][directionSon0],
      _likelihoodData[idSon1][directionSon1],
      rootSpId, _speciesIDs[idSon0][directionSon0],
      _speciesIDs[idSon1][directionSon1],
      rootDupData, _dupData[idSon0][directionSon0],
      _dupData[idSon1][directionSon1], true);
  while ( _LksToNodes.find ( rootLikelihood ) != _LksToNodes.end() ) {
    rootLikelihood+=SMALLPROBA;
  }
  _LksToNodes[rootLikelihood] = node1;
}



/*****************************************************************************
 * This function performs a preorder tree traversal in order to find likelihoods for rootings.
 * When used after the postorder tree traversal function, likelihoods for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 *
 ****************************************************************************/
void FastReconciliationTools::computeSubtreeLikelihoodPreorder (
    Node * node,
    int sonNumber)

{
  computeSubtreeLikelihoodPreorderIter(node, sonNumber);
}

void FastReconciliationTools::computeSubtreeLikelihoodPreorderIter (Node * node,
    int sonNumber)
{
  computeRootingLikelihood(node, sonNumber);
  if ( node->isLeaf() ) {
    return;
  }
  Node * son;
  if ( sonNumber == 1 ) {
    son = node->getSon ( 1 );
  }
  else {
    son = node->getSon ( 0 );
  }
  for ( unsigned int j =0; j<son->getNumberOfSons(); j++ ) {
    computeSubtreeLikelihoodPreorderIter (son, j);
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
double FastReconciliationTools::computeBranchProbability (int branch, int numberOfLineages ) {
  double lossProbability = _lossRates[branch];
  double duplicationProbability = _duplicationRates[branch];
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
double FastReconciliationTools::computeLogBranchProbability (int branch, int numberOfLineages ) {
  if (numberOfLineages > 2) {
      return  ( log ( computeBranchProbability (branch, numberOfLineages)));
  }
  double res = _logBranchProbabilities[numberOfLineages][branch];
    if (res == 0.0) {
      res =  ( log ( computeBranchProbability (branch, numberOfLineages)));
      _logBranchProbabilities[numberOfLineages][branch] = res;
    }
    return res;
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
void FastReconciliationTools::recoverLosses(Node *& node, 
    int &a, 
    int b, 
    int &olda, 
    int a0,
    double &likelihoodCell)
{
  olda = a;
  node = node->getFather();
  a = node->getId();
  Node *son0 = node->getSon(0);
  Node *son1 = node->getSon(1);
  int son0Id = son0->getId();
  int son1Id = son1->getId();
  
  if (son1Id != olda && b != a && (!isDescendant(son1, b)))
    likelihoodCell += computeLogBranchProbability(son1Id, 0);
  else if (son0Id != olda && b != a && !isDescendant(son0, b))
    likelihoodCell += computeLogBranchProbability(son0Id, 0);

  if ((olda!=a0)) 
    likelihoodCell += computeLogBranchProbability(olda, 1);
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
  likelihoodCell += computeLogBranchProbability ( lostNode->getId(), 0 );
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
    num0lineages[lostNodeId]+=1;
  }
  if ( ( olda!=a0 ) ) {
    num1lineages[olda]+=1;
  }
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
  num0lineages[lostNode->getId()]+=1;
  return;
}





/****************************************************************************/


/****************************************************************************/



