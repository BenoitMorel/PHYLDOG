/* This file contains various functions useful for reconciliations, such as reconciliation computation, printing of trees with integer indexes, search of a root with reconciliation...*/

#include "ReconciliationTools.h"


/**************************************************************************
* Assign arbitrary branch lengths (0.1) to branches if they do not have branch lengths
 **************************************************************************/


void assignArbitraryBranchLengths(TreeTemplate<Node> & tree){
			std::vector <Node *> ns = tree.getNodes();
			for (std::vector <Node *>::iterator it = ns.begin() ; it!= ns.end(); it++) {
				if ((*it)->hasFather()) {
					if (!(*it)->hasDistanceToFather()) {
						(*it)->setDistanceToFather(1);
					}
				}
			}


}

/**************************************************************************
 * This function re-numbers nodes in a binary tree with a pre-order tree traversal, so that the root has got index 0 and the underlying nodes have higher indexes.
 **************************************************************************/

void reNumber (TreeTemplate<Node> & tree, Node * noeud, int & index) {
  noeud->setId(index);
  index=index+1;
  if (! noeud->isLeaf()) {
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);
    reNumber(tree, son0, index);
    reNumber(tree, son1, index);
  }
}



void reNumber (TreeTemplate<Node> & tree) {
  Node * root = tree.getRootNode();
  int index=0;
  root->setId(index);
  index=index+1;
  Node * son0=root->getSon(0);
  Node * son1=root->getSon(1);
  reNumber(tree, son0, index);
  reNumber(tree, son1, index);
}

/**************************************************************************
 * This function re-numbers nodes with a breadth-first traversal.
 * 0 =white, 1 = grey, 2 = black
 * It also updates std::vectors of duplication and loss probabilities.
 * It should also return a std::map giving the correspondence between depth levels and node Ids
 **************************************************************************/
std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <double> dupProba=duplicationProbabilities;
  std::vector <double> lossProba=lossProbabilities;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for (unsigned int i = 0; i< nodes.size() ; i++) {
    color.insert(std::pair <Node *,int>(nodes[i],0));
  }
  std::queue <Node *> toDo;
  toDo.push(tree.getRootNode());
  color[tree.getRootNode()] = 1;
  double dupValue=duplicationProbabilities[tree.getRootNode()->getId()];
  double lossValue=lossProbabilities[tree.getRootNode()->getId()];
  dupProba[index]=dupValue;
  lossProba[index]=lossValue;
  tree.getRootNode()->setId(index);
  std::vector <int> v;
  DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
  DepthToIds[0].push_back(index);
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while(!toDo.empty()) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
      sons.push_back(u->getSon(j));
    }
    for (unsigned int j = 0; j< sons.size() ; j++) {
      if (color[sons[j]]==0) {
	color[sons[j]]=1;
	dupValue=duplicationProbabilities[sons[j]->getId()];
	lossValue=lossProbabilities[sons[j]->getId()];
	dupProba[index]=dupValue;
	lossProba[index]=lossValue;
	sons[j]->setId(index);
	if (DepthToIds.count(fatherDepth+1)==0) {
	  DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
	}
	DepthToIds[fatherDepth+1].push_back(index); 
	IdsToDepths[index] = fatherDepth+1;
	index++;
	toDo.push(sons[j]);
      }
    }
    color[u]=2;
  }
  duplicationProbabilities=dupProba;
  lossProbabilities=lossProba;
  return DepthToIds;
}

/**************************************************************************
 * This function re-numbers nodes with a breadth-first traversal.
 * 0 =white, 1 = grey, 2 = black
 * It also updates std::vectors of duplication and loss probabilities and of coalescent branch lengths.
 * It should also return a std::map giving the correspondence between depth levels and node Ids
 **************************************************************************/
std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities, std::vector <double> & coalBl) {
    int index = 0;
    std::map<Node *, int> color ;
    std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
    std::map <int, int > IdsToDepths;
    std::vector <double> dupProba=duplicationProbabilities;
    std::vector <double> lossProba=lossProbabilities;
    std::vector <double> cBl=coalBl;
    std::vector <Node * > nodes = tree.getNodes();
    //All nodes white
    for (unsigned int i = 0; i< nodes.size() ; i++) {
        color.insert(std::pair <Node *,int>(nodes[i],0));
    }
    std::queue <Node *> toDo;
    toDo.push(tree.getRootNode());
    color[tree.getRootNode()] = 1;
    double dupValue=duplicationProbabilities[tree.getRootNode()->getId()];
    double lossValue=lossProbabilities[tree.getRootNode()->getId()];
    double cValue = coalBl[tree.getRootNode()->getId()];
    dupProba[index]=dupValue;
    lossProba[index]=lossValue;
    cBl[index] = cValue;
    tree.getRootNode()->setId(index);
    std::vector <int> v;
    DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
    DepthToIds[0].push_back(index);
    IdsToDepths[index] = 0;
    index++;
    Node * u;
    while(!toDo.empty()) {
        u = toDo.front();
        toDo.pop();
        int fatherDepth = IdsToDepths[u->getId()];
        std::vector <Node *> sons;
        for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
            sons.push_back(u->getSon(j));
        }
        for (unsigned int j = 0; j< sons.size() ; j++) {
            if (color[sons[j]]==0) {
                color[sons[j]]=1;
                dupValue=duplicationProbabilities[sons[j]->getId()];
                lossValue=lossProbabilities[sons[j]->getId()];
                cValue = coalBl[sons[j]->getId()];
                dupProba[index]=dupValue;
                lossProba[index]=lossValue;
                cBl[index] = cValue;
                sons[j]->setId(index);
                if (DepthToIds.count(fatherDepth+1)==0) {
                    DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
                }
                DepthToIds[fatherDepth+1].push_back(index); 
                IdsToDepths[index] = fatherDepth+1;
                index++;
                toDo.push(sons[j]);
            }
        }
        color[u]=2;
    }
    duplicationProbabilities=dupProba;
    lossProbabilities=lossProba;
    coalBl = cBl;
    return DepthToIds;
}


/**************************************************************************
 * This function re-numbers nodes with a breadth-first traversal.
 * 0 =white, 1 = grey, 2 = black
 * It also updates a std::vector of coalescent branch lengths.
 * It should also return a std::map giving the correspondence between depth levels and node Ids
 **************************************************************************/
std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree, std::vector <double> & coalBl) {
    int index = 0;
    std::map<Node *, int> color ;
    std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
    std::map <int, int > IdsToDepths;
    std::vector <double> cBl=coalBl;
    std::vector <Node * > nodes = tree.getNodes();
    //All nodes white
    for (unsigned int i = 0; i< nodes.size() ; i++) {
        color.insert(std::pair <Node *,int>(nodes[i],0));
    }
    std::queue <Node *> toDo;
    toDo.push(tree.getRootNode());
    color[tree.getRootNode()] = 1;
    double cValue = coalBl[tree.getRootNode()->getId()];
    cBl[index] = cValue;
    tree.getRootNode()->setId(index);
    std::vector <int> v;
    DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
    DepthToIds[0].push_back(index);
    IdsToDepths[index] = 0;
    index++;
    Node * u;
    while(!toDo.empty()) {
        u = toDo.front();
        toDo.pop();
        int fatherDepth = IdsToDepths[u->getId()];
        std::vector <Node *> sons;
        for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
            sons.push_back(u->getSon(j));
        }
        for (unsigned int j = 0; j< sons.size() ; j++) {
            if (color[sons[j]]==0) {
                color[sons[j]]=1;
                cValue = coalBl[sons[j]->getId()];
                cBl[index] = cValue;
                sons[j]->setId(index);
                if (DepthToIds.count(fatherDepth+1)==0) {
                    DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
                }
                DepthToIds[fatherDepth+1].push_back(index); 
                IdsToDepths[index] = fatherDepth+1;
                index++;
                toDo.push(sons[j]);
            }
        }
        color[u]=2;
    }
    coalBl = cBl;
    return DepthToIds;
}


//There we do not update probabilities

std::map <int, std::vector <int> > breadthFirstreNumber (TreeTemplate<Node> & tree) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for (unsigned int i = 0; i< nodes.size() ; i++) {
    color.insert(std::pair <Node *,int>(nodes[i],0));
  }
  std::queue <Node *> toDo;
  toDo.push(tree.getRootNode());
  color[tree.getRootNode()] = 1;

  tree.getRootNode()->setId(index);

  std::vector <int> v;
  DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
  DepthToIds[0].push_back(index);
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while(!toDo.empty()) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
      sons.push_back(u->getSon(j));
    }
    for (unsigned int j = 0; j< sons.size() ; j++) {
      if (color[sons[j]]==0) {
	color[sons[j]]=1;
	sons[j]->setId(index);
	if (DepthToIds.count(fatherDepth+1)==0) {
	  DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
	}
	DepthToIds[fatherDepth+1].push_back(index); 
	IdsToDepths[index] = fatherDepth+1;
	index++;
	toDo.push(sons[j]);
      }
    }
    color[u]=2;
  }
  return DepthToIds;
}

//There we do not update probabilities, and we reset "toComp" node properties
std::map <int, std::vector <int> > breadthFirstreNumberAndResetProperties (TreeTemplate<Node> & tree) {
    int index = 0;
    std::map<Node *, int> color ;
    std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
    std::map <int, int > IdsToDepths;
    std::vector <Node * > nodes = tree.getNodes();
    //All nodes white
    for (unsigned int i = 0; i< nodes.size() ; i++) {
        color.insert(std::pair <Node *,int>(nodes[i],0));
        if (nodes[i]->hasNodeProperty("toComp")) {
            nodes[i]->setNodeProperty("toComp", BppString( "y" ) );
        }
    }
    std::queue <Node *> toDo;
    toDo.push(tree.getRootNode());
    color[tree.getRootNode()] = 1;
    tree.getRootNode()->setId(index);
    std::vector <int> v;
    DepthToIds.insert(std::pair <int, std::vector<int> > (0,v));
    DepthToIds[0].push_back(index);
    IdsToDepths[index] = 0;
    index++;
    Node * u;
    while(!toDo.empty()) {
        u = toDo.front();
        toDo.pop();
        int fatherDepth = IdsToDepths[u->getId()];
        std::vector <Node *> sons;
        for (unsigned int j = 0 ; j< u->getNumberOfSons() ; j++) {
            sons.push_back(u->getSon(j));
        }
        for (unsigned int j = 0; j< sons.size() ; j++) {
            if (color[sons[j]]==0) {
                color[sons[j]]=1;
                sons[j]->setId(index);
                if (DepthToIds.count(fatherDepth+1)==0) {
                    DepthToIds.insert(std::pair <int, std::vector<int> > (fatherDepth+1,v));
                }
                DepthToIds[fatherDepth+1].push_back(index); 
                IdsToDepths[index] = fatherDepth+1;
                index++;
                toDo.push(sons[j]);
            }
        }
        color[u]=2;
    }
    return DepthToIds;
}




/**************************************************************************
 * These functions reset the LOSSES and DUPLICATIONS values on the species tree and in the associated std::vectors.
 **************************************************************************/

void changeNodeProperty(Node & noeud, const std::string & name, const Clonable & property) { 
  if (noeud.hasNodeProperty(name)) {
    noeud.deleteNodeProperty(name);
    noeud.setNodeProperty(name, property);
  }
  else {
    noeud.setNodeProperty(name, property);
  }
}

void changeBranchProperty(Node & noeud, const std::string & name, const Clonable & property) { 
  if (noeud.hasBranchProperty(name)) {
    noeud.deleteBranchProperty(name);
    noeud.setBranchProperty(name, property);
  }
  else {
    noeud.setBranchProperty(name, property);
  }
}



void resetLossesAndDuplications(TreeTemplate<Node> & tree, /*std::vector <int> &lossNumbers, */std::vector <double> &lossProbabilities, /*std::vector <int> &duplicationNumbers, */std::vector <double> &duplicationProbabilities) {
  std::vector< int > nodesIds = tree.getNodesId ();
//  Number<int> zero = Number<int>(0);
  BppString zero = BppString(TextTools::toString(0));
  for (unsigned int i=0; i<nodesIds.size(); i++) {
    changeBranchProperty(*(tree.getNode(nodesIds[i])),LOSSES, zero);
    changeBranchProperty(*(tree.getNode(nodesIds[i])),DUPLICATIONS, zero);
    /*lossNumbers[nodesIds[i]] = 0;
    duplicationNumbers[nodesIds[i]] = 0;*/
  }
}

void resetLossesDuplicationsSpeciations(TreeTemplate<Node> & tree, std::vector <int> &lossNumbers, std::vector <double> &lossProbabilities, std::vector <int> &duplicationNumbers, std::vector <double> &duplicationProbabilities, std::vector <int> &branchNumbers) {
  std::vector< int > nodesIds = tree.getNodesId ();
 // Number<int> zero = Number<int>(0);
  BppString zero = BppString(TextTools::toString(0));

  for (unsigned int i=0; i<nodesIds.size(); i++) {
    changeBranchProperty(*(tree.getNode(nodesIds[i])),LOSSES, zero);
    changeBranchProperty(*(tree.getNode(nodesIds[i])),DUPLICATIONS, zero);
    lossNumbers[nodesIds[i]] = 0;
    duplicationNumbers[nodesIds[i]] = 0;
    branchNumbers[nodesIds[i]] = 0;
  }
}



void resetLossesDuplicationsSpeciationsForGivenNodes(TreeTemplate<Node> & tree, std::vector <int> & lossNumbers, std::vector <double> & lossProbabilities, std::vector <int> & duplicationNumbers, std::vector <double> & duplicationProbabilities, std::vector <int> & branchNumbers, std::vector <int> nodesToUpdate, std::map <int, std::vector<int> > & geneNodeIdToLosses, std::map <int, int > & geneNodeIdToDuplications, std::map <int, std::vector<int> > & geneNodeIdToSpeciations){
    for (unsigned int i=0; i<nodesToUpdate.size(); i++) {
      unsigned int size = geneNodeIdToLosses[nodesToUpdate[i]].size();
      for (unsigned int j = 0 ; j< size ; j++) {
        lossNumbers[geneNodeIdToLosses[nodesToUpdate[i]][j]]--;
        changeBranchProperty(*(tree.getNode(geneNodeIdToLosses[nodesToUpdate[i]][j])),LOSSES, BppString(TextTools::toString((lossNumbers[geneNodeIdToLosses[nodesToUpdate[i]][j]]))));
      }
      size = geneNodeIdToSpeciations[nodesToUpdate[i]].size();
      for (unsigned int j = 0 ; j< size ; j++) {
	branchNumbers[geneNodeIdToSpeciations[nodesToUpdate[i]][j]]--;
      }
      if (geneNodeIdToDuplications[nodesToUpdate[i]] != -1) {
        duplicationNumbers[geneNodeIdToDuplications[nodesToUpdate[i]]]--;
        changeBranchProperty(*(tree.getNode(geneNodeIdToDuplications[nodesToUpdate[i]])), DUPLICATIONS, BppString(TextTools::toString((duplicationNumbers[geneNodeIdToDuplications[nodesToUpdate[i]]]))));
      }
    }


    //We do not reset the 3 std::maps, because it is currently not useful. It might be useful provided the program is changed in some way.
}

/**************************************************************************
 * This function sets all values in a std::vector<int> to 0.
 **************************************************************************/

void printVector(std::vector<int> & v) {
  unsigned int temp=v.size();
  for (unsigned int i=0;i<temp ; i++) {
    std::cout <<  "i : "<<i<< " vi : "<< v[i]<<std::endl;
  }
}

void printVectorLine(std::vector<int> & v) {
  unsigned int temp=v.size();
  for (unsigned int i=0;i<temp ; i++) {
    std::cout <<  "i : "<<i<< " vi : "<< v[i]<<" ";
  }
  std::cout <<std::endl;
}

void resetVector(std::vector<unsigned int> & v) {
    unsigned int temp=v.size();
    for (unsigned int i=0;i<temp ; i++) {
        v[i]=0;
    }
}

void resetVector(std::vector<int> & v) {
  unsigned int temp=v.size();
  for (unsigned int i=0;i<temp ; i++) {
    v[i]=0;
  }
}

void resetVector(std::vector<double> & v) {
  unsigned int temp=v.size();
  for (unsigned int i=0;i<temp ; i++) {
    v[i]=0.0;
  }
}


void resetVectorForGivenNodes(std::vector<int> & v, std::vector <int> nodesToUpdate) {
  unsigned int temp=nodesToUpdate.size();
  for (unsigned int i=0;i<temp ; i++) {
    v[nodesToUpdate[i]]=0;
  }
}

/**************************************************************************
 * This function resets species Ids in a tree.
 **************************************************************************/
void resetSpeciesIds (TreeTemplate<Node> & tree) {
  std::vector <Node *> nodes = tree.getNodes();
  for (unsigned int i = 0; i< nodes.size() ; i++) {
    changeNodeProperty(*(nodes[i]), SPECIESID, Number<int>(-1));
  }
}

void resetSpeciesIdsForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds) {
  for (unsigned int i = 0; i< nodesToUpdate.size() ; i++) {
    removedNodeIds.push_back((dynamic_cast<const Number<int> *>(tree.getNode(nodesToUpdate[i])->getNodeProperty(SPECIESID))->getValue())); 
    changeNodeProperty(*(tree.getNode(nodesToUpdate[i])), SPECIESID, Number<int>(-1));
  }
}

 void resetSpeciesIdsForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate) {
  for (unsigned int i = 0; i< nodesToUpdate.size() ; i++) {
    changeNodeProperty(*(tree.getNode(nodesToUpdate[i])), SPECIESID, Number<int>(-1));
  }
}
 
void resetSpeciesIdsAndLiks (TreeTemplate<Node> & tree) {
  std::vector <Node *> nodes = tree.getNodes();
  for (unsigned int i = 0; i< nodes.size() ; i++) {
    changeNodeProperty(*(nodes[i]), SPECIESID, Number<int>(-1));
    changeNodeProperty(*(nodes[i]), LOWLIK, Number<double>(1.0));  
    changeBranchProperty(*(nodes[i]),EVENTSPROBA, Number<double>(1.0)); 
   }
}

void resetSpeciesIdsAndLiksForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds) {
  for (unsigned int i = 0; i< nodesToUpdate.size() ; i++) {
    removedNodeIds.push_back((dynamic_cast<const Number<int> *>(tree.getNode(nodesToUpdate[i])->getNodeProperty(SPECIESID))->getValue())); 
    changeNodeProperty(*(tree.getNode(nodesToUpdate[i])), SPECIESID, Number<int>(-1));
    changeNodeProperty( *(tree.getNode(nodesToUpdate[i])), LOWLIK, Number<double>(1.0));  
    changeBranchProperty(*(tree.getNode(nodesToUpdate[i])), EVENTSPROBA, Number<double>(1.0));
  }
}

 void resetSpeciesIdsAndLiksForGivenNodes (TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate) {
  for (unsigned int i = 0; i< nodesToUpdate.size() ; i++) {
    changeNodeProperty(*(tree.getNode(nodesToUpdate[i])), SPECIESID, Number<int>(-1));
    changeNodeProperty( *(tree.getNode(nodesToUpdate[i])), LOWLIK, Number<double>(1.0));  
    changeBranchProperty(*(tree.getNode(nodesToUpdate[i])), EVENTSPROBA, Number<double>(1.0));
    }
}
 

/**************************************************************************
 * This function uses a std::map giving the link between levels and nodes ids and sends the maximum id corresponding to a limit level.
 **************************************************************************/

int getLimitIdFromNumberOfLevels(std::map<int, std::vector <int> > DepthToIds, int levelsToRootAt){
  return(  VectorTools::max(DepthToIds[levelsToRootAt]));
}



/**************************************************************************
 * Computes the rates of duplication and loss of the birth-death process along a 
 * particular branch given counts of times where, at the end of the branch, 
 * there were 0 gene (i), 1 gene (j) or 2 genes (k). 
 * This is therefore an approximate formula (normally, 
 * one should consider all possible counts, from 0 to +infinity).
 **************************************************************************/
void computeDuplicationAndLossProbabilities (int i, int j, int k, 
                                             double & lossProbability, 
                                             double & duplicationProbability) {
  double id=double(i);
  double jd=double(j);
  double kd=double(k);

  //Linear Birth-Death process, from Bartholomay 1958
  double sum = id +jd +kd;
  double X = id/sum;
  double Y = jd/sum;
  double Z = kd/sum;
  //num2Lineages (Z) corresponds to the sum of cases where there were 2 or 
  //more than 2 lineages at the end of a branch.
  double ZoverYPlusZ = Z/(Y+Z);
  double ZoverXTimesYPlusZ = ZoverYPlusZ * 1/X;
  
  //expected number of losses over a branch of length 1.
  double mu = log( (X-1)/(ZoverYPlusZ-1) ) / (ZoverXTimesYPlusZ -1);
  
  //expected number of duplications over a branch of lenth 1.
  double lambda = ZoverXTimesYPlusZ * mu;
  
/* If we truncate 
  double Y2Z = pow(Y, 2)/Z;
  double XYZ = X*Y/Z;
  double oneXYZ = 1/XYZ;
  
  //expected number of losses over a branch of length 1.
  double mu = (log(XYZ - X + Y2Z) - log(1 - X + Y2Z) ) / (1- oneXYZ);
  
  //expected number of duplications over a branch of lenth 1.
  double lambda = oneXYZ * mu;
  */
  if (!(std::isnan(lambda)||std::isinf(lambda))) { 
    if (lambda<=0.000001) {
      duplicationProbability = 0.000001;//SMALLPROBA;
    }
    else {  
      duplicationProbability = lambda;
    }
  } 
  else if (std::isnan(lambda)) {
    duplicationProbability = 0.000001;
  }
  
  
  
  if (!(std::isnan(mu)||std::isinf(mu))) {
    if (mu<=0.0000011) {
      lossProbability = 0.0000011;//SMALLPROBA;
    }
    else {
      lossProbability = mu;
    }
  }  
  else if (std::isnan(mu)) {
    std::cerr<<"A branch loss probability could not be properly estimated. Set to 0.000001 by default. X:"<<X<<" Y:"<<Y<<" Z:"<<Z<<std::endl;
    lossProbability = 0.0000011;
  }
  
  
  
  
  
  /*Old Formulas as of 18 02 2010

  while((id < 1.0)||(jd < 1.0)||(kd < 1.0 )) {
    id = id*100;
    jd = jd*100;
    kd = kd*100;
  }
  
  
  if (id==kd) {
    kd=kd-0.00001;
  }
  
  double denom = kd*id-jd*kd-kd*kd+jd*id;
  if (denom==0) {
    denom=0.00000001;
  }
  double ln ;
  ln = log ((jd+2*kd)/(id+jd+kd));
  
   
  double temp = -ln*(jd+2*kd)*id/denom;
 
  if (!(std::isnan(temp)||std::isinf(temp))) {
    if (temp<=0) {
      lossProbability = 0.0001;//SMALLPROBA;
    }
     else {
      lossProbability = temp;
    }
  }
  temp = -ln*kd*(id+jd+kd)/denom;  
 
  if (!(std::isnan(temp)||std::isinf(temp))) { 
    if (temp<=0) {
      duplicationProbability = 0.0001;//SMALLPROBA;
    }
    else {
      duplicationProbability = temp;
    }
  }
  if (lossProbability <0.0001) {
    lossProbability = 0.0001;//1e-6;
  }
  if (duplicationProbability <0.0001) {
    duplicationProbability = 0.0001;//1e-6;
  }
  */
/*  
  //New formulas as of 18 02 2010
  double sum = id +jd +kd;
  double X = id/sum;
  double Y = jd/sum;
  double Z = kd/sum;
  
  double K = 1 - sqrt(1-4*Z/(1-X));
  
  double temp = (1/(1 - 2*X/K )) * log( ((X-1)/(Y)) * ((2*X/K)*(1-Y/(1-X))-1));
  
  if (!(std::isnan(temp)||std::isinf(temp))) { 
    if (temp<=0) {
      duplicationProbability = 0.0001;//SMALLPROBA;
    }
    else {  
      duplicationProbability = temp;
    }
  }
  temp = 2*X*temp/K;
  
  if (!(std::isnan(temp)||std::isinf(temp))) {
    if (temp<=0) {
      lossProbability = 0.0001;//SMALLPROBA;
    }
    else {
      lossProbability = temp;
          }
  }  
 
  if (lossProbability <0.0001) {
    lossProbability = 0.0001;//1e-6;
  }
  if (duplicationProbability <0.0001) {
    duplicationProbability = 0.0001;//1e-6;
  }
  */
    
  
  
 /* 
  std::cout<<"0: "<<i<<" 1: "<<j<<" 2: "<<k<<std::endl;
std::cout<<"X: "<<X<<" Y: "<<Y<<" Z: "<<Z<<std::endl;
  std::cout<<" lossProbability: "<<lossProbability <<" duplicationProbability: "<<duplicationProbability<<std::endl;
  */
 // std::cout <<"in computeDuplicationAndLossProbabilities 2"<<std::endl;

  return;
}

/*************************************************
 * For the root branch, a special trick is used. At the root, we never count families that have losses, because these families do not have any gene from the species involved in the study. So we set the number of losses at the average number observed in other branches. 
*************************************************/



void computeDuplicationAndLossProbabilitiesForAllBranches (std::vector <int> numOGenes, std::vector <int> num1Genes, std::vector <int> num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities) {
  //The trick for the root, already done in alterLineages...:
 /*
  int totNum0=0, totNum12=0;
  for (int i =1 ; i< lossProbabilities.size() ; i++) {
    totNum0+=numOGenes[i];
    totNum12+=num1Genes[i]+num2Genes[i];
  }
  //At branch 0, by definition, we never count cases where there has been a loss, so we set it to an average value.
  numOGenes[0] = (int)floor( ((double)totNum0/(double)totNum12)*((double)num1Genes[0]+(double)num2Genes[0]));
  */
  for (unsigned int i =0 ; i< lossProbabilities.size() ; i++) {
    computeDuplicationAndLossProbabilities (numOGenes[i], num1Genes[i], num2Genes[i], lossProbabilities[i], duplicationProbabilities[i]);
    }
  return;
}


void computeAverageDuplicationAndLossProbabilitiesForAllBranches (std::vector <int> numOGenes, std::vector <int> num1Genes, std::vector <int> num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities) {
  //The trick for the root:
 /* int totNum0=0, totNum12=0;
  for (int i =1 ; i< lossProbabilities.size() ; i++) {
    totNum0+=numOGenes[i];
    totNum12+=num1Genes[i]+num2Genes[i];
  }
  numOGenes[0] = (int)floor( ((double)totNum0/(double)totNum12)*((double)num1Genes[0]+(double)num2Genes[0]));
  */
  int sumOGene = VectorTools::sum(numOGenes);
  int sum1Gene = VectorTools::sum(num1Genes);
  int sum2Gene = VectorTools::sum(num2Genes);
  double dupProba, lossProba;
  computeDuplicationAndLossProbabilities(sumOGene, sum1Gene, sum2Gene, lossProba, dupProba);
  int size = lossProbabilities.size();
  if (!(std::isnan(lossProba)||std::isinf(lossProba)||(lossProba<0.0))) {
    lossProbabilities.assign (size, lossProba);
  }
  if (!(std::isnan(dupProba)||std::isinf(dupProba)||(dupProba<0.0))) {
    duplicationProbabilities.assign (size, dupProba);
  }
  return;
}



/**************************************************************************
 * Computes the probability of a given number of lineages numberOfLineages at the end of a branch given that there was only one at the beginning of the branch. 
 **************************************************************************/
double computeBranchProbability (double duplicationProbability, double lossProbability, int numberOfLineages) {
  //Using a linear birth-death process (Bartholomay, 1958)
  int min ;
  if (numberOfLineages == 0) {
    min = 0; 
  }
  else {
    min = 1;
  }
  
  double elm = exp(duplicationProbability - lossProbability);
  double A = elm -1;
  double B = duplicationProbability * elm - lossProbability;
  double C = lossProbability * elm - duplicationProbability;
  
  double res =0;
  
  if (numberOfLineages == 0) {
    res = lossProbability * A / B;
  }
  else {
    res = pow(duplicationProbability, numberOfLineages) * lossProbability * pow(A,1+numberOfLineages) * pow(B, -(1+numberOfLineages)) - pow(duplicationProbability, numberOfLineages - 1) *  pow(A,numberOfLineages - 1) *  pow(B, -numberOfLineages) * C ;
  }
  if ( (res == 0) || std::isnan(res) )
    {
     res = NumConstants::VERY_TINY ; 
    }
  return(res);

  
 /* Old computations as of 21 02 2010
  if (duplicationProbability==lossProbability) {
    lossProbability+=0.0000001;
  }
  double beta = (1-exp(duplicationProbability-lossProbability))/(lossProbability-duplicationProbability*exp(duplicationProbability-lossProbability));
  if (std::isnan(beta)) {
    std::cout <<"ISNAN beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
  }
  if (numberOfLineages==0) {
    return (lossProbability*beta);
  }
  else {
    double temp = (1-duplicationProbability*beta)*pow((duplicationProbability*beta),(numberOfLineages-1)) * (1-lossProbability*beta);
    if (std::isnan(temp)) {
      std::cout <<"ISNAN2 beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
    }
    else if (temp == 0) {
     // std::cerr <<"TEMP equals 0; beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
      temp = NumConstants::VERY_TINY ;
    }
    return(temp);
  }
  */
}


/**************************************************************************
 * Computes the log of the probability. 
 **************************************************************************/
double computeLogBranchProbability (double duplicationProbability, double lossProbability, int numberOfLineages) {
  //UNDONE TEST 1409
/*  if (numberOfLineages<=1) {
    return 0;
  }
  else { 
    return(-numberOfLineages+1);*/
   return(log(computeBranchProbability(duplicationProbability, lossProbability, numberOfLineages)));
// }
}


/**************************************************************************
 * Computes the probability of a given number of lineages numberOfLineages at the end of a branch given that there was 0 at the beginning of the branch. 
 **************************************************************************/
double computeBranchProbabilityAtRoot (double duplicationProbability, double lossProbability, int numberOfLineages) {
  if (duplicationProbability==lossProbability) {
    lossProbability+=0.0000001;
  }
  double beta = (1-exp(duplicationProbability-lossProbability))/(lossProbability-duplicationProbability*exp(duplicationProbability-lossProbability));
  if (std::isnan(beta)) {
    std::cout <<"ISNAN beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
  }
  if (numberOfLineages==0) {
    std::cerr <<"Error in computeBranchProbabilityAtRoot: cannot compute P(0 lineage)!"<<std::endl;
      MPI::COMM_WORLD.Abort(1);
    exit(-1);
//    return (lossProbability*beta);
  }
  else {
    double temp = (1-duplicationProbability*beta)*pow((duplicationProbability*beta),(numberOfLineages-1));
    if (std::isnan(temp)) {
      std::cout <<"ISNAN2 beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
    }
   else if (temp == 0) {
    // std::cerr <<"TEMP equals 0; beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
     temp = NumConstants::VERY_TINY ;
   }
    return(temp);
  }
}


/**************************************************************************
 * Computes the log of the above probability (if it had not been changed not to!):
 * Update as of January 28th, 2010.
 * Now we do not consider gene originations anymore.
 * Therefore events at the root are considered as on any other branch.
 **************************************************************************/
double computeLogBranchProbabilityAtRoot (double duplicationProbability, double lossProbability, int numberOfLineages) {
  //UNDONE TEST 1409
 /* if (numberOfLineages<=1) {
    return 0;
  }
  else {
   return(-numberOfLineages+1);
  }*/
  //Update as of January 28th, 2010.
  //Now we do not consider gene originations anymore.
  //Therefore, we do not make a special case when we are at the root.
// return(log(computeBranchProbabilityAtRoot(duplicationProbability, lossProbability, numberOfLineages)));
  //instead, we use:
  return(log(computeBranchProbability(duplicationProbability, lossProbability, numberOfLineages)));
}



/**************************************************************************
 * This function computes the score of a scenario in a gene tree given a species tree. 
 * It makes a post-order tree-traversal to annotate each node of the tree with the conditional likelihood of the underlying subtree, and each branch is annotated with a double giving the product of terms along this branch (i.e. product of duplication and loss probabilities along that branch). 
 * Getting the conditional likelihood at the root of the tree gives the scenario likelihood. 
 * The conditional likelihoods are useful because they can be reused when another root is tried. 
 * Moreover, the function also outputs counts of branches, duplications and losses.
 * Lower conditional likelihoods are noted LOWLIK.
 * Doubles associated with each branch and giving probabilities of all events along that branch are noted EVENTSPROBA.
 **************************************************************************/

void computeScenarioScore (TreeTemplate<Node> & tree, 
                           TreeTemplate<Node> & geneTree, 
                           Node * noeud, 
                           std::vector<int> &branchNumbers, 
                           std::map <int, std::vector <int> > &geneNodeIdToSpeciations, 
                           std::vector<double>duplicationProbabilities, 
                           std::vector <double> lossProbabilities, 
                           std::vector <int> &num0lineages, 
                           std::vector <int> &num1lineages, 
                           std::vector <int> &num2lineages) {
  int num = (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue()); 
  Node * spNode = tree.getNode(num);
  if (noeud->isLeaf()) {
    int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
    double probaBranch ;
    if (*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D") {
      //do nothing, we already counted the probability of the duplication event on one of the upper branches !
      probaBranch=0.0;
      changeNodeProperty(*(noeud), NUMGENES, Number<int>(1)); 
      if (num!=fatherNum) {
        //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
        probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], 1);
        num1lineages[num]++;
        branchNumbers[num]++;
      }
    }
    else {
      //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
      probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], 1);
      changeNodeProperty(*(noeud), NUMGENES, Number<int>(1));
      num1lineages[num]++;
      branchNumbers[num]++;
    }
    //if not at the root, there may be a long story of losses from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.
    
    if (spNode->hasFather()) {
      int tempNum = num;
      Node * spTempNode = spNode;
      //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
      int tempNumBrother;
      while ((tempNum!=fatherNum) && (spTempNode->hasFather()) && (spTempNode->getFather()->getId()!=fatherNum)) {
        //There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication
        if (spTempNode->getFather()->getSon(0)->getId()==tempNum) {
          tempNumBrother = spTempNode->getFather()->getSon(1)->getId();
        }
        else {
          tempNumBrother = spTempNode->getFather()->getSon(0)->getId();
        }
        
        probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0);
        num0lineages[tempNumBrother]++;
        branchNumbers[tempNumBrother]++; 
        spTempNode = spTempNode->getFather();
        tempNum = spTempNode->getId();
        //No duplication on branch leading to node now named tempNum, and no loss either
        probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNum], lossProbabilities[tempNum], 1);
        num1lineages[tempNum]++;
        branchNumbers[tempNum]++;
      }
      if ((*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D")&&(num!=fatherNum)) {
        if (spTempNode->getFather()->getSon(0)->getId()==tempNum) {
          tempNumBrother = spTempNode->getFather()->getSon(1)->getId();
        }
        else {
          tempNumBrother = spTempNode->getFather()->getSon(0)->getId();
        }
        probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0);	
        num0lineages[tempNumBrother]++;
        branchNumbers[tempNumBrother]++;
      }
    }
    
    changeBranchProperty(*(noeud), EVENTSPROBA, Number<double>(probaBranch));
    changeNodeProperty(*(noeud), LOWLIK, Number<double>(probaBranch));
  }
  
  else { //noeud is not a leaf
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);  
    if ((dynamic_cast<const Number<double> *>(son0->getNodeProperty(LOWLIK))->getValue())==1.0)
      {
        computeScenarioScore (tree, geneTree, son0, branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages);
      }
    if ((dynamic_cast<const Number<double> *>(son1->getNodeProperty(LOWLIK))->getValue())==1.0)
      {
        computeScenarioScore (tree, geneTree, son1, branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages);
      }
    int num = (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue());
    int num0 = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue());
    int num1 = (dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue());
    double probaBranch ;
    if(num==num0||num==num1) { //there is a duplication at this branch
      if (noeud->hasFather()) {
        int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
        int brotherNum;
        if (noeud->getFather()->getSon(0)->getId()==noeud->getId()) {
          brotherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getSon(1)->getNodeProperty(SPECIESID))->getValue());
        }
        else {
          brotherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getSon(0)->getNodeProperty(SPECIESID))->getValue());
        }
        
        if (fatherNum == num){
          //this is not the first duplication: we cannot compute the score yet
          int numGenes = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(NUMGENES))->getValue())+(dynamic_cast<const Number<int> *>(son1->getNodeProperty(NUMGENES))->getValue());
          changeNodeProperty(*(noeud),NUMGENES, Number<int>(numGenes)); 
          probaBranch = 0.0;
        }
        else {
          //this is the first duplication, we compute the probability  
          branchNumbers[num]++;
          int numGenes = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(NUMGENES))->getValue())+(dynamic_cast<const Number<int> *>(son1->getNodeProperty(NUMGENES))->getValue());
          probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], numGenes);
          if (numGenes==2) {
            num2lineages[num]++;
          }
          changeNodeProperty(*(noeud),NUMGENES, Number<int>(1));
        }
      }
      else {//A duplication at the root
        //this is the first duplication, we compute the probability  
        branchNumbers[num]++;
        int numGenes = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(NUMGENES))->getValue())+(dynamic_cast<const Number<int> *>(son1->getNodeProperty(NUMGENES))->getValue());
        probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], numGenes);
        if (numGenes==2) {
          num2lineages[num]++;
        }
        changeNodeProperty(*(noeud),NUMGENES, Number<int>(1));
      }
    }
    else { //no duplication, no loss on the branch leading to noeud
      
      if (noeud->hasFather()) {
        int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
        if (*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D") {
          probaBranch = 0.0;
          changeNodeProperty(*(noeud),NUMGENES, Number<int>(1)); 
          if (num!=fatherNum) {
            //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
            probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], 1);
            num1lineages[num]++;
            branchNumbers[num]++; 
          }
        }
        else {
          branchNumbers[num]++;
          probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], 1);
          num1lineages[num]++;
          changeNodeProperty(*(noeud),NUMGENES, Number<int>(1));
        }
      }
      else {
        branchNumbers[num]++;
        probaBranch = computeLogBranchProbability(duplicationProbabilities[num], lossProbabilities[num], 1);
        num1lineages[num]++;
        changeNodeProperty(*(noeud),NUMGENES, Number<int>(1));
      }
    }
    
    //if not at the root, there may be a long story from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.
    if (noeud->hasFather() && tree.getNode(num)->hasFather()) {
      int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
      if (num!=fatherNum) {
        int tempNum = num;
        Node * spTempNode = spNode;
        //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
        int tempNumBrother;
        while ((tempNum!=fatherNum) && (spTempNode->hasFather()) && (spTempNode->getFather()->getId()!=fatherNum)) {
          //There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication
          
          if (spTempNode->getFather()->getSon(0)->getId()==tempNum) {
            tempNumBrother = spTempNode->getFather()->getSon(1)->getId();
          }
          else {
            tempNumBrother = spTempNode->getFather()->getSon(0)->getId();
          }
          probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0);   	
          num0lineages[tempNumBrother]++;
          branchNumbers[tempNumBrother]++; 
          spTempNode = spTempNode->getFather();
          tempNum = spTempNode->getId();
          //No duplication on branch leading to node now named tempNum, and no loss either
          probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNum], lossProbabilities[tempNum], 1);
          num0lineages[tempNum]++;
          branchNumbers[tempNum]++;
        }
        if ((*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D")&&(num!=fatherNum)) {
          if (spTempNode->getFather()->getSon(0)->getId()==tempNum) {
            tempNumBrother = spTempNode->getFather()->getSon(1)->getId();
          }
          else {
            tempNumBrother = spTempNode->getFather()->getSon(0)->getId();
          }
          probaBranch+=computeLogBranchProbability(duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0);
          num0lineages[tempNumBrother]++;
          branchNumbers[tempNumBrother]++; 
        }
      }
    }
    
    changeBranchProperty(*(noeud), EVENTSPROBA, Number<double>(probaBranch));
    double lowlik = probaBranch+(dynamic_cast<const Number<double> *>(son0->getNodeProperty(LOWLIK))->getValue())+(dynamic_cast<const Number<double> *>(son1->getNodeProperty(LOWLIK))->getValue());
    changeNodeProperty(*(noeud), LOWLIK, Number<double>(lowlik));
  }
}



/**************************************************************************
 * This function writes a tree with a property whose format is integer.
 **************************************************************************/

std::string nodeToParenthesisWithIntNodeValues(const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName) throw (NodeNotFoundException)
{ 
  if(!tree.hasNode(nodeId)) throw NodeNotFoundException("nodeToParenthesisWithIntNodeValues", nodeId);
  std::ostringstream s;
  if(tree.isLeaf(nodeId))
  {
    s << tree.getNodeName(nodeId) << " "<< (dynamic_cast<const Number<int> *>(tree.getBranchProperty(nodeId, propertyName))->getValue());
  }
  else
  {
    s << "(";
    std::vector<int> sonsId = tree.getSonsId(nodeId);
    s << nodeToParenthesisWithIntNodeValues(tree, sonsId[0], bootstrap, propertyName);
    for(unsigned int i = 1; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesisWithIntNodeValues(tree, sonsId[i], bootstrap, propertyName);
    }
    s << ")";
   
    if(bootstrap)
    {
      if(tree.hasBranchProperty(nodeId, TreeTools::BOOTSTRAP))
        s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(nodeId, TreeTools::BOOTSTRAP))->getValue());
    }
    else
    {
      if(tree.hasBranchProperty(nodeId, propertyName))
        s << (dynamic_cast<const Number<int> *>(tree.getBranchProperty(nodeId, propertyName))->getValue());
    }
  }
  if(tree.hasDistanceToFather(nodeId)) s << ":" << tree.getDistanceToFather(nodeId);
  return s.str();  
}




std::string treeToParenthesisWithIntNodeValues(const Tree & tree, bool bootstrap, const std::string & propertyName)
{
  std::ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  std::vector<int> sonsId = tree.getSonsId(rootId);
  if(tree.isLeaf(rootId))
  {
    s << tree.getNodeName(rootId);
    for(unsigned int i = 0; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesisWithIntNodeValues(tree, sonsId[i], bootstrap, propertyName);
    }
  }
  else
  {
    s << nodeToParenthesisWithIntNodeValues(tree, sonsId[0], bootstrap, propertyName);
    for(unsigned int i = 1; i < sonsId.size(); i++)
    {
      s << "," << nodeToParenthesisWithIntNodeValues(tree, sonsId[i], bootstrap, propertyName);
    }
  }
  s << ")";
  if(bootstrap)
  {
    if(tree.hasBranchProperty(rootId, TreeTools::BOOTSTRAP))
      s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(rootId, TreeTools::BOOTSTRAP))->getValue());
  }
  else
  {
    if(tree.hasBranchProperty(rootId, propertyName))
      s << (dynamic_cast<const Number<int> *>(tree.getBranchProperty(rootId, propertyName))->getValue());
  }
  s << ";" << std::endl;
  return s.str();  
}

/**************************************************************************
 * This function writes a tree with a property whose format is double.
 **************************************************************************/

std::string nodeToParenthesisWithDoubleNodeValues(const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName) throw (NodeNotFoundException)
{ 
  if(!tree.hasNode(nodeId)) throw NodeNotFoundException("nodeToParenthesisWithDoubleNodeValues", nodeId);
  std::ostringstream s;
  if(tree.isLeaf(nodeId))
    {
      s << tree.getNodeName(nodeId) << " "<< (dynamic_cast<const Number<double> *>(tree.getBranchProperty(nodeId, propertyName))->getValue());
    }
  else
    {
      s << "(";
      std::vector<int> sonsId = tree.getSonsId(nodeId);
      s << nodeToParenthesisWithDoubleNodeValues(tree, sonsId[0], bootstrap, propertyName);
      for(unsigned int i = 1; i < sonsId.size(); i++)
        {
          s << "," << nodeToParenthesisWithDoubleNodeValues(tree, sonsId[i], bootstrap, propertyName);
        }
      s << ")";
      
      if(bootstrap)
        {
          if(tree.hasBranchProperty(nodeId, TreeTools::BOOTSTRAP))
            s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(nodeId, TreeTools::BOOTSTRAP))->getValue());
        }
      else
        {
          if(tree.hasBranchProperty(nodeId, propertyName))
            s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(nodeId, propertyName))->getValue());
        }
    }
  if(tree.hasDistanceToFather(nodeId)) s << ":" << tree.getDistanceToFather(nodeId);
  return s.str();  
}




std::string treeToParenthesisWithDoubleNodeValues(const Tree & tree, bool bootstrap, const std::string & propertyName)
{
  std::ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  std::vector<int> sonsId = tree.getSonsId(rootId);
  if(tree.isLeaf(rootId))
    {
      s << tree.getNodeName(rootId);
      for(unsigned int i = 0; i < sonsId.size(); i++)
        {
          s << "," << nodeToParenthesisWithDoubleNodeValues(tree, sonsId[i], bootstrap, propertyName);
        }
    }
  else
    {
      s << nodeToParenthesisWithDoubleNodeValues(tree, sonsId[0], bootstrap, propertyName);
      for(unsigned int i = 1; i < sonsId.size(); i++)
        {
          s << "," << nodeToParenthesisWithDoubleNodeValues(tree, sonsId[i], bootstrap, propertyName);
        }
    }
  s << ")";
  if(bootstrap)
    {
      if(tree.hasBranchProperty(rootId, TreeTools::BOOTSTRAP))
        s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(rootId, TreeTools::BOOTSTRAP))->getValue());
    }
  else
    {
      if(tree.hasBranchProperty(rootId, propertyName))
        s << (dynamic_cast<const Number<double> *>(tree.getBranchProperty(rootId, propertyName))->getValue());
    }
  s << ";" << std::endl;
  return s.str();  
}



/**************************************************************************
 * This function computes the score of a scenario of duplications and losses.
 * For this purpose, an exponential law is used to model change probabilities on branches.
 **************************************************************************/

/*
double computeScenarioScore (std::vector<int> lossNumbers, std::vector< double> lossProbabilities, std::vector< int> duplicationNumbers, std::vector < double> duplicationProbabilities)
{
  double score = 0;
  for (int i =0; i< lossNumbers.size() ; i++) {
    if (lossNumbers[i] >0){
      // std::cout << "i : "<<i<<" lossNumber : "<<lossNumbers[i] << "lossProbabilities[i] : "<<lossProbabilities[i]<< std::endl;
      score += lossNumbers[i]*log(lossProbabilities[i]);
    }
    else {
      score += log (1-2*lossProbabilities[i])-log (1-lossProbabilities[i]);
    }

    if (duplicationNumbers[i] >0){
      // std::cout << "i : "<<i<<" duplicationNumber : "<<duplicationNumbers[i] << "duplicationProbabilities[i] : "<<duplicationProbabilities[i]<<std::endl;
      score += duplicationNumbers[i]*log(duplicationProbabilities[i]);
    }
    else {
      score += log (1-2*duplicationProbabilities[i])-log (1-duplicationProbabilities[i]);
    }

  }
  
  return score;
}

*/

/**************************************************************************
 * This function sets the LOSSES and DUPLICATIONS values on the species tree, given std::vectors of duplication and loss numbers
 **************************************************************************/
void setLossesAndDuplications(TreeTemplate<Node> & tree, 
                              std::vector <int> &lossNumbers, 
                              std::vector <int> &duplicationNumbers) {
  std::vector< int > nodesIds = tree.getNodesId ();

  for (unsigned int i=0; i<nodesIds.size(); i++) {
    if(tree.hasBranchProperty(nodesIds[i], LOSSES)) {
      tree.getNode(nodesIds[i])->deleteBranchProperty(LOSSES);
    }
    tree.getNode(nodesIds[i])->setBranchProperty(LOSSES, BppString(TextTools::toString((lossNumbers[nodesIds[i]]))));
    if(tree.hasBranchProperty(nodesIds[i], DUPLICATIONS)) {
      tree.getNode(nodesIds[i])->deleteBranchProperty(DUPLICATIONS);
    }
    tree.getNode(nodesIds[i])->setBranchProperty(DUPLICATIONS, BppString(TextTools::toString((duplicationNumbers[nodesIds[i]]))));
  }
}

/**************************************************************************
 * This function sets the NUM0LINEAGES, NUM1LINEAGES and NUM2LINEAGES values on the species tree, given std::vectors
 **************************************************************************/
void assignNumLineagesOnSpeciesTree(TreeTemplate<Node> & tree, 
                              std::vector <int> &num0Lineages, 
                              std::vector <int> &num1Lineages, 
                              std::vector <int> &num2Lineages) {
  std::vector< int > nodesIds = tree.getNodesId ();
  std::string nums = "";
  for (unsigned int i=0; i<nodesIds.size(); i++) {
    if(tree.hasBranchProperty(nodesIds[i], NUMLINEAGES)) {
      tree.getNode(nodesIds[i])->deleteBranchProperty(NUMLINEAGES);
    }
    if(tree.hasNodeProperty(nodesIds[i], NUMLINEAGES)) {
      tree.getNode(nodesIds[i])->deleteNodeProperty(NUMLINEAGES);
    }
    nums = TextTools::toString(num0Lineages[nodesIds[i]])+"_"+TextTools::toString(num1Lineages[nodesIds[i]])+"_"+TextTools::toString(num2Lineages[nodesIds[i]]);
    tree.getNode(nodesIds[i])->setBranchProperty(NUMLINEAGES, BppString(nums));
    if (tree.getNode(nodesIds[i])->isLeaf() ) {
      tree.getNode(nodesIds[i])->setName( tree.getNode(nodesIds[i])->getName() + "_"+nums) ;
    }
  }
}


/*****************************************************************************
 * This function returns the speciesID assigned to a leaf.
 * 
 ****************************************************************************/


int assignSpeciesIdToLeaf(Node * node,  const std::map<std::string, std::string > & seqSp,  
                          const std::map<std::string, int > & spID) {
  std::map<std::string, std::string >::const_iterator seqtosp;
  seqtosp=seqSp.find(node->getName());
  if (seqtosp!=seqSp.end()){
    std::map<std::string, int >::const_iterator sptoid;
    sptoid = spID.find(seqtosp->second);
    if (sptoid!=spID.end()) {
      return(sptoid->second);
    }
    else {
      std::cerr <<"Error in assignSpeciesIdToLeaf: "<< seqtosp->second <<" not found in std::map spID for sequence "<< seqtosp->first<<std::endl;
        MPI::COMM_WORLD.Abort(1);
      exit(-1);
    }
  }
  else {
    std::cerr <<"Error in assignSpeciesIdToLeaf: "<< node->getName() <<" not found in std::map seqSp"<<std::endl;
      MPI::COMM_WORLD.Abort(1);
    exit(-1);
  }
}

/*****************************************************************************
 * This function recovers gene losses by comparing a subtree in a gene tree to
 * a species tree (tree).
 * 
 ****************************************************************************/

void recoverLosses(Node *& node, int & a, const int & b, int & olda, int & a0, 
                   const TreeTemplate<Node> & tree, 
                   double & likelihoodCell, 
                   const std::vector< double> & lossRates, 
                   const std::vector< double> & duplicationRates, 
                   int & dupData) {
 // const Node* const node = tree.getNode(a);
  olda=a;
  Node* nodeA;
  if (node->hasFather()) {
    nodeA = node->getFather();
  }
  else {
   std::cout <<"Problem in recoverLosses, nodeA has no father"<<std::endl; 
  }
  //a = node->getFather()->getId();
 // std::cout <<"a "<<a <<std::endl;
  a = nodeA->getId();
 // std::cout <<"a "<<a <<std::endl;

  node = nodeA;
 //  std::cout <<"a "<<a <<std::endl;
  //std::cout <<"here 10"<<std::endl;

 /* std::cout <<"nodeagetson0getId: "<<nodeA->getSon(0)->getId()<<std::endl;
  std::cout <<"nodeagetson1getId: "<<nodeA->getSon(1)->getId()<<std::endl;

   std::vector <int> ids = tree.getNodesId();
  for (int i =0; i<ids.size() ; i++) {
    std::cout <<"ids :"<<ids[i]<<std::endl;
  }

  std::cout <<"hereheh\n"<<std::endl;*/

  std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(*(nodeA->getSon(0)));
  nodesIds0.push_back(nodeA->getSon(0)->getId());
  std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(*(nodeA->getSon(1)));
  nodesIds1.push_back(nodeA->getSon(1)->getId());
  int lostNodeId = -1;
 // std::cout <<"here 9"<<std::endl;

  if ((nodeA->getSon(0)->getId()==olda)&&(!(VectorTools::contains(nodesIds1, b)))&&(b!=a))
    {
    //  std::cout <<"here 8"<<std::endl;
      lostNodeId=nodeA->getSon(1)->getId();
      //UNDONE TEST 1009
    //  if (!nodeA->getSon(1)->isLeaf()) {
      likelihoodCell += computeLogBranchProbability(duplicationRates[lostNodeId], lossRates[lostNodeId], 0);
     // }
     // std::cout <<"recoverLosses 1 lostNodeId"<<lostNodeId<<std::endl;
    }
  else  if ((nodeA->getSon(1)->getId()==olda)&&(!(VectorTools::contains(nodesIds0, b)))&&(b!=a))
    {
    //  std::cout <<"here 7"<<std::endl;

      lostNodeId=nodeA->getSon(0)->getId();
      //UNDONE TEST 1009
     // if (!nodeA->getSon(0)->isLeaf()) {
      likelihoodCell += computeLogBranchProbability(duplicationRates[lostNodeId], lossRates[lostNodeId], 0);
     // }
    //  std::cout <<"recoverLosses 2 lostNodeId"<<lostNodeId<<std::endl;
    }
 // if (lostNodeId!= -1) {
    
  //  std::cout << "A loss on branch "<<lostNodeId<<std::endl;
//  }
/*  if ((dupData>0)&&(olda==a0)) {
    likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], dupData);
    std::cout <<dupData<<" B genes on branch "<<olda<<std::endl;
   // dupData = 0; //Resetting the dupdata value
  }
  else {*/
  if ((olda!=a0)) {
    likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);

   // std::cout <<"1 I genes on branch "<<olda<<std::endl;
  }
 // std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
 
  return;
 }

/*****************************************************************************
 * This function recovers gene losses by comparing a subtree in a gene tree to
 * a species tree (tree), when a duplication has affected the subtree.
 * 
 ****************************************************************************/

void recoverLossesWithDuplication(const Node * nodeA, 
                                  const int &a, 
                                  const int &olda, 
                                  const TreeTemplate<Node> & tree,
                                  double & likelihoodCell, 
                                  const std::vector< double> & lossRates, 
                                  const std::vector< double> & duplicationRates) {
  //The loss has occured before a0
 // const Node * nodeA = tree.getNode(a);
  const Node * nodeOldA ;
  const Node * lostNode;
  if (nodeA->getSon(0)->getId() == olda) {
    nodeOldA = nodeA->getSon(0);
    lostNode=nodeOldA->getFather()->getSon(1);
  }
  else {
    nodeOldA = nodeA->getSon(1);
    lostNode=nodeOldA->getFather()->getSon(0);
  }
  //We need to place the loss event in the right lineage  
  //UNDONE TEST 1009
 // if(!lostNode->isLeaf()) {
  likelihoodCell += computeLogBranchProbability(duplicationRates[lostNode->getId()], lossRates[lostNode->getId()], 0);
 // }
// std::cout <<"recoverLossesWithDuplication loss on branch "<<lostNode->getId()<<std::endl; 
 /* likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
  std::cout << "H 1 gene on branch "<<olda<<std::endl;*/
//  std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
  return;
}

/*****************************************************************************
 * This function computes the lower conditional likelihood of a subtree and 
 * assigns its summit node a species ID. Notations are influenced by 
 * Zmasek and Eddy algorithm (2001).
 * 
 ****************************************************************************/


double computeConditionalLikelihoodAndAssignSpId(TreeTemplate<Node> & tree,
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
                                                 bool atRoot) {
  if (rootLikelihood == 0.0) {
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
   /* std::cout <<"son0spid : "<<son0SpId<<std::endl;
    std::cout <<"son1spid : "<<son1SpId<<std::endl;
*/
    
  /*  std::vector <int> ids = tree.getNodesId();
    for (int i =0; i<ids.size() ; i++) {
      std::cout <<"ids :"<<ids[i]<<std::endl;
    }
    */
    
    Node * temp0 = tree.getNode(son0SpId);
    Node * temp1 = tree.getNode(son1SpId);

    while (a!=b) { //There have been losses !
      if (a>b) {

      //  std::cout <<"before recoverLosses temp0id: "<<temp0.getId()<<std::endl;
        recoverLosses(temp0, a, b, olda, a0, tree, rootLikelihood, lossRates, duplicationRates, son0DupData);
      /*  std::cout <<"HERE"<<std::endl;
        std::cout <<"after recoverLosses temp0id: "<<temp0.getId()<<std::endl;*/
      }
      else {

/*         std::cout <<"before recoverLosses temp1id: "<<temp1.getId()<<std::endl;
        std::cout <<"HEHEHEEHEH\n\n"<<std::endl;*/
        recoverLosses(temp1, b, a, oldb, b0, tree, rootLikelihood, lossRates, duplicationRates, son1DupData);
        // std::cout <<"after recoverLosses temp1id: "<<temp1.getId()<<std::endl;

      }
    }
    rootSpId = a;
   // std::cout <<"rootspid : "<<rootSpId<<std::endl;
    if ((a==a0) || (b==b0)) //There has been a duplication !
      {
        if ((a==a0) && (b==b0)) {
          rootDupData += son0DupData+son1DupData;
          rootLikelihood-=(computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData) + 
                           computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData));
        }//there has been no loss, here
        else if (b==b0) { //The loss has occured before a0
          rootDupData += son1DupData+1;
          rootLikelihood-=computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
          recoverLossesWithDuplication(temp0, a, olda, tree, rootLikelihood, lossRates, duplicationRates);

        //  if (son0DupData>0) {
          /*  std::cout <<"C On branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
            son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);*/
         // }
        }
        else { //The loss has occured before b0

          rootDupData += son0DupData+1;
          rootLikelihood-=computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
          recoverLossesWithDuplication(temp1, b, oldb, tree, rootLikelihood, lossRates, duplicationRates);
        //  if (son1DupData>0) {
          /*  std::cout <<"D On branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
            son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);*/
        //  }
        }
        //Counting the duplication(s) event(s)
      /*  if (a==a0) {
          rootDupData += son0DupData+1;
        }
        else if (son0DupData>1) {
          std::cout <<"duplication on branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
          son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
        }
        if (b==b0) {
          rootDupData += son1DupData+1;
        }
        else if (son1DupData>1) {
          std::cout <<"duplication on branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
          son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
        }*/
        //if at the root, we need to compute the contribution of the 
        //duplication event to the lower conditional likelihood now.
      //  if (atRoot){
         // std::cout <<"E  duplication on branch "<<a<< "number of genes :"<<rootDupData<<std::endl;
        if(atRoot) {
          rootLikelihood += computeLogBranchProbabilityAtRoot(duplicationRates[a], lossRates[a], rootDupData);

        }
        else {
          rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], rootDupData);
         // std::cout<<"rootLikelihood5: "<< rootLikelihood<<std::endl;  

        }
      //  }
      }
    else //there was no duplication
      {
        // rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], 1);
        //Perhaps there have been duplications that need to be counted in son nodes        
      //  if (son1DupData>0) {
         // std::cout << "On branch "<<b0<<"num dup"<<son1DupData<<std::endl;
       //   son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
      //  }
      //  if (son0DupData>0) {
         // std::cout << "On branch "<<a0<<"num dup"<<son0DupData<<std::endl;
       //   son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
       // }
        rootDupData = 1;
       // if (atRoot){
         // std::cout <<"F duplication on branch "<<a<< "number of genes :"<<rootDupData<<std::endl;
        if (atRoot) {
          rootLikelihood += computeLogBranchProbabilityAtRoot(duplicationRates[a], lossRates[a], rootDupData);

        }else {
          rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], rootDupData);
        }

       // }
      }
    //Setting the lower conditional likelihood for the node of interest.
    rootLikelihood += son0Likelihood + son1Likelihood;

   // std::cout <<"HERErootlikelihood "<<rootLikelihood<<std::endl;
  }
  else {
    //std::cout << "Error in computeLowerConditionalLikelihoodAndAssignSpId: initial conditional likelihood != 0."<<std::endl;
  }

  return(rootLikelihood);
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

double computeSubtreeLikelihoodPostorder(TreeTemplate<Node> & spTree, 
                                         TreeTemplate<Node> & geneTree, 
                                         Node * node, 
                                         const std::map<std::string, std::string > & seqSp, 
                                         const std::map<std::string, int > & spID, 
                                         std::vector <std::vector<double> > & likelihoodData, 
                                         const std::vector< double> & lossRates, 
                                         const std::vector < double> & duplicationRates, 
                                         std::vector <std::vector<int> > & speciesIDs, 
                                         std::vector <std::vector<int> > & dupData) {
	int id=node->getId();
 	if (node->isLeaf()) {

		if (likelihoodData[id][0]==0.0) {
      speciesIDs[id][0]=speciesIDs[id][1]=speciesIDs[id][2]=assignSpeciesIdToLeaf(node, seqSp, spID);
      likelihoodData[id][0]=likelihoodData[id][1]=likelihoodData[id][2]=computeLogBranchProbability(duplicationRates[speciesIDs[id][0]], lossRates[speciesIDs[id][0]], 1);
      dupData[id][0] = dupData[id][1] = dupData[id][2] = 1;
 
          //  std::cout <<"leafLk id "<< id << " lk: "<<likelihoodData[id][0]<<std::endl;
    
        }
  /*  std::cout <<"at leaf "<<node->getName()<<std::endl;
    std::cout <<"lk "<<likelihoodData[id][0]<<std::endl;*/

    return(likelihoodData[id][0]);
  }
  else {

    std::vector <Node *> sons = node->getSons();
#ifdef _OPENMP
#pragma omp single
	  {
#pragma omp task untied
      computeSubtreeLikelihoodPostorder(spTree, geneTree, 
                                        sons[0], seqSp, 
                                        spID, likelihoodData, 
                                        lossRates, duplicationRates, 
                                        speciesIDs, dupData);
#pragma omp task untied
      computeSubtreeLikelihoodPostorder(spTree, geneTree, 
                                        sons[1], seqSp, 
                                        spID, likelihoodData, 
                                        lossRates, duplicationRates, 
                                        speciesIDs, dupData);
#pragma omp taskwait
	  }
#else	  
    for (unsigned int i = 0; i< sons.size(); i++){
      computeSubtreeLikelihoodPostorder(spTree, geneTree, 
                                        sons[i], seqSp, 
                                        spID, likelihoodData, 
                                        lossRates, duplicationRates, 
                                        speciesIDs, dupData);

    }
#endif
    
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0, directionSon1;
    std::vector <Node *> neighbors = sons[0]->getNeighbors();
    for (unsigned int i=0; i<neighbors.size(); i++) {
      if (neighbors[i]==node) {
        directionSon0 = i;
      }
    }
    neighbors = sons[1]->getNeighbors();
    for (unsigned int i=0; i<neighbors.size(); i++) {
      if (neighbors[i]==node) {
        directionSon1 = i;
      }
    }

  /*  neighbors = node->getNeighbors();
    for (unsigned int i=0; i<neighbors.size(); i++) {
      if (neighbors[i]==node) {
        directionSon1 = i;
      }
    }
    
      else if (neighbors[i]==sons[1]) {
        directionSon1 = i;
      }
      else if (neighbors[i]==node->getFather()) {
        directionFather = i;
      }
    }
   
    
    
    std::cout << "fatherDirection "<<directionFather<<std::endl;*/
 /*   
  std::cout << "son 0 lk "<<likelihoodData[idSon0][directionSon0]<< " directionSon0 "<<  directionSon0<<std::endl;
   std::cout << "son 1 lk "<<likelihoodData[idSon1][directionSon1]<<" directionSon1 "<<  directionSon1<<std::endl;
    std::cout <<"node ID "<<id<<"isRoot? "<<TreeTemplateTools::isRoot(*node)<<std::endl;*/
      
    computeConditionalLikelihoodAndAssignSpId(spTree, sons, 
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
                                              TreeTemplateTools::isRoot(*node));

     // std::cout <<"father lk "<< likelihoodData[id][0]<<std::endl;

      return(likelihoodData[id][0]);
	}
	
	
}





/*****************************************************************************
 * This function computes the likelihood of a rooting. 
 * It is called by the preorder tree traversal.
 * "direction" determines the branch we're on: it is the branch leading to the
 * "direction"th son of node. direction = sonNumber+1
 * 
 ****************************************************************************/
void computeRootingLikelihood(TreeTemplate<Node> & spTree, 
                              Node * node, 
                              std::vector <std::vector<double> > & likelihoodData, 
                              const std::vector< double> & lossRates, 
                              const std::vector < double> & duplicationRates, 
                              std::vector <std::vector<int> > & speciesIDs, 
                              std::vector <std::vector<int> > & dupData, 
                              int sonNumber, 
                              std::map <double, Node*> & LksToNodes) {
  int geneNodeId = node->getId();
   
  int directionForFather;
  std::vector <Node*> nodes;
  nodes.push_back(node->getFather());
  if (sonNumber==0) { //If sonNumber==0, the subtree we're interested in is composed of son 1 and father of node.
    nodes.push_back(node->getSon(1));
  }
  else { //If sonNumber==1, the subtree we're interested in is composed of son 0 and father of node.
    nodes.push_back(node->getSon(0));
  }
  if (node->getFather()->getSon(0)==node) {
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

/*  std::cout <<"computeRootingLikelihood : nodeID "<<geneNodeId<<"; sonID "<<idNode1<<std::endl;

  std::cout <<"likelihoodData[Father] "<< likelihoodData[idNode0][directionNode0]<<std::endl;
  std::cout <<"likelihoodData[Son] "<< likelihoodData[idNode1][directionNode1]<<std::endl;
*/
  computeConditionalLikelihoodAndAssignSpId(spTree, nodes, 
                                            likelihoodData[geneNodeId][sonNumber+1], 
                                            likelihoodData[idNode0][directionNode0], 
                                            likelihoodData[idNode1][directionNode1], 
                                            lossRates, duplicationRates, 
                                            speciesIDs[geneNodeId][sonNumber+1], 
                                            speciesIDs[idNode0][directionNode0], 
                                            speciesIDs[idNode1][directionNode1], 
                                            dupData[geneNodeId][sonNumber+1], 
                                            dupData[idNode0][directionNode0], 
                                            dupData[idNode1][directionNode1], false);
  //Now we have the conditional likelihood of the upper subtree, 
  //as well as the conditional likelihood of the lower subtree (which we already had)
  //We can thus compute the total likelihood of the rooting.
  
  std::vector <Node*> sons;
  sons.push_back(node);
 
  sons.push_back(node->getSon(sonNumber));
  int idSon0 = geneNodeId;
  int idSon1 = sons[1]->getId();
  unsigned int directionSon0, directionSon1;
  directionSon0 = sonNumber+1;
  directionSon1 = 0;
   
  double rootLikelihood = 0.0;
  int rootSpId;
  int rootDupData = 0;
  
  computeConditionalLikelihoodAndAssignSpId(spTree, sons, 
                                            rootLikelihood, 
                                            likelihoodData[idSon0][directionSon0], 
                                            likelihoodData[idSon1][directionSon1], 
                                            lossRates, duplicationRates, 
                                            rootSpId, speciesIDs[idSon0][directionSon0], 
                                            speciesIDs[idSon1][directionSon1], 
                                            rootDupData, dupData[idSon0][directionSon0], 
                                            dupData[idSon1][directionSon1], true);
 // std::cout <<"LK FOUND "<<rootLikelihood<<std::endl;
  while (LksToNodes.find(rootLikelihood)!=LksToNodes.end()) {
   // std::cout <<"changing rootLikelihood !!!!!!!!!!!!!!!!!!!"<<std::endl;
    rootLikelihood+=SMALLPROBA;
  }
  LksToNodes[rootLikelihood]=node->getSon(sonNumber);
}





/*****************************************************************************
 * This function performs a preorder tree traversal in order to find likelihoods for rootings. 
 * When used after the postorder tree traversal function, likelihoods for all rootings are computed.
 * likelihoodData contains all lower conditional likelihoods for all nodes.
 * speciesIDs contains all species IDs for all nodes.
 * 
 ****************************************************************************/

void computeSubtreeLikelihoodPreorder(TreeTemplate<Node> & spTree, 
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
                                      std::map <double, Node*> & LksToNodes) {
	computeRootingLikelihood(spTree, node, 
                           likelihoodData, lossRates, 
                           duplicationRates, speciesIDs, 
                           dupData, sonNumber, LksToNodes);
  if (node->isLeaf()) {
    return; 
  }
  Node * son;
  if (sonNumber==1) {
    son= node->getSon(1);
  }
  else {
    son= node->getSon(0);
  }
//  for (int i = 0; i< sons.size(); i++){
#ifdef _OPENMP
#pragma omp single
	{
#endif
    for (unsigned int j =0; j<son->getNumberOfSons(); j++) {
#ifdef _OPENMP
#pragma omp task
#endif
      computeSubtreeLikelihoodPreorder(spTree, geneTree, 
                                       son, seqSp, spID, 
                                       likelihoodData, 
                                       lossRates, duplicationRates, 
                                       speciesIDs, dupData, j, LksToNodes);
    }
# ifdef _OPENMP
#pragma omp taskwait
	}
#endif
//  }
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




void recoverLossesAndLineages(Node *& node, int & a, const int & b, int & olda, int & a0, 
                   const TreeTemplate<Node> & tree, 
                   int & dupData, std::vector<int> &num0lineages, std::vector<int> &num1lineages) {
 // const Node* const node = tree.getNode(a);
 
  olda=a;
  Node* nodeA;
  if (node->hasFather()) {
    nodeA = node->getFather();
  }
  else {
    std::cout <<"Problem in recoverLossesAndLineages, nodeA has no father"<<std::endl; 
  }
  a = nodeA->getId();
  node = nodeA;
  std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(*(nodeA->getSon(0)));
  nodesIds0.push_back(nodeA->getSon(0)->getId());
  std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(*(nodeA->getSon(1)));
  nodesIds1.push_back(nodeA->getSon(1)->getId());
  int lostNodeId = -1;
    
  if ((nodeA->getSon(0)->getId()==olda)&&(!(VectorTools::contains(nodesIds1, b)))&&(b!=a))
    {
    
      lostNodeId=nodeA->getSon(1)->getId();
    }
  else  if ((nodeA->getSon(1)->getId()==olda)&&(!(VectorTools::contains(nodesIds0, b)))&&(b!=a))
    {
      
      lostNodeId=nodeA->getSon(0)->getId();
    }
  if (lostNodeId!= -1) {
    //likelihoodCell += computeLogBranchProbability(duplicationRates[lostNodeId], lossRates[lostNodeId], 0);
    num0lineages[lostNodeId]+=1;
    //  std::cout << "A loss on branch "<<lostNodeId<<std::endl;
  }
  /*  if ((dupData>0)&&(olda==a0)) {
   likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], dupData);
   std::cout <<dupData<<" B genes on branch "<<olda<<std::endl;
   // dupData = 0; //Resetting the dupdata value
   }
   else {*/
  if ((olda!=a0)) {
    //likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
    num1lineages[olda]+=1;
    // std::cout <<"1 I genes on branch "<<olda<<std::endl;
  }
  // std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
  
  return;
}




/****************************************************************************/



void recoverLossesAndLineagesWithDuplication(const Node * nodeA, 
                                             const int &a, 
                                             const int &olda, 
                                             const TreeTemplate<Node> & tree, 
                                             std::vector <int> &num0lineages) {
  //The loss has occured before a0
//  const Node * nodeA = tree.getNode(a);
  const Node * nodeOldA ;
  const Node * lostNode;
  if (nodeA->getSon(0)->getId() == olda) {
    nodeOldA = nodeA->getSon(0);
    lostNode=nodeOldA->getFather()->getSon(1);
  }
  else {
    nodeOldA = nodeA->getSon(1);
    lostNode=nodeOldA->getFather()->getSon(0);
  }
  //We need to place the loss event in the right lineage  
  //likelihoodCell += computeLogBranchProbability(duplicationRates[lostNode->getId()], lossRates[lostNode->getId()], 0);
  num0lineages[lostNode->getId()]+=1;
  // std::cout <<"G loss on branch "<<lostNode->getId()<<std::endl; 
  /* likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
   std::cout << "H 1 gene on branch "<<olda<<std::endl;*/
  //  std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
  return;
}





/****************************************************************************/




void computeNumbersOfLineagesInASubtree(TreeTemplate<Node> & tree,
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
                                          std::set <int> &branchesWithDuplications) {
  int a, a0, olda;
  int b, b0, oldb;
  a = a0 = olda = son0SpId;
  b = b0 = oldb = son1SpId;
  
  Node * temp0 = tree.getNode(son0SpId);
  Node * temp1 = tree.getNode(son1SpId);

  
  while (a!=b) { //There have been losses !
    if (a>b) {
      recoverLossesAndLineages(temp0, a, b, olda, a0, tree, son0DupData, num0lineages, num1lineages);
    }
    else {
      recoverLossesAndLineages(temp1, b, a, oldb, b0, tree, son1DupData, num0lineages, num1lineages);
    }
  }
  rootSpId = a;
  if ((a==a0) || (b==b0)) //There has been a duplication !
    {
   // std::cout <<"INSERTING"<<std::endl;
      branchesWithDuplications.insert(sons[0]->getFather()->getId());
      branchesWithDuplications.insert(sons[0]->getId()); //We also include the two son nodes themselves
      branchesWithDuplications.insert(sons[1]->getId());
      branchesWithDuplications.insert(sons[0]->getId()); 
      if (sons[1]->getNumberOfSons () > 0 ) //and the grandsons, if present.
        {
        branchesWithDuplications.insert(sons[1]->getSon(0)->getId() );
        branchesWithDuplications.insert(sons[1]->getSon(1)->getId() );
        }
      if (sons[0]->getNumberOfSons () > 0 )
        {
        branchesWithDuplications.insert(sons[0]->getSon(0)->getId() );
        branchesWithDuplications.insert(sons[0]->getSon(1)->getId() );
        }
    
      if ((a==a0) && (b==b0)) {
        rootDupData += son0DupData+son1DupData;
        /* rootLikelihood-=(computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData) + 
         computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData));*/
        if (son0DupData==1) {
          num1lineages[a0]=num1lineages[a0]-1;
        }
        else if (son0DupData>=2) {//All branches with 2 or more lineages
          num2lineages[a0]=num2lineages[a0]-1;
        }
        if (son1DupData==1) {
          num1lineages[b0]=num1lineages[b0]-1;
        }
        else if (son1DupData>=2) {//All branches with 2 or more lineages
          num2lineages[b0]=num2lineages[b0]-1;
        }
        
      }//there has been no loss, here
      else if (b==b0) { //The loss has occured before a0
        rootDupData += son1DupData+1;
        //rootLikelihood-=computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
        if (son1DupData==1) {
          num1lineages[b0]=num1lineages[b0]-1;
        }
        else if (son1DupData>=2) {//All branches with 2 or more lineages
          num2lineages[b0]=num2lineages[b0]-1;
        }
        recoverLossesAndLineagesWithDuplication(temp0, a, olda, tree, num0lineages);
       // if (son0DupData>0) {
          /*  std::cout <<"C On branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
           son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);*/
      //  }
      }
      else { //The loss has occured before b0
        rootDupData += son0DupData+1;
        //rootLikelihood-=computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
        if (son0DupData==1) {
          num1lineages[a0]=num1lineages[a0]-1;
        }
        else if (son0DupData>=2) {//All branches with 2 or more lineages
          num2lineages[a0]=num2lineages[a0]-1;
        }
        recoverLossesAndLineagesWithDuplication(temp1, b, oldb, tree, num0lineages);
      //  if (son1DupData>0) {
          /*  std::cout <<"D On branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
           son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);*/
      //  }
      }
      //Counting the duplication(s) event(s)
      /*  if (a==a0) {
       rootDupData += son0DupData+1;
       }
       else if (son0DupData>1) {
       std::cout <<"duplication on branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
       son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
       }
       if (b==b0) {
       rootDupData += son1DupData+1;
       }
       else if (son1DupData>1) {
       std::cout <<"duplication on branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
       son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
       }*/
      //if at the root, we need to compute the contribution of the 
      //duplication event to the lower conditional likelihood now.
      //  if (atRoot){
      // std::cout <<"E  duplication on branch "<<a<< "number of genes :"<<rootDupData<<std::endl;
      //rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], rootDupData);
      if (rootDupData==1) {
        num1lineages[rootSpId]+=1;
      }
      else if (rootDupData>=2){//All branches with 2 or more lineages
        num2lineages[rootSpId]+=1;
      }
      //  }
    }
  else //there was no duplication
    {
      // rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], 1);
      //Perhaps there have been duplications that need to be counted in son nodes        
     // if (son1DupData>0) {
        // std::cout << "On branch "<<b0<<"num dup"<<son1DupData<<std::endl;
        //   son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);
     // }
     // if (son0DupData>0) {
        // std::cout << "On branch "<<a0<<"num dup"<<son0DupData<<std::endl;
        //   son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);
    //  }
      rootDupData = 1;
      // if (atRoot){
      // std::cout <<"F duplication on branch "<<a<< "number of genes :"<<rootDupData<<std::endl;
      // rootLikelihood += computeLogBranchProbability(duplicationRates[a], lossRates[a], rootDupData);
      num1lineages[rootSpId]+=1;
      // }
    }
  //Setting the lower conditional likelihood for the node of interest.
  // rootLikelihood += son0Likelihood + son1Likelihood;
  // std::cout <<"HERErootlikelihood "<<rootLikelihood<<std::endl;

  return;
}



/****************************************************************************/


void computeNumbersOfLineagesFromRoot(TreeTemplate<Node> * spTree, 
                                      TreeTemplate<Node> * geneTree, 
                                      Node * node, 
                                      std::map<std::string, std::string > seqSp,
                                      std::map<std::string, int > spID,
                                      std::vector <int> &num0lineages, 
                                      std::vector <int> &num1lineages, 
                                      std::vector <int> &num2lineages, 
                                      std::vector <std::vector<int> > & speciesIDs,
                                      std::vector <std::vector<int> > & dupData, 
                                      std::set <int> & branchesWithDuplications) {
  int id=node->getId();
 	if (node->isLeaf()) {
    speciesIDs[id][0]=speciesIDs[id][1]=speciesIDs[id][2]=assignSpeciesIdToLeaf(node, seqSp, spID);
    num1lineages[speciesIDs[id][0]]+=1;
    dupData[id][0] = dupData[id][1] = dupData[id][2] = 1;
    return;
  }
  else {
    std::vector <Node *> sons = node->getSons();
#ifdef _OPENMP
#pragma omp single 
	  {
#endif
    for (unsigned int i = 0; i< sons.size(); i++){
#ifdef _OPENMP
#pragma omp task untied
#endif
		computeNumbersOfLineagesFromRoot(spTree, geneTree, sons[i], 
                                       seqSp, spID, num0lineages, 
                                       num1lineages, num2lineages, 
                                       speciesIDs, dupData, 
                                       branchesWithDuplications);
    }
#ifdef _OPENMP
#pragma omp wait untied
	  }
#endif	  
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0, directionSon1;
    std::vector <Node *> neighbors = sons[0]->getNeighbors();
    for (unsigned int i=0; i<neighbors.size(); i++) {
      if (neighbors[i]==node) {
        directionSon0 = i;
      }
    }
    neighbors = sons[1]->getNeighbors();
    for (unsigned int i=0; i<neighbors.size(); i++) {
      if (neighbors[i]==node) {
        directionSon1 = i;
      }
    }

    computeNumbersOfLineagesInASubtree(*spTree, sons, 
                                       speciesIDs[id][0], 
                                       speciesIDs[idSon0][directionSon0], 
                                       speciesIDs[idSon1][directionSon1], 
                                       dupData[id][0], dupData[idSon0][directionSon0], 
                                       dupData[idSon1][directionSon1], 
                                       TreeTemplateTools::isRoot(*node), 
                                       num0lineages, num1lineages, 
                                       num2lineages, branchesWithDuplications);
      return;
	}
  
  
  
}


/*****************************************************************************
  * This function aims at finding the most likely reconciliation, 
  * using a double recursive tree traversal. 
  * The first traversal is post-order, and then the second traversal is pre-order.
  * This is a modification of an algorithm quickly explained in 
  * Chen, Durand, Farach-Colton, J. Comp. Biol. pp429-447, 2000.
  * Conditional likelihoods are recorded in a table. 
  * This table has (number of nodes) elements, and for each node, 
  * contains three conditional likelihoods. 
  * The table is thus (number of nodes)*3 cells. For each node i, 
  * likelihoodData[i][j] contains the conditional likelihood of the subtree 
  * having its root in subtree opposite neighbour j of node i.
  * Node species IDs are also recorded in a (number of nodes)*3 cells table.
  * The boolean "fillTables" is here to tell whether we want to update the vectors num*lineages.
  ****************************************************************************/

double findMLReconciliationDR (TreeTemplate<Node> * spTree, 
                               TreeTemplate<Node> * geneTree, 
                               std::map<std::string, std::string > seqSp,
                               std::map<std::string, int > spID,
                               std::vector< double> lossRates, 
                               std::vector < double> duplicationRates, 
                               int & MLindex, 
                               std::vector <int> &num0lineages, 
                               std::vector <int> &num1lineages, 
                               std::vector <int> &num2lineages, 
                               std::set <int> &nodesToTryInNNISearch, 
                               bool fillTables)
{
/*    std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis (*spTree, true)<<std::endl;
  */  
    
	if (!geneTree->isRooted()) {
		std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
		std::cout <<"!!!!!!gene tree is not rooted in findMLReconciliationDR !!!!!!"<<std::endl;
        MPI::COMM_WORLD.Abort(1);
		exit(-1);
      }		
	std::vector <double> nodeData(3, 0.0);
	std::vector <std::vector<double> > likelihoodData(geneTree->getNumberOfNodes(), nodeData);

  std::vector <int> nodeSpId(3, 0);
  std::vector <std::vector<int> > speciesIDs(geneTree->getNumberOfNodes(), nodeSpId);
  std::vector <std::vector<int> > dupData = speciesIDs;
  
  double initialLikelihood;
  //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
  std::map <double, Node*> LksToNodes;
  
// std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
//  std::cout << "IN findMLReconciliationDR "<<TreeTemplateTools::treeToParenthesis (*spTree, true)<<std::endl;
  
//	std::cout <<"CLOCKBEFOREPOSTORDER "<<clock()<<std::endl;
	Node * geneRoot = geneTree->getRootNode();
    
/*  std::cout <<"root number :"<<geneRoot->getId()<<std::endl;
  std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;*/
# ifdef _OPENMP
	double START, END;
	START= omp_get_wtime();
#pragma omp parallel
	{
		#pragma omp single 
		{
#endif

  initialLikelihood = computeSubtreeLikelihoodPostorder(*spTree, *geneTree, 
                                                        geneRoot, seqSp, spID, 
                                                        likelihoodData, lossRates, 
                                                        duplicationRates, speciesIDs, dupData);
# ifdef _OPENMP
		}
	}
	END= omp_get_wtime();
	std::cout << "computeSubtreeLikelihoodPostorder took: "<< END-START <<std::endl;
#endif

	
	/*    VectorTools::print(duplicationRates);
    std::cout << TreeTemplateTools::treeToParenthesis (*spTree, true)<<std::endl;
    std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
    std::cout << "initialLikelihood: "<< initialLikelihood << std::endl;*/
 //  std::cout <<"CLOCKAFTERPOSTORDER "<<clock()<<std::endl;
//std::cout <<"Postorder tree traversal over: Initial likelihood: "<<initialLikelihood<<std::endl;
	//computeSubtreeLikelihoodPrefix then computes the other conditional likelihoods, and also returns the best rooting.
  std::vector <Node *> sons = geneRoot->getSons();

  if (sons.size()!=2) {
    std::cerr <<"Error: "<<sons.size()<< "sons at the root!"<<std::endl; 
  }
  LksToNodes[initialLikelihood]=sons[0];
  //We fill the likelihood and species ID data for the root node.
  //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
  //with other nodes.
  likelihoodData[geneRoot->getId()][1] = likelihoodData[geneRoot->getSon(1)->getId()][0];
  likelihoodData[geneRoot->getId()][2] = likelihoodData[geneRoot->getSon(0)->getId()][0];
  speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon(1)->getId()][0];
  speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon(0)->getId()][0];
  dupData[geneRoot->getId()][1] = dupData[geneRoot->getSon(1)->getId()][0];
  dupData[geneRoot->getId()][2] = dupData[geneRoot->getSon(0)->getId()][0];

    
/*  std::cout <<"likelihoodData[geneRoot->getId()][1]"<<likelihoodData[geneRoot->getId()][1]<<std::endl;
  std::cout <<"likelihoodData[geneRoot->getId()][2]"<<likelihoodData[geneRoot->getId()][2]<<std::endl;
  
  for (int i=0 ; i<geneTree->getNumberOfNodes(); i++) {
    std::cout <<"ID "<<i<<"likelihoodData[ID][0]"<<likelihoodData[i][0]<<std::endl;
  }
  */

# ifdef _OPENMP
	START= omp_get_wtime();
#pragma omp parallel
	{
#pragma omp single
		{
#endif
  for (unsigned int i = 0; i< sons.size(); i++){
    for (unsigned int j =0; j<sons[i]->getNumberOfSons(); j++) {
      computeSubtreeLikelihoodPreorder(*spTree, *geneTree, 
                                       sons[i], seqSp, spID, 
                                       likelihoodData, 
                                       lossRates, duplicationRates, 
                                       speciesIDs, dupData, j, LksToNodes);
    }
  }
# ifdef _OPENMP
		}
	}
	END= omp_get_wtime();
	std::cout << "computeSubtreeLikelihoodPreorder took: "<< END-START <<std::endl;
#endif
	

    
    
  /* 
  std::cout <<"Printing all rooting likelihoods as found by the DR tree traversal"<<std::endl;
  std::map<double, Node*>::iterator it;

  for ( it=LksToNodes.begin() ; it != LksToNodes.end(); it++ )
    std::cout << (*it).second->getId() << " => " << (*it).first << std::endl;
   */
  //  std::cout << "IN findMLReconciliationDR: "<<TreeTemplateTools::treeToParenthesis (*geneTree, false)<<std::endl;

    vector<Node*> nodes = geneTree->getNodes();
    for (unsigned int i = 0 ; i < nodes.size() ; i++ ) {
        if (nodes[i]->hasNodeProperty("outgroupNode") ) {
            nodes[i]->deleteNodeProperty("outgroupNode");
            break;
        }
    }
//    std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, false)<<std::endl;

    
    LksToNodes.rbegin()->second->setNodeProperty("outgroupNode", BppString("here") );
//    std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, false)<<std::endl;

    //TEST02052012
//    geneTree->newOutGroup(LksToNodes.rbegin()->second); //uncomment that if you want to keep gene family trees fixed except for the root
//    std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, false)<<std::endl;

  if (fillTables) {
  
	//Now the best root has been found. I can thus run a function with this best root to fill all the needed tables. This additional tree traversal could be avoided.
	//To this end, the needed tables should be filled by the postfix and prefix traversals. This has not been done yet.
  //resetting
  speciesIDs = std::vector<std::vector<int> > (geneTree->getNumberOfNodes(), nodeSpId);

  dupData = speciesIDs;
  // Getting a well-rooted tree
  TreeTemplate<Node > * tree = geneTree->clone();
  
	  //Test25 05 2012: I think we need to reroot the gene tree with the best root
  tree->newOutGroup(LksToNodes.rbegin()->second->getId());
  
  
 //std::cout << TreeTemplateTools::treeToParenthesis (*tree, true)<<std::endl;

  nodesToTryInNNISearch.clear();

  //Resetting numLineages std::vectors
  resetVector(num0lineages);
  resetVector(num1lineages);
  resetVector(num2lineages);

  
// std::cout <<"HERE_rooted_tree "<<TreeTemplateTools::treeToParenthesis (*tree, true)<<std::endl;
//  std::cout <<"HERE_rooted_tree2 "<<TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
//  std::cout <<"HERE_SP_tree "<<TreeTemplateTools::treeToParenthesis (*spTree, true)<<std::endl;
# ifdef _OPENMP
	  START = omp_get_wtime();
#pragma omp parallel
	  {
#pragma omp single
		  {
#endif
	  computeNumbersOfLineagesFromRoot(spTree, tree, 
                                   tree->getRootNode(), 
                                   seqSp, spID, 
                                   num0lineages, num1lineages, 
                                   num2lineages, speciesIDs, 
                                   dupData, nodesToTryInNNISearch);
# ifdef _OPENMP
		  }
	  }
	  END = omp_get_wtime();
	  std::cout << "computeNumbersOfLineagesFromRoot took: "<< END-START <<std::endl;
#endif

  delete tree;

  }
  
/*
  std::cout <<"num0Lineages :"<<std::endl;
  VectorTools::print(num0lineages);
  std::cout <<"num1Lineages :"<<std::endl;
  VectorTools::print(num1lineages);
  std::cout <<"num2Lineages :"<<std::endl;
  VectorTools::print(num2lineages);
  std::cout <<std::endl;
  */
  
  //We return the best likelihood
  MLindex = LksToNodes.rbegin()->second->getId();
// std::cout <<"Bestlikelihood"<< LksToNodes.rbegin()->first<<std::endl;
    
	return LksToNodes.rbegin()->first;
  
  

}


/**************************************************************************/
//We compute the average loss proportion observed among all species tree branches.
double computeAverageLossProportion(std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages) {
  double totNum0 = (double)VectorTools::sum(num0lineages);
  double totNum1 = (double)VectorTools::sum(num1lineages);
  double totNum2= (double)VectorTools::sum(num2lineages);
  if (totNum0+totNum1+totNum2==0) {//This happens at the beginning of the program
  /*  for (int i =0; i<num0lineages.size(); i++) {
      num0lineages[i]=267;
      num1lineages[i]=723;
      num2lineages[i]=10;
    }*/    
    resetLineageCounts(num0lineages, num1lineages, num2lineages);
  }
  totNum0 = (double)VectorTools::sum(num0lineages);
  totNum1 = (double)VectorTools::sum(num1lineages);
  totNum2= (double)VectorTools::sum(num2lineages);
  double prop = totNum0/(totNum0+totNum1+totNum2);
  return prop;
}


void extractSubVectorsWithInternalLineages(std::vector <int> & num0lineages, 
                                           std::vector <int> & num1lineages, 
                                           std::vector <int> & num2lineages, 
                                           std::map <std::string, int> & genomeMissing, 
                                           TreeTemplate<Node> & tree, 
                                           std::vector <int> & num0, 
                                           std::vector <int> & num1, 
                                           std::vector <int> & num2) {
  std::vector <int> branchesToDiscard = tree.getLeavesId();
  //We also discard the root branch, which by definition does not have any loss
  branchesToDiscard.push_back(0);
  
  sort(branchesToDiscard.begin(), branchesToDiscard.end());
  
  num0 = num0lineages;
  num1 = num1lineages;
  num2 = num2lineages;
  
  
  for(std::vector<int>::reverse_iterator it = branchesToDiscard.rbegin(); it != branchesToDiscard.rend(); it++){
    num0.erase(num0.begin()+*it);
    num1.erase(num1.begin()+*it);
    num2.erase(num2.begin()+*it);
  }
  return;
  
}




//We compute the average loss proportion observed among species tree branches that
//are not found in "genomeMissing" to have missing data.

void extractSubVectorsWithCompletelySequencedLineages(std::vector <int> & num0lineages, 
                                                      std::vector <int> & num1lineages, 
                                                      std::vector <int> & num2lineages, 
                                                      std::map <std::string, int> & genomeMissing, 
                                                      TreeTemplate<Node> & tree, 
                                                      std::vector <int> & num0, 
                                                      std::vector <int> & num1, 
                                                      std::vector <int> & num2) {
  std::vector <int> branchesToDiscard;
  for(std::map<std::string, int >::iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++){
    if (it->second != 0) {
      branchesToDiscard.push_back(tree.getLeafId(it->first));
    }
  }
  //We also discard the root branch, which by definition does not have any loss
  branchesToDiscard.push_back(0);
  
  sort(branchesToDiscard.begin(), branchesToDiscard.end());
  
  num0 = num0lineages;
  num1 = num1lineages;
  num2 = num2lineages;
  
  for(std::vector<int>::reverse_iterator it = branchesToDiscard.rbegin(); it != branchesToDiscard.rend(); it++){
    num0.erase(num0.begin()+*it);
    num1.erase(num1.begin()+*it);
    num2.erase(num2.begin()+*it);
  }
  return;
  
}




double computeAverageLossProportionOnCompletelySequencedLineages(std::vector <int> & num0lineages, 
                                                                 std::vector <int> & num1lineages, 
                                                                 std::vector <int> & num2lineages, 
                                                                 std::map <std::string, int> & genomeMissing, 
                                                                 TreeTemplate<Node> & tree) {
   
  std::vector  <int> num0; 
  std::vector  <int> num1; 
  std::vector  <int> num2; 
  
  //sort (myvector.begin()+4, myvector.end()
  
  extractSubVectorsWithCompletelySequencedLineages(num0lineages, 
                                                   num1lineages, 
                                                   num2lineages, 
                                                   genomeMissing, 
                                                   tree, 
                                                   num0, 
                                                   num1, 
                                                   num2);
  
  double totNum0 = (double)VectorTools::sum(num0);
  double totNum1 = (double)VectorTools::sum(num1);
  double totNum2= (double)VectorTools::sum(num2);
  if (totNum0+totNum1+totNum2==0) {//This happens at the beginning of the program
  /*  for (int i =0; i<num0lineages.size(); i++) {
      num0lineages[i]=267;
      num1lineages[i]=723;
      num2lineages[i]=10;
    }*/
    totNum0 = 267;
    totNum1 = 723;
    totNum2 = 10;
  }
 /* totNum0 = (double)VectorTools::sum(num0lineages);
  totNum1 = (double)VectorTools::sum(num1lineages);
  totNum2= (double)VectorTools::sum(num2lineages);*/
  double prop = totNum0/(totNum0+totNum1+totNum2);
 // std::cout <<"prop : "<<prop<<std::endl;
  return prop;
}


 /**************************************************************************/
//We increase num0lineages in some external branches to account for the low sequence coverage of some genomes. 
//We also set num0Lineages in the root branch.

/*
void alterLineageCountsWithCoverages(std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages, std::map <std::string, int> & genomeMissing, TreeTemplate<Node> & tree) {
  //double propLoss = computeAverageLossProportion(num0lineages, num1lineages, num2lineages);
  double propLoss = computeAverageLossProportionOnCompletelySequencedLineages(num0lineages, num1lineages, num2lineages, genomeMissing, tree);
  
  //At branch 0, by definition, we never count cases where there has been a loss, so we set it to an average value.
  num0lineages[0] = (int)floor( (propLoss)*((double)num1lineages[0]+(double)num2lineages[0]));
  std::cout <<"num0Genes "<< num0lineages[0] <<std::endl;
  
  
  
//  std::cout <<"propLoss "<<propLoss<<std::endl;
  for(std::map<std::string, int >::iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++){
  //   std::cout <<"it first "<<it->first<<std::endl;
    int id = tree.getLeafId(it->first);
    // std::cout <<"After getLeafId 3"<<std::endl;
    int percent = it->second;
//  std::cout <<"it second "<<it->second<<std::endl;
     
    // std::cout <<"obs : "<<obs<<std::endl;
    
    double percentd = (double)percent/100.0;
    double totLoss = propLoss +percentd;
  
    if (totLoss>0.99) {
      totLoss=0.99;
    }
    std::cout <<"percentd "<<percentd<<" propLoss "<<propLoss<<" totLoss "<<totLoss<<std::endl;
    
    while((num1lineages[id] < 10)||(num2lineages[id] < 10)) {
      num1lineages[id] = num1lineages[id]*100;
      num2lineages[id] = num2lineages[id]*100;
      if (num1lineages[id] ==0) {
        num1lineages[id] =1;
      }
      if (num2lineages[id] ==0) {
        num2lineages[id] =1;
      }
    }
    
 //   std::cout <<"Before :"<<num0lineages[id]<< " num1lineages[id] "<<num1lineages[id]<<" num2lineages[id] "<<num2lineages[id]<<std::endl; 
    
    num0lineages[id] = (int)(totLoss * ((double)num1lineages[id]+(double)num2lineages[id]) / (1.0-totLoss));
 //   std::cout <<" after "<<num0lineages[id]<<std::endl;
      
  }
}
*/


/**************************************************************************/
//We increase num0lineages in some external branches to account for the low sequence coverage of some genomes. 
//We also set num0Lineages in the root branch.
void alterLineageCountsWithCoverages(std::vector <int> & num0lineages, 
                                     std::vector <int> & num1lineages, 
                                     std::vector <int> & num2lineages, 
                                     std::map <std::string, int> & genomeMissing, 
                                     TreeTemplate<Node> & tree, bool average) {
  std::vector  <int> num0; 
  std::vector  <int> num1; 
  std::vector  <int> num2; 
  extractSubVectorsWithCompletelySequencedLineages(num0lineages, 
                                                   num1lineages, 
                                                   num2lineages, 
                                                   genomeMissing, 
                                                   tree, 
                                                   num0, num1, num2);
  
  double avg0d = VectorTools::mean<int, double>(num0);
  double avg1d = VectorTools::mean<int, double>(num1);
  double avg2d = VectorTools::mean<int, double>(num2);

  if (avg0d+avg1d+avg2d==0) {//This shouldn't happen but if it does...
    std::cerr<<"WARNING: No event on average on all branches. Setting default values."<<std::endl;
    avg0d = 267;
    avg1d = 723;
    avg2d = 10;
  }

  double propLoss = avg0d/(avg0d+avg1d+avg2d);
  if (propLoss > 0.99) {
    propLoss = 0.99;
  }
  
  int avg0 = (int)avg0d;
  int avg1 = (int)avg1d;
  int avg2 = (int)avg2d;
    
  if (average) {
    for (unsigned int i =0; i<num0lineages.size(); i++) {
      num0lineages[i]=avg0;
      num1lineages[i]=avg1;
      num2lineages[i]=avg2;
    }
  }
  else {
    //At branch 0 and branches surrounding the root, by definition, we never count cases where there has been a loss, so we set it to an average value.
    //We also put all num0, num1, and num2 values at the average values computed above.
    for (unsigned int i =0; i<num0lineages.size(); i++) {
      if (num0lineages[i] == 0) {
        num0lineages[i]=avg0;
        num1lineages[i]=avg1;
        num2lineages[i]=avg2;
      }
    }
  }
    
  //  Now we apply corrections for poorly sequenced genomes
  for(std::map<std::string, int >::iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++){
    int id = tree.getLeafId(it->first);
    int percent = it->second;
    if (percent>0) {
      double percentd = (double)percent/100.0;
      double totLoss = propLoss +percentd;     
      if (totLoss>0.99) {
        totLoss=0.99;
      }
      
      while((num1lineages[id] < 10)||(num2lineages[id] < 10)) {
        num1lineages[id] = num1lineages[id]*100;
        num2lineages[id] = num2lineages[id]*100;
        if (num1lineages[id] < 1) {
          num1lineages[id] =1;
        }
        if (num2lineages[id] < 1) {
          num2lineages[id] =1;
        }
      }
      num0lineages[id] = (int)(totLoss * ((double)num1lineages[id]+(double)num2lineages[id]) / (1.0-totLoss));
    }
  }
}

/**************************************************************************/
//We average all count values, and increase num0lineages in some external branches to account for the low sequence coverage of some genomes. 
/*
void alterLineageCountsWithCoveragesAverage(std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages, std::map <std::string, int> & genomeMissing, TreeTemplate<Node> & tree) {
  //Averaging the std::vectors

  std::vector  <int> num0; 
  std::vector  <int> num1; 
  std::vector  <int> num2; 
  
  extractSubVectorsWithCompletelySequencedLineages(num0lineages, num1lineages, num2lineages, genomeMissing, tree, num0, num1, num2);
  
  
  double avg0d = VectorTools::mean<int, double>(num0);
  double avg1d = VectorTools::mean<int, double>(num1);
  double avg2d = VectorTools::mean<int, double>(num2);
  

//  std::cout <<"BEFORE "<<avg0d<<" "<<avg1d<<" "<<avg2d<<std::endl;
  while((avg0d < 10.0)||(avg1d < 10.0)||(avg2d < 10.0 )) {
    avg0d = avg0d*100;
    avg1d = avg1d*100;
    avg2d = avg2d*100;
    if (avg0d ==0.0) {
      avg0d =1;
    }
    if (avg1d ==0.0) {
      avg1d =1;
    }
    if (avg2d ==0.0) {
      avg2d =1;
    }
  }
  int avg0 = (int)avg0d;
  int avg1 = (int)avg1d;
  int avg2 = (int)avg2d;
  
//  std::cout <<"avg0 : "<<avg0<<" avg1 : "<<avg1<< "avg2: "<<avg2<<std::endl;
  
  for (int i =0; i<num0lineages.size(); i++) {
    num0lineages[i]=avg0;
    num1lineages[i]=avg1;
    num2lineages[i]=avg2;
  }
  alterLineageCountsWithCoverages(num0lineages, num1lineages, num2lineages, genomeMissing, tree);


}
*/
/**************************************************************************/
//We set all duplication and loss values to values obtained on a former dataset (loss rate 0.313, duplication rate 0.0159, which corresponds to num0lineages=267, num1lineages=723, num2lineages=10), and increase num0lineages in some external branches to account for the low sequence coverage of some genomes. 
void alterLineageCountsWithCoveragesInitially(std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages, std::map <std::string, int> & genomeMissing, TreeTemplate<Node> & tree) {
  //Setting initial values for the std::vectors
  resetLineageCounts(num0lineages, num1lineages, num2lineages);
  
  alterLineageCountsWithCoverages(num0lineages, num1lineages, num2lineages, genomeMissing, tree, false);
}


void resetLineageCounts(std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages) {  
  for (unsigned int i =0; i<num0lineages.size(); i++) {
    num0lineages[i]=267;
    num1lineages[i]=723;
    num2lineages[i]=1;
   /* num0lineages[i]=267*(i+1);
    num1lineages[i]=723;
    num2lineages[i]=10;*/
  }
}

/**************************************************************************/
// This function removes properties associated to a tree.



void deleteSubtreeProperties(Node &node) {
 for(unsigned int i = 0; i < node.getNumberOfSons(); i++)
   {
     Node * son = node.getSon(i);
     deleteSubtreeProperties(* son);
     son->deleteBranchProperties();
     son->deleteNodeProperties();
   }
}



void deleteTreeProperties(TreeTemplate<Node> & tree) {
  deleteSubtreeProperties(*(tree.getRootNode()));
}




/**************************************************************************/


std::map <std::string, int> computeSpeciesNamesToIdsMap (TreeTemplate<Node> & tree) {
  std::map <std::string, int> spId;
  std::vector <Node *> nodes = tree.getNodes();
  for (unsigned int i = 0; i< nodes.size() ; i++) {
    if (nodes[i]->isLeaf()) {
      spId[nodes[i]->getName()] = nodes[i]->getId(); 
    }
  }
  return spId;
}

/**************************************************************************/



void computeDuplicationAndLossRatesForTheSpeciesTree (std::string &branchProbaOptimization, 
                                                      std::vector <int> & num0Lineages, 
                                                      std::vector <int> & num1Lineages, 
                                                      std::vector <int> & num2Lineages, 
                                                      std::vector<double> & lossExpectedNumbers, 
                                                      std::vector<double> & duplicationExpectedNumbers, 
                                                      std::map <std::string, int> & genomeMissing, 
                                                      TreeTemplate<Node> & tree) {
  //outputting trees with branch lengths in numbers of events, before correction.
  computeDuplicationAndLossProbabilitiesForAllBranches (num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers);
  std::cout << "Species tree with expected numbers of duplications as branch lengths:"<<std::endl;
  for (unsigned int i =0; i<num0Lineages.size() ; i++ ) 
    {
    tree.getNode(i)->setBranchProperty("DUPLICATIONS", Number<double>( duplicationExpectedNumbers[i]));
    if (tree.getNode(i)->hasFather()) 
      {
      tree.getNode(i)->setDistanceToFather(duplicationExpectedNumbers[i]);
      }
    }
  std::cout << treeToParenthesisWithDoubleNodeValues(tree, false, "DUPLICATIONS")<<std::endl;
  std::cout << "Species tree with expected numbers of losses as branch lengths:"<<std::endl;
  for (unsigned int i =0; i<num0Lineages.size() ; i++ ) 
    {
    tree.getNode(i)->setBranchProperty("LOSSES", Number<double>(lossExpectedNumbers[i]));
    if (tree.getNode(i)->hasFather()) 
      {
      tree.getNode(i)->setDistanceToFather(lossExpectedNumbers[i]);
      }
    }
  std::cout << treeToParenthesisWithDoubleNodeValues(tree, false, "LOSSES")<<std::endl;
  
  //Doing the correction:
  if (branchProbaOptimization=="average") {
	  std::cout << "Averaging branch-wise estimates."<<std::endl;
    alterLineageCountsWithCoverages(num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree, true);
  }
  else {
    alterLineageCountsWithCoverages(num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree, false);
  }
  
  computeDuplicationAndLossProbabilitiesForAllBranches (num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers);

}


/**************************************************************************/
void computeDuplicationAndLossRatesForTheSpeciesTreeInitially (std::string &branchProbaOptimization, 
                                                               std::vector <int> & num0Lineages, 
                                                               std::vector <int> & num1Lineages, 
                                                               std::vector <int> & num2Lineages, 
                                                               std::vector<double> & lossProbabilities, 
                                                               std::vector<double> & duplicationProbabilities, 
                                                               std::map <std::string, int> & genomeMissing, 
                                                               TreeTemplate<Node> & tree) {
  std::cout <<"Computing Initial Rates"<<std::endl; 
  alterLineageCountsWithCoveragesInitially(num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree);
  computeDuplicationAndLossProbabilitiesForAllBranches (num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities);
}


/**************************************************************************/
/*Removes a leaf in a tree*/
void removeLeaf(TreeTemplate<Node> & tree, std::string toRemove){
 // std::cout <<"in removeLeaf 1"<<std::endl;
  Node * NToRemove  = tree.getNode(toRemove);
  if (!NToRemove->hasFather()) {
   // std::cout <<"node is root !!!!!"<<std::endl; 
    Node * father = NToRemove->getSon(0);
    Node * son1 = father->getSon(1);
    father->removeFather();
    tree.newOutGroup(son1->getId());
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
 // std::cout <<"in removeLeaf 2"<<std::endl;
  Node * father = NToRemove->getFather();
 // std::cout <<"in removeLeaf 3"<<std::endl;
  Node * brother;
  for(unsigned int i=0;i<father->getNumberOfSons();i++)
    if(father->getSon(i)!=NToRemove){brother=father->getSon(i); break;}
 // std::cout <<"in removeLeaf 4"<<std::endl;
  double distBro;
  try{
    distBro = brother->getDistanceToFather();
  } catch (std::exception) {
    distBro = 0.000001;
  }
 // std::cout <<"in removeLeaf 5"<<std::endl;
  if (!father->hasFather()) {
   // std::cout <<"father is root !!!!!"<<std::endl; 
   // brother->removeFather();
    tree.rootAt(brother->getId());
   // std::cout <<"After rootAt"<<std::endl;
    //tree.newOutGroup(brother->getSon(0)->getId());
    for (unsigned int i = 0; i<brother->getNumberOfSons(); i++) {
      if ( brother->getSon(i)==NToRemove) {
        brother->removeSon(i);
        break;
      }
    }
    tree.newOutGroup(brother->getSon(0)->getId());
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
  double distFa;
  try{
    distFa = father->getDistanceToFather();
  } catch (std::exception) {
    distFa = 0.000001;
  }
 // std::cout <<"in removeLeaf 6"<<std::endl;
  Node * grandFather = father->getFather();
  grandFather->addSon(brother);
  brother->setDistanceToFather(distBro+distFa);
  for(unsigned int i=0;i<grandFather->getNumberOfSons();i++)
    if(grandFather->getSon(i)==father){grandFather->removeSon(i); break;}
  /*if (!grandFather->hasFather()) {
    tree.newOutGroup(grandFather->getSon(0)->getId());
  }*/
  tree.resetNodesId();
/*  delete NToRemove;
  delete father;*/
}

/**************************************************************************/

/******************************************************************************/
//This function cleans a vector of options read from a file.
/******************************************************************************/


void cleanVectorOfOptions (std::vector<std::string> & listOptions, bool sizeConstraint)
{
  int maxStrSize=0;
  std::vector <int> toRemove;
  int i=0;
  for(std::vector<std::string >::iterator it = listOptions.begin(); it != listOptions.end(); it++)
    {
    std::string arg = removeComments(*it, std::string("#"), std::string("\n"));//remove shell comments.
    arg = removeComments(arg, std::string("//"), std::string("\n"));//remove C simple comments.
    arg = removeComments(arg, std::string("/*"), std::string("*/"));//remove C multiple comments.
    arg = TextTools::removeWhiteSpaces(arg);
    *it = arg;
    int temp=it->length();
    if (temp==0) {
      toRemove.push_back(i);
    }
    else {
      if (temp>maxStrSize){maxStrSize=temp;};
    }
    i++;
    }
  for(std::vector<int>::reverse_iterator it = toRemove.rbegin(); it != toRemove.rend(); it++)
    {
    listOptions.erase(listOptions.begin()+*it);
    }
  if (sizeConstraint) 
    {
    if(maxStrSize>MAXFILENAMESIZE) 
      {
      std::cout << "\nNames are too long, please abbreviate ! File names (including the path) need to be < "<< MAXFILENAMESIZE <<" letters long.\n"<<std::endl;
          MPI::COMM_WORLD.Abort(1);
      exit(-1);
      }
    }
  return;
}

/**************************************************************************/


/**************************************************************************
 * This function distributes gene families so as to homogenize computational load between clients.
 **************************************************************************/

bool sortMaxFunction (std::pair <std::string, double> i, std::pair <std::string, double> j) 
{ 
  if (i.second > j.second ) 
    {
    return (true);
    }
  else 
    {
    return (false);
    }
}

bool sortMinFunction (std::pair <std::vector<std::string>, double> i, std::pair <std::vector<std::string>, double> j) 
{ 
  if (i.second < j.second ) 
    {
    return (true);
    }
  else 
    {
    return (false);
    }
}




void generateListOfOptionsPerClient(std::vector <std::string> listOptions, int size, std::vector <std::vector<std::string> > &listOfOptionsPerClient, std::vector <unsigned int> &numbersOfGenesPerClient) 
{
  //Here, two alternatives: either we do have information regarding the gene family sizes, or we don't.
  //Is there size information? == Is there a ":" in the first line?
  if (TextTools::hasSubstring(listOptions[0],":")) 
    {
    std::vector <std::pair <std::string, double> > elements;
    for (unsigned int i = 0; i<listOptions.size() ; i++) {
      StringTokenizer st1(listOptions[i], ":", true);
      elements.push_back(std::pair <std::string, double>(st1.getToken(0), TextTools::toDouble(st1.getToken(1))) );
    }
    
    //Now sort the gene families by their size, in descending order
    sort(elements.begin(), elements.end(), sortMaxFunction);
    //Now we assign gene families to nodes.
    //We start with big families, and then go in decreasing order.
    //First we assign the first families
    std::vector <std::pair <std::vector<std::string>, double> > listOfOptionsAndTotSizePerClient;
    
    std::vector<std::string> temp2;
    std::pair <std::vector<std::string>, double> temp3;
    unsigned int j = 0;
    
    for (int i = 0; i<size-1 ; i++) {
      listOfOptionsAndTotSizePerClient.push_back(temp3);
      listOfOptionsAndTotSizePerClient[i].first.push_back(elements[j].first);
      listOfOptionsAndTotSizePerClient[i].second = elements[j].second;  
      j = j+1;
    }
    
    while (j<listOptions.size()) {
      //We sort listOfOptionsAndTotSizePerClient in increasing function
      sort(listOfOptionsAndTotSizePerClient.begin(), listOfOptionsAndTotSizePerClient.end(), sortMinFunction);
      listOfOptionsAndTotSizePerClient[0].first.push_back(elements[j].first);
      listOfOptionsAndTotSizePerClient[0].second = listOfOptionsAndTotSizePerClient[0].second + elements[j].second;
      j= j+1;
    }
    
    //Now all gene families must have been assigned to nodes.
    // std::vector <std::vector<std::string> > listOfOptionsPerClient;
    listOfOptionsPerClient.push_back(temp2);
    listOfOptionsPerClient[0].push_back(std::string("####"));//For the root node
    numbersOfGenesPerClient.push_back(0);
    
    //We print the result of the assignment and fill listOfOptionsPerClient:
    for ( int i = 0; i<size-1 ; i++) {
      std::cout <<"Client "<<i<<" is in charge of "<< listOfOptionsAndTotSizePerClient[i].first.size()<<" gene families; Total Weight : "<<listOfOptionsAndTotSizePerClient[i].second<<std::endl;
      listOfOptionsPerClient.push_back(listOfOptionsAndTotSizePerClient[i].first);
      numbersOfGenesPerClient.push_back(listOfOptionsAndTotSizePerClient[i].first.size());
    }
    return ;
    //generateListOfOptionsPerClient(listOptions, size, listOfOptionsPerClient, numbersOfGenesPerClient);
    }
  else 
    {
    int numberOfGenesPerClient = (int)(listOptions.size()) / (size -1); 
    int reste = (listOptions.size()) % (size -1);
    std::cout <<"Number of genes per client : "<<numberOfGenesPerClient<< " and extra "<<reste<<std::endl;
    for (int i = 0 ; i< size-1 ; i++ ) 
      {
      if (i<reste) 
        {
        numbersOfGenesPerClient.push_back(numberOfGenesPerClient+1);
        }
      else 
        {
        numbersOfGenesPerClient.push_back(numberOfGenesPerClient);
        }
      } 
    std::vector<unsigned int>::iterator it2 = numbersOfGenesPerClient.begin();
    numbersOfGenesPerClient.insert( it2, int(0) ); //For the server, we insert a "dumb" option file at the beginning of the std::vector, so only clients compute the reconciliation
    int currentFile = 0;
    std::vector<std::string> temp2;
    for (int i = 0 ; i< size ; i++ ) 
      {
      listOfOptionsPerClient.push_back(temp2);
      if (i == 0) 
        {
        listOfOptionsPerClient[i].push_back(std::string("####"));
        }
      else 
        {
        for (unsigned int j = 0 ; j < numbersOfGenesPerClient[i] ; j++) 
          {
          listOfOptionsPerClient[i].push_back(listOptions[currentFile]);
          currentFile++;
          }
        }
      }
    return;
    }
  
}





/**************************************************************************
 * Utilitary fonction
 *************************************************************************/

std::string removeComments(
                           const std::string & s,
                           const std::string & begin,
                           const std::string & end)
{
  std::string r = s;
  std::string::size_type last = 0;
	do
    {
    std::string::size_type first = r.find(begin, last);
		if(first == std::string::npos) return r; //No shell comment.
		//else:  
		last = r.find(end, first);
		if(last == std::string::npos)
      {
			r.erase(r.begin() + first, r.end());
      }
		else
      {
			r.erase(r.begin() + first, r.begin() + last);
      }
    } while(last != std::string::npos);
	return r;
}




/**************************************************************************
 * Annotate gene tree with duplication events. Adds D=Y to node properties where
 * there was a duplication, and D=N where there was no duplication.
 * Also adds species node id for each node of the gene tree, using :S=spId.
 *************************************************************************/

void annotateGeneTreeWithDuplicationEvents (TreeTemplate<Node> & spTree, 
                                            TreeTemplate<Node> & geneTree, 
                                            Node * node, 
                                            std::map<std::string, std::string > seqSp,
                                            std::map<std::string, int > spID) 
{
 	if (node->isLeaf()) {
        node->setNodeProperty("S", BppString(TextTools::toString(assignSpeciesIdToLeaf(node, seqSp, spID))));
        node->setBranchProperty("Ev", BppString("S"));
        return;
    }
    else {
        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            annotateGeneTreeWithDuplicationEvents(spTree, geneTree, sons[i], seqSp, spID);
        }
        
        int a = TextTools::toInt((dynamic_cast<const BppString*>(sons[0]->getNodeProperty("S")))->toSTL());
        int b = TextTools::toInt((dynamic_cast<const BppString*>(sons[1]->getNodeProperty("S")))->toSTL());
        
        int aold = a;
        int bold = b;
        int lossA = 0;
        int lossB = 0;
        
        while (a!=b) 
        {
            if (a>b){
                a = spTree.getNode(a)->getFather()->getId(); 
                lossA = lossA +1;
            }
            else {
                b = spTree.getNode(b)->getFather()->getId(); 
                lossB = lossB + 1;
            }
        }
        sons[0]->setBranchProperty("L", BppString(TextTools::toString(lossA)));
        sons[1]->setBranchProperty("L", BppString(TextTools::toString(lossB)));
        node->setNodeProperty("S", BppString(TextTools::toString(a)));
        if ((a == aold ) || (a == bold))
        {
            node->setBranchProperty("Ev", BppString("D"));
        }
        else 
        {
            node->setBranchProperty("Ev", BppString("S"));
        }
        
        return;
	}
}


/**************************************************************************
 * Annotate gene tree with duplication events. Adds D=Y to node properties where
 * there was a duplication, and D=N where there was no duplication.
 * Also adds species node id for each node of the gene tree, using :S=spId.
 *************************************************************************/

void annotateGeneTreeWithScoredDuplicationEvents (TreeTemplate<Node> & spTree, 
                                            TreeTemplate<Node> & geneTree, 
                                            Node * node, 
                                            std::map<std::string, std::string > seqSp,
                                            std::map<std::string, int > spID) 
{
  //    std::cout << "annotateGeneTreeWithScoredDuplicationEvents "<<std::endl;
 	if (node->isLeaf()) {
	  //	   std::cout << "annotateGeneTreeWithScoredDuplicationEvents Leaf"<<std::endl;

		int sp = assignSpeciesIdToLeaf(node, seqSp, spID);
        node->setNodeProperty("S", BppString(TextTools::toString(sp)));
		node->setNodeProperty("spPresentInSubtree", BppVector<int> (1, sp));
        node->setBranchProperty("Ev", BppString("S"));
	//	std::cout << "annotateGeneTreeWithScoredDuplicationEvents Leaf 2"<<std::endl;

        return;
    }
    else {
      //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf "<<std::endl;

        std::vector <Node *> sons = node->getSons();
        for (unsigned int i = 0; i< sons.size(); i++){
            annotateGeneTreeWithScoredDuplicationEvents(spTree, geneTree, sons[i], seqSp, spID);
        }
        
	//std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 2"<<std::endl;

        int a = TextTools::toInt((dynamic_cast<const BppString*>(sons[0]->getNodeProperty("S")))->toSTL());
        int b = TextTools::toInt((dynamic_cast<const BppString*>(sons[1]->getNodeProperty("S")))->toSTL());
	//	std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 3"<<std::endl;

        int aold = a;
        int bold = b;
        int lossA = 0;
        int lossB = 0;
        
        while (a!=b) 
        {
            if (a>b){
                a = spTree.getNode(a)->getFather()->getId(); 
                lossA = lossA +1;
            }
            else {
                b = spTree.getNode(b)->getFather()->getId(); 
                lossB = lossB + 1;
            }
        }
	//	std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 4"<<std::endl;

        sons[0]->setBranchProperty("L", BppString(TextTools::toString(lossA)));
        sons[1]->setBranchProperty("L", BppString(TextTools::toString(lossB)));
        node->setNodeProperty("S", BppString(TextTools::toString(a)));
		//VectorTools::print( (dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL() );
		//VectorTools::print( (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL() );

		std::vector<int> vec = VectorTools::vectorUnion((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
		node->setNodeProperty("spPresentInSubtree", BppVector<int> (vec.begin(), vec.end() ));
		//		std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 5"<<std::endl;

        if ((a == aold ) || (a == bold))
        {
            node->setBranchProperty("Ev", BppString("D"));
			node->setNodeProperty("Ev", BppString("D"));
			//Need to compute the intersection (and its size) and divide by the size of the union (computed above)
			//to compute the score of the duplication event.
			std::vector<int> vec2 = VectorTools::vectorIntersection((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
			node->setNodeProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() )));
			node->setBranchProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() ))); //TEMPORARY, JUST FOR DEBUGGING
        }
        else 
        {
            node->setBranchProperty("Ev", BppString("S"));
        }
	//	std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 6"<<std::endl;

        return;
	}
}

/*
void editDuplicationNodesMuffato2(TreeTemplate<Node> & spTree, 
								 TreeTemplate<Node> & geneTree,
								 Node * node,
								 double editionThreshold) {
 	if (node->isLeaf()) { }
	else {
		if (node->hasNodeProperty("Score")) {
			double dist = 0.1;
			double score = TextTools::toDouble((dynamic_cast<const BppString*>(node->getNodeProperty("Score")))->toSTL());
			std::cout << "editDuplicationNodesMuffato 2: score: "<< score  <<std::endl;					std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(geneTree, true) << std::endl;					std::cout << "node id: "<< TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL()) <<std::endl;					std::cout << "son 0 id: "<< TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty("S")))->toSTL()) <<std::eendl; std::cout << "son 1 id: "<< TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty("S")))->toSTL()) <<std::endl; std::cout << "editionThreshold "<< editionThreshold <<std::endl;
					if (score <= editionThreshold && ( ( TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL()) != TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty("S")))->toSTL() ) ) || (TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL()) != TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty("S")))->toSTL() ) ) ) ) {
				//we edit the node
				int ancestorSpecies = TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL());
				Node *ancestor = spTree.getNode(ancestorSpecies);
				int son0Id = ancestor->getSon(0)->getId();
				int son1Id = ancestor->getSon(1)->getId();	
				Node *son0;
				Node *son1;
				if (! node->getSon(0)->isLeaf() ) {
					son0= node->getSon(0);
				}
				else {
					//Need to create a father for node son0
					//And need to rename this father node0
					Node* grandSon0 = node->getSon(0);
					unsigned int temp = 0;
					node->removeSon(temp);
					son0 = new Node();
					son0->addSon( grandSon0 );
					node->addSon( son0 );
					son0->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					grandSon0->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
				}
				if (! node->getSon(1)->isLeaf() ) {	
					son1 = node->getSon(1);
				}
				else {
					//Need to create a father for node son1
					//And need to rename this father node1
					Node* grandSon1 = node->getSon(1);
					unsigned int temp = 1;
					node->removeSon(temp);
					son1 = new Node();
					son1->addSon( grandSon1 );
					node->addSon( son1 );
					son1->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					grandSon1->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
				}
				//now we need to move the grandsons around to make sure that everything is in the right place
				//First, we get the node ids in the subtrees from the species tree
				std::vector<int> son0UnderlyingIds = TreeTemplateTools::getNodesId( *(ancestor->getSon(0) ) );
				//std::vector<int> son1UnderlyingIds = ancestor->getSon(1)->getNodesId();
				//Now, we move, in the gene tree, the grandsons to make sure that everything is in the right place
				std::map<Node*, int> sonsOrGrandSons;
				son0->setNodeProperty("S", BppString(TextTools::toString(son0Id)));
				for (unsigned int i = 0 ; i < son0->getNumberOfSons() ; i++) {
					sonsOrGrandSons[son0->getSon(i)] = TextTools::toInt((dynamic_cast<const BppString*>(son0->getSon(i)->getNodeProperty("S")))->toSTL());
				}
				son1->setNodeProperty("S", BppString(TextTools::toString(son1Id)));
				for (unsigned int i = 0 ; i < son1->getNumberOfSons() ; i++) {
					sonsOrGrandSons[son1->getSon(i)] = TextTools::toInt((dynamic_cast<const BppString*>(son1->getSon(i)->getNodeProperty("S")))->toSTL());
				}
				for (std::map<Node*, int >::iterator it = sonsOrGrandSons.begin(); it != sonsOrGrandSons.end(); it++){
					if (VectorTools::contains(son0UnderlyingIds, it->second) ) {
					
						std::vector<Node*> ns = son0->getSons();
						if (! VectorTools::contains (ns, it->first) ) {
						  makeSPR(geneTree, it->first->getId(), son0->getSon(0)->getId(), false, false );
						  //							son1->removeSon( it->first );
						 //	son0->addSon( it->first );
						//	it->first->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
						  
						}
					}
					else {
						std::vector<Node*> ns = son1->getSons();
						if (! VectorTools::contains (ns, it->first) ) {
						  makeSPR(geneTree, it->first->getId(), son1->getSon(0)->getId(), false, false );
						  // son0->removeSon( it->first );
						//	son1->addSon( it->first );
						//	it->first->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
						  
						}
					}
				}
				//Removing nodes of out-degree 1
				if (son0->getNumberOfSons() == 1) {
					Node* father =  son0->getFather();
					Node* son = son0->getSon(0);
					son0->removeSon(son);
					father->removeSon(son0);
					father->addSon(son);
					son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					delete son0;
				}
				if (son1->getNumberOfSons() == 1) {
					Node* father =  son1->getFather();
					Node* son = son1->getSon(0);
					son1->removeSon(son);
					father->removeSon(son1);
					father->addSon(son);
					son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					delete son1;
				}
				//Now we need to update the scores of the sons
				if (!node->getSon(0)->isLeaf()) {
				std::vector <Node *> sons = node->getSon(0)->getSons();
				std::vector<int> vec = VectorTools::vectorUnion((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
				std::vector<int> vec2 = VectorTools::vectorIntersection((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					node->getSon(0)->setNodeProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() )));
				}
				if (!node->getSon(1)->isLeaf()) {
					std::vector <Node *> sons = node->getSon(1)->getSons();
					std::vector<int> vec = VectorTools::vectorUnion((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					std::vector<int> vec2 = VectorTools::vectorIntersection((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					node->getSon(1)->setNodeProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() )));
				}
				std::cout <<"Edited tree from ReconciliationTools: \n"<< TreeTemplateTools::treeToParenthesis(geneTree, false, "Score") << "\n" << std::endl;
			}
        }
		//				std::cout << "editDuplicationNodesMuffato 9"  <<std::endl;		std::cout << TreeTemplateTools::treeToParenthesis(spTree, true)<<std::endl;		std::cout << "editDuplicationNodesMuffato 10"  <<std::endl;		std::cout << TreeTemplateTools::treeToParenthesis(geneTree, true)<<std::endl;		std::cout << "editDuplicationNodesMuffato 11"  <<std::endl;		std::cout << "editDuplicationNodesMuffato: id"<< node->getSon(0)->getId()  <<std::endl;		std::cout << "editDuplicationNodesMuffato: id"<< node->getSon(1)->getId()  <<std::endl;		
		//continue editing
		editDuplicationNodesMuffato(spTree, 
									geneTree,
									node->getSon(0), 
									editionThreshold);
		editDuplicationNodesMuffato(spTree, 
									geneTree,
									node->getSon(1), 
									editionThreshold);
	}
	//	std::cout <<"Edited tree from ReconciliationTools: \n"<< TreeTemplateTools::treeToParenthesis(geneTree, false, "Score") << "\n" << std::endl;
}
*/

/* Better version of the Muffato algorithm: we can edit deeper in the tree, 
 for instance to get one rogue sequence lost in a large tree.*/

/*
void editDuplicationNodesMuffato(TreeTemplate<Node> & spTree, 
								 TreeTemplate<Node> & geneTree,
								 Node * node,
								 double editionThreshold) {
 	if (node->isLeaf()) { }
	else {
		//		std::cout << "editDuplicationNodesMuffato "  <<std::endl;
		
		if (node->hasNodeProperty("Score")) {
			double dist = 0.1;
			double score = TextTools::toDouble((dynamic_cast<const BppString*>(node->getNodeProperty("Score")))->toSTL());
			//		std::cout << "editDuplicationNodesMuffato 2: score: "<< score  <<std::endl;
			//		std::cout <<"Gene Tree: \n" <<    TreeTemplateTools::treeToParenthesis(geneTree, true) << std::endl;
			
			if (score <= editionThreshold && TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL()) != TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(0)->getNodeProperty("S")))->toSTL()) && TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL()) != TextTools::toInt((dynamic_cast<const BppString*>(node->getSon(1)->getNodeProperty("S")))->toSTL()) ) {
				//we edit the node
				int ancestorSpecies = TextTools::toInt((dynamic_cast<const BppString*>(node->getNodeProperty("S")))->toSTL());
				//		std::cout << "editDuplicationNodesMuffato 2: spid: "<< ancestorSpecies <<" and id: "<< node->getId() <<std::endl;
				
				Node *ancestor = spTree.getNode(ancestorSpecies);
				int son0Id = ancestor->getSon(0)->getId();
			//	int son1Id = ancestor->getSon(1)->getId();	
				Node *son0;
				Node *son1;
				//				std::cout << "editDuplicationNodesMuffato 3"  <<std::endl;
				
				if (! node->getSon(0)->isLeaf() ) {
					son0= node->getSon(0);
				}
				else {
					//Need to create a father for node son0
					//And need to rename this father node0
					Node* grandSon0 = node->getSon(0);
					unsigned int temp = 0;
					node->removeSon(temp);
					son0 = new Node();
					son0->addSon( grandSon0 );
					node->addSon( son0 );
					son0->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					grandSon0->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
				}
				//				std::cout << "editDuplicationNodesMuffato 4"  <<std::endl;
				
				if (! node->getSon(1)->isLeaf() ) {	
					son1 = node->getSon(1);
				}
				else {
					//Need to create a father for node son1
					//And need to rename this father node1
					Node* grandSon1 = node->getSon(1);
					unsigned int temp = 1;
					node->removeSon(temp);
					son1 = new Node();
					son1->addSon( grandSon1 );
					node->addSon( son1 );
					son1->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					grandSon1->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
					
				}
				//				std::cout << "editDuplicationNodesMuffato 5"  <<std::endl;
				
				//now we need to move the grandsons around to make sure that everything is in the right place
				//First, we get the node ids in the subtrees from the species tree
				std::vector<int> son0UnderlyingIds = TreeTemplateTools::getNodesId( *(ancestor->getSon(0) ) );
				//std::vector<int> son1UnderlyingIds = ancestor->getSon(1)->getNodesId();
				//Now, we move, in the gene tree, the grandsons to make sure that everything is in the right place
				std::map<Node*, int> sonsOrGrandSonsOrBeyond;
				son0->setNodeProperty("S", BppString(TextTools::toString(son0Id)));
				std::vector<Node*> nodesFromSon0 = TreeTemplateTools::getNodes( *son0  );
				unsigned int siz = nodesFromSon0.size();
				for (unsigned int i = 0 ; i < siz ; ++i) {
					sonsOrGrandSonsOrBeyond[ nodesFromSon0[i] ] =  TextTools::toInt((dynamic_cast<const BppString*>(nodesFromSon0[i]->getNodeProperty("S")))->toSTL());
				}
				std::vector< Node*> nodesFromSon1 = TreeTemplateTools::getNodes( *son1  );
				siz = nodesFromSon1.size();
				for (unsigned int i = 0 ; i < siz ; ++i) {
					sonsOrGrandSonsOrBeyond[ nodesFromSon1[i] ] =  TextTools::toInt((dynamic_cast<const BppString*>(nodesFromSon1[i]->getNodeProperty("S")))->toSTL());
				}
				//Now we have obtained the species ids for all nodes in the two subtrees descending from node.
				//We put the nodes in the right subtree: nodes that should descend from species 0 will go there, 
				//and nodes that should descent from species 1 will go there.
				std::map<Node*, int >::iterator end = sonsOrGrandSonsOrBeyond.end();
				std::set< Node* > blackListed ;
				for (std::map<Node*, int >::iterator it = sonsOrGrandSonsOrBeyond.begin(); it != end; ++it){
					if ( blackListed.find(it->first) == blackListed.end() ) { // if the node has not been blacklisted yet
						if (VectorTools::contains(son0UnderlyingIds, it->second) ) {
							std::vector<Node*> ns = son0->getSons();
							if (! VectorTools::contains (ns, it->first) ) {
								it->first->getFather()->removeSon( it->first );
								son0->addSon( it->first );
								it->first->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
								std::vector < Node* > nodesToBlackList = TreeTemplateTools::getNodes( *(it->first) );
								unsigned int siz2 = nodesToBlackList.size();
								for (unsigned int j = 0 ; j < siz2 ; j++ ) {
									blackListed.insert( nodesToBlackList[j] );
								}
							}
						}
						else {
							std::vector<Node*> ns = son1->getSons();
							if (! VectorTools::contains (ns, it->first) ) {
								it->first->getFather()->removeSon( it->first );
								son1->addSon( it->first );
								it->first->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
							}
						}
					}
				}

				//Removing nodes of out-degree 1, and putting nodes of out-degree > 2 to 2. (TODO)
				nodesFromSon0 = TreeTemplateTools::getNodes( *son0  );
				siz = nodesFromSon0.size();
				for (unsigned int i = 0 ; i < siz ; ++i) {
					if (nodesFromSon0[i]->getNumberOfSons() == 1) {
						Node* father =  nodesFromSon0[i]->getFather();
						Node* son = nodesFromSon0[i]->getSon(0);
						nodesFromSon0[i]->removeSon(son);
						father->removeSon(nodesFromSon0[i]);
						father->addSon(son);
						son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
						delete nodesFromSon0[i];
					}
					if (nodesFromSon0[i]->getNumberOfSons() > 2) {
						Node* father =  nodesFromSon0[i]->getFather();
						Node* son = nodesFromSon0[i]->getSon(0);
						nodesFromSon0[i]->removeSon(son);
						father->removeSon(nodesFromSon0[i]);
						father->addSon(son);
						son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
						delete nodesFromSon0[i];
					}

				}
				nodesFromSon1 = TreeTemplateTools::getNodes( *son1  );
				siz = nodesFromSon1.size();
				for (unsigned int i = 0 ; i < siz ; ++i) {
					if (nodesFromSon1[i]->getNumberOfSons() == 1) {
						Node* father =  nodesFromSon1[i]->getFather();
						Node* son = nodesFromSon1[i]->getSon(0);
						nodesFromSon1[i]->removeSon(son);
						father->removeSon(nodesFromSon1[i]);
						father->addSon(son);
						son->setDistanceToFather(dist);// BY DEFAULT RIGHT NOW. MAY NEED TO CHANGE IN THE FUTURE
						delete nodesFromSon1[i];
					}
				}
				//	std::cout << "editDuplicationNodesMuffato 8"  <<std::endl;
				
				//Now we need to update the scores of the sons
				if (!node->getSon(0)->isLeaf()) {
					std::vector <Node *> sons = node->getSon(0)->getSons();
					std::vector<int> vec = VectorTools::vectorUnion((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					std::vector<int> vec2 = VectorTools::vectorIntersection((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					node->getSon(0)->setNodeProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() )));
				}
				if (!node->getSon(1)->isLeaf()) {
					std::vector <Node *> sons = node->getSon(1)->getSons();
					std::vector<int> vec = VectorTools::vectorUnion((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					std::vector<int> vec2 = VectorTools::vectorIntersection((dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL(), (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL());
					node->getSon(1)->setNodeProperty("Score", BppString ( TextTools::toString( vec2.size() / vec.size() )));
				}
			}
        }
		
		//continue editing
		editDuplicationNodesMuffato(spTree, 
									geneTree,
									node->getSon(0), 
									editionThreshold);
		editDuplicationNodesMuffato(spTree, 
									geneTree,
									node->getSon(1), 
									editionThreshold);
	}
}
*/




//To sort in descending order
bool cmp( int a, int b ) {
    return a > b;
}  

//To sort in ascending order
bool anticmp( int a, int b ) {
    return a < b;
}  













/*******************************************************************************************************
//////////////////////////////////////////////RUBBISH///////////////////////////////////////////////////
*******************************************************************************************************/









/**************************************************************************
 * This function searches for the best reconciliation while only trying a subset of all roots.
 **************************************************************************/
/*
double findMLReconciliationSafeHeuristics (TreeTemplate<Node> * tree, TreeTemplate<Node> * geneTreeSafe, std::map<std::string, std::string > seqSp, std::vector<int> & lossNumbers, std::vector< double> lossProbabilities, std::vector< int> & duplicationNumbers, std::vector < double> duplicationProbabilities, int & MLindex, std::vector<int> &branchNumbers, int counter, std::vector <int> previousRoots)
{
  TreeTemplate<Node> * geneTree = geneTreeSafe->clone();
  
  std::cout <<"Counter : "<<counter<<std::endl;

  TreeTemplate<Node> * bestGeneTree = NULL; //where we store the most likely gene scenario
  TreeTemplate<Node> * bestSpeciesTree = NULL; //where we store the number of losses and duplications along the species tree

  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  //  geneTree->resetNodesId();
  
 
  std::vector <int> geneNodes;
  double MLRooting = UNLIKELY;
  double currentScore = UNLIKELY;
  // MLindex = 0;
  std::vector <int> bestLossNumbers = lossNumbers;
  std::vector <int> bestDuplicationNumbers = duplicationNumbers;
  std::vector <int> bestBranchNumbers = branchNumbers;
  std::cout << "In findMLReconciliationSafeHeuristics"<<std::endl;
  //We do not want to try all rootings, only rootings close to the current root
  if (MLindex==-1) {
     geneNodes = geneTree->getNodesId();
  }
  else { 
    //    geneNodes = geneTree->getNodesId();

    geneNodes.push_back(MLindex);
    
    for (int i = 0 ; i < geneTree->getNumberOfNodes() ; i++){
      if (!geneTree->getNode(geneTree->getNode(i)->getId())->isLeaf()) {
	geneNodes.push_back(geneTree->getNode(geneTree->getNode(i)->getId())->getId());
      }
    }
  }
 
    
  
  std::cout << "In findMLReconciliationSafeHeuristics : geneNodes.size() : "<< geneNodes.size()<<std::endl;
  //In this loop we try all rootings
  for (int i = 0; i< geneNodes.size(); i++) {
   
    geneTree = geneTreeSafe->clone();

    std::cout <<"Rooting at "<< geneNodes[i]<<std::endl;

    //std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
    // std::cout <<"before geneTree->newOutGroup(i)"<<std::endl;
    if (geneNodes[i] != geneTree->getRootNode()->getId()){
    
      geneTree->newOutGroup(geneNodes[i]);// geneTree->resetNodesId(); 
       std::cout <<"after geneTree->newOutGroup(i)"<<std::endl; std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
     
      resetLossesAndDuplications(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      resetVector(branchNumbers);
      // std::cout << "#######################BEFORE RECONCILE##########################"<<std::endl;
      reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, branchNumbers); 
      // std::cout << "#######################AFTER RECONCILE##########################"<<std::endl;
      currentScore = computeScenarioScore (lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      // std::cout <<"Score : "<< currentScore <<"\n"<<std::endl;
   
      if (currentScore>MLRooting) {
	MLRooting = currentScore;
	MLindex = geneNodes[i];
	bestGeneTree = geneTree->clone();
	bestSpeciesTree = tree->clone();

	bestLossNumbers = lossNumbers;
	bestDuplicationNumbers =duplicationNumbers;
	bestBranchNumbers = branchNumbers;
	int totalLoss = VectorTools::sum(bestLossNumbers);
	int totalDup = VectorTools::sum(bestDuplicationNumbers);
	std::cout <<"New scenario with "<<totalLoss<<" losses and "<<totalDup <<" duplications"<<std::endl;
	
	//  std::cout <<"#######################CHANGE######################"<<std::endl;
      }
    }
  }
  

 
  lossNumbers = bestLossNumbers;
  duplicationNumbers = bestDuplicationNumbers;
  branchNumbers = bestBranchNumbers;
  return MLRooting;
  
  
  
}
*/


/**************************************************************************
 * This function tries all roots on a tree and computes scenarios of duplications and losses and their scores for each root.
 * For this purpose, an exponential law is used to model change probabilities on branches.
 * In the end, the ML value is returned. 
 **************************************************************************/

/*
double findMLReconciliation (TreeTemplate<Node> * tree, TreeTemplate<Node> * geneTree, std::map<std::string, std::string > seqSp, std::vector<int> & lossNumbers, std::vector< double> lossProbabilities, std::vector< int> & duplicationNumbers, std::vector < double> duplicationProbabilities, int & MLindex, std::vector<int> &branchNumbers)
{
  TreeTemplate<Node> * bestGeneTree = NULL; //where we store the most likely gene scenario
  TreeTemplate<Node> * bestSpeciesTree = NULL; //where we store the number of losses and duplications along the species tree
  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  // geneTree->resetNodesId();
  std::vector <int> geneNodes = geneTree->getNodesId();
  double MLRooting = UNLIKELY;
  double currentScore = UNLIKELY;
  // MLindex = 0;
  std::vector <int> bestLossNumbers;
  std::vector <int> bestDuplicationNumbers;
  std::vector <int> bestBranchNumbers;
  //In this loop we try all rootings
  for (int i = 0; i< geneNodes.size(); i++) {
    if (geneTree->isRooted()) {
      geneTree->unroot();
    }
    //    geneTree->resetNodesId();
    //  std::cout <<"Rooting at "<< i<<std::endl;
    if (geneNodes[i] != geneTree->getRootNode()->getId()){
      geneTree->newOutGroup(geneNodes[i]);// geneTree->resetNodesId(); 
      // std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl; 
      resetLossesAndDuplications(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities); 
      resetVector(branchNumbers); 
      reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, branchNumbers); 
      currentScore = computeScenarioScore (lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      // std::cout <<"Score : "<< currentScore <<"\n"<<std::endl;
    
      if (currentScore>MLRooting) {
	MLRooting = currentScore;
	MLindex = geneNodes[i];
	bestGeneTree = geneTree->clone();
	bestSpeciesTree = tree->clone();

	bestLossNumbers = lossNumbers;
	bestDuplicationNumbers =duplicationNumbers;
	bestBranchNumbers = branchNumbers;
	
	//std::cout <<"#######################CHANGE######################"<<std::endl;
      }
    }
  }
  


  *geneTree = *bestGeneTree;
  *tree = *bestSpeciesTree;

  lossNumbers = bestLossNumbers;
  duplicationNumbers = bestDuplicationNumbers;
  branchNumbers = bestBranchNumbers;
  
  return MLRooting;
  
  
  
}
*/

 /**************************************************************************/
 /*
double findMLReconciliation (TreeTemplate<Node> * tree, TreeTemplate<Node> * geneTree, std::map<std::string, std::string > seqSp, std::vector< double> lossProbabilities, std::vector < double> duplicationProbabilities) {
  
  std::vector <int> lossNumbers;
  std::vector <int> duplicationNumbers;
  std::vector <int> branchNumbers;
  
  int numSpNodes = tree->getNumberOfNodes(); 
  for (int i=0; i<numSpNodes; i++) {
    lossNumbers.push_back(0);
    duplicationNumbers.push_back(0);
    branchNumbers.push_back(0);
  }


  int MLindex = 0;
  TreeTemplate<Node> * bestGeneTree = NULL; //where we store the most likely gene scenario
  TreeTemplate<Node> * bestSpeciesTree = NULL; //where we store the number of losses and duplications along the species tree

  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  //  geneTree->resetNodesId();
  

  std::vector <int> geneNodes = geneTree->getNodesId(); 
 
  double MLRooting = UNLIKELY;
  double currentScore = UNLIKELY;
  MLindex = 0;
 
 
  //In this loop we try all rootings
  for (int i = 0; i< geneNodes.size(); i++) {
    if (geneTree->isRooted()) {
      geneTree->unroot();
    }
    //    geneTree->resetNodesId();
    // std::cout <<"Rooting at "<< i<<std::endl; 
    if (geneNodes[i] != geneTree->getRootNode()->getId()){
      geneTree->newOutGroup(geneNodes[i]);// geneTree->resetNodesId(); 
      // std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
      resetLossesAndDuplications(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities); 
      resetVector(branchNumbers);
      reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, branchNumbers);
      currentScore = computeScenarioScore (lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      // std::cout <<"Score : "<< currentScore <<"\n"<<std::endl;
     
      if (currentScore>MLRooting) {
	MLRooting = currentScore;
	MLindex = geneNodes[i];
	bestGeneTree = geneTree->clone();
	bestSpeciesTree = tree->clone();

	
	//std::cout <<"#######################CHANGE######################"<<std::endl;
      }
    }
  }
  

 
  *geneTree = *bestGeneTree;
  *tree = *bestSpeciesTree;

  
  return MLRooting;
  


}
 */

 /**************************************************************************/
 /*
double findMLReconciliationSafe (TreeTemplate<Node> * tree, TreeTemplate<Node> * geneTreeSafe, std::map<std::string, std::string > seqSp, std::vector<int> & lossNumbers, std::vector< double> lossProbabilities, std::vector< int> & duplicationNumbers, std::vector < double> duplicationProbabilities, int & MLindex, std::vector<int> &branchNumbers)
{
  TreeTemplate<Node> * geneTree = geneTreeSafe->clone();
  

  TreeTemplate<Node> * bestGeneTree = NULL; //where we store the most likely gene scenario
  TreeTemplate<Node> * bestSpeciesTree = NULL; //where we store the number of losses and duplications along the species tree

  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  // geneTree->resetNodesId();
  
 
  std::vector <int> geneNodes = geneTree->getNodesId();
  double MLRooting = UNLIKELY;
  double currentScore = UNLIKELY;
  // MLindex = 0;
  std::vector <int> bestLossNumbers = lossNumbers;
  std::vector <int> bestDuplicationNumbers = duplicationNumbers;
  std::vector <int> bestBranchNumbers = branchNumbers;
  //In this loop we try all rootings
  for (int i = 0; i< geneNodes.size(); i++) {
    if (geneTree->isRooted()) {
      geneTree->unroot();
    }
    // geneTree->resetNodesId();
    // std::cout <<"Rooting at "<< i<<std::endl;
    // std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
    // std::cout <<"before geneTree->newOutGroup(i)"<<std::endl;
    if ((geneNodes[i] != geneTree->getRootNode()->getId()) ){
      geneTree->newOutGroup(geneNodes[i]);// geneTree->resetNodesId(); 
      //std::cout <<"after geneTree->newOutGroup(i)"<<std::endl;
     
      resetLossesAndDuplications(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      resetVector(branchNumbers);
      // std::cout << "#######################BEFORE RECONCILE##########################"<<std::endl;
      reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, branchNumbers); 
      // std::cout << "#######################AFTER RECONCILE##########################"<<std::endl;
      currentScore = computeScenarioScore (lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      // std::cout <<"Score : "<< currentScore <<"\n"<<std::endl;
   
      if (currentScore>MLRooting) {
	MLRooting = currentScore;
	MLindex =geneNodes[i];
	bestGeneTree = geneTree->clone();
	bestSpeciesTree = tree->clone();

	bestLossNumbers = lossNumbers;
	bestDuplicationNumbers =duplicationNumbers;
	bestBranchNumbers = branchNumbers;
	int totalLoss = VectorTools::sum(bestLossNumbers);
	int totalDup = VectorTools::sum(bestDuplicationNumbers);
	//	std::cout <<"New scenario with "<<totalLoss<<" losses and "<<totalDup <<" duplications"<<std::endl;
	
	//  std::cout <<"#######################CHANGE######################"<<std::endl;
      }
    }
  }
  


  lossNumbers = bestLossNumbers;
  duplicationNumbers = bestDuplicationNumbers;
  branchNumbers = bestBranchNumbers;
  return MLRooting;
  
  
  
}

 */

 /**************************************************************************/

  /*

double findMLReconciliationSafe (TreeTemplate<Node> * tree, TreeTemplate<Node> * geneTreeSafe, std::map<std::string, std::string > seqSp, std::vector< double> lossProbabilities, std::vector < double> duplicationProbabilities, int & MLindex) {
  
  TreeTemplate<Node> * geneTree = geneTreeSafe->clone();

  std::vector <int> lossNumbers;
  std::vector <int> duplicationNumbers; 
  std::vector <int> branchNumbers; 
  int numSpNodes = tree->getNumberOfNodes(); 
  for (int i=0; i<numSpNodes; i++) {
    lossNumbers.push_back(0);
    duplicationNumbers.push_back(0);
    branchNumbers.push_back(0);
  }

  //int MLindex = 0;
  TreeTemplate<Node> * bestGeneTree = NULL; //where we store the most likely gene scenario
  TreeTemplate<Node> * bestSpeciesTree = NULL; //where we store the number of losses and duplications along the species tree

  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  //  geneTree->resetNodesId();

  
 
  std::vector <int> geneNodes = geneTree->getNodesId();
  double MLRooting = UNLIKELY;
  double currentScore = UNLIKELY;
  // MLindex = 0;
 
  //In this loop we try all rootings
  for (int i = 0; i< geneNodes.size(); i++) {
    if (geneTree->isRooted()) {
      geneTree->unroot();
    }
    //    geneTree->resetNodesId();
    
    // std::cout <<"Rooting at "<< i<<std::endl;
    if (geneNodes[i] != geneTree->getRootNode()->getId()){
      geneTree->newOutGroup(geneNodes[i]);// geneTree->resetNodesId();
      // std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
      resetLossesAndDuplications(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      resetVector(branchNumbers); 
      reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, branchNumbers);
      currentScore = computeScenarioScore (lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities);
      // std::cout <<"Score : "<< currentScore <<"\n"<<std::endl;
    
      if (currentScore>MLRooting) {
	MLRooting = currentScore;
	MLindex = geneNodes[i];
	bestGeneTree = geneTree->clone();
	bestSpeciesTree = tree->clone();

	
	//  std::cout <<"#######################CHANGE######################"<<std::endl;
      }
    }
  }
  
 

 
 
 
  *geneTree = *bestGeneTree;
  *tree = *bestSpeciesTree;

  
  return MLRooting;
  


}


  */







  /*
  if (noeud->isLeaf()) {
    double probaBranch ;
    int num = (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue());
    if (*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D") {
    //do nothing, we already counted the probability of the duplication event on the upper branch !
    std::cout <<"branch "<<noeud->getName()<<" had a duplication before"<<std::endl; 
    probaBranch=0.0;
    }
    else {
       //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
      probaBranch = log(exp(-duplicationExpectations[num])*(1-lossProbabilities[num]));
      branchNumbers[num]++;
    }

    //if not at the root, there may be a long story from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.
    if (noeud->hasFather() && tree.getNode(num)->hasFather()) {
      int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
      int tempNum = num;
      //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
      while ((tree.getNode(tempNum)->hasFather()) && (tree.getNode(tempNum)->getFather()->getId()!=fatherNum)) {
	//There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication
	int tempNumBrother;
	if (tree.getNode(tempNum)->getFather()->getSon(0)->getId()==tempNum) {
	  tempNumBrother = tree.getNode(tempNum)->getFather()->getSon(1)->getId();
	}
	else {
	  tempNumBrother = tree.getNode(tempNum)->getFather()->getSon(0)->getId();
	}
	probaBranch+=log(exp(-duplicationExpectations[tempNumBrother])*(lossProbabilities[tempNumBrother])); 
	//std::cout <<"PROBABRANCH : "<<probaBranch<<" "<<duplicationExpectations[tempNumBrother]<<" "<<lossProbabilities[tempNumBrother]<<"  "<<log(exp(-duplicationExpectations[tempNumBrother])*(lossProbabilities[tempNumBrother]))<<std::endl;
	 branchNumbers[tempNumBrother]++; 
	tempNum = tree.getNode(tempNum)->getFather()->getId();
	//No duplication on branch leading to node now named tempNum, and no loss either
	probaBranch+=log(exp(-duplicationExpectations[tempNum])*(1-lossProbabilities[tempNum]));
	branchNumbers[tempNum]++;
	std::cout << "tempNum "<<tempNum<<std::endl;
      }
    }
    noeud->setBranchProperty(EVENTSPROBA, Number<double>(probaBranch));
    //std::cout << "num : "<<num<<"Probabranch : "<<probaBranch<<std::endl;
    noeud->setNodeProperty(LOWLIK, Number<double>(probaBranch));
  }

  else { //noeud is not a leaf
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);  
    // if ((dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue())==-1)
      // {
	computeScenarioScore (tree, geneTree, son0, branchNumbers, geneNodeIdToSpeciations, duplicationExpectations, lossProbabilities);
	//   }
   // if ((dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue())==-1)
   // {
	computeScenarioScore (tree, geneTree, son1, branchNumbers, geneNodeIdToSpeciations, duplicationExpectations, lossProbabilities);
	// }
    int num = (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue());
    int num0 = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue());
    int num1 = (dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue());
    double probaBranch ;
    if(num==num0||num==num1) {
      int numDup = 0; 
      if (noeud->hasFather()) {
	 Node * temp = noeud;
         while (temp->hasFather() && num==(dynamic_cast<const Number<int> *>(temp->getFather()->getNodeProperty(SPECIESID))->getValue())) {
           numDup++;
           temp = temp->getFather();
         }
      }
      if (numDup==0) {
        //First duplication
        probaBranch = log(exp(-duplicationExpectations[num])*duplicationExpectations[num]*(1-lossProbabilities[num]));
	branchNumbers[num]++;
      }
      else {
	std::cout <<"case where numDup!=0 : "<<numDup<<std::endl;
	//there were duplications before
        probaBranch = log(duplicationExpectations[num]/(numDup+1));	
	//	std::cout <<"before '' " << branchNumbers[num]<<std::endl;
	//  branchNumbers[num]++;
	//  std::cout <<"after '' " << branchNumbers[num]<<std::endl;
      }
    }
   else { //no duplication, no loss on the branch leading to noeud
    
     if (noeud->hasFather()) {
       if (*(dynamic_cast < const std::string* >(noeud->getFather()->getBranchProperty(EVENT)))=="D") {
	 probaBranch = 0.0;
       }
       else {
	 branchNumbers[num]++;
	 probaBranch = log(exp(-duplicationExpectations[num])*(1-lossProbabilities[num]));
       }
     }
     else {
       branchNumbers[num]++;
       probaBranch = log(exp(-duplicationExpectations[num])*(1-lossProbabilities[num]));
     }
   }
    
    //if not at the root, there may be a long story from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.
    if (noeud->hasFather() && tree.getNode(num)->hasFather()) {
      int fatherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getNodeProperty(SPECIESID))->getValue());
      int tempNum = num;
      //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
      std::cout <<"before tempNum "<<tempNum<<std::endl;
      std::cout <<"before fatherNum "<<fatherNum<<std::endl;
      while ((tree.getNode(tempNum)->hasFather()) && (tree.getNode(tempNum)->getFather()->getId()!=fatherNum)) {
	//There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication
	int tempNumBrother;
	if (tree.getNode(tempNum)->getFather()->getSon(0)->getId()==tempNum) {
	  tempNumBrother = tree.getNode(tempNum)->getFather()->getSon(1)->getId();
	}
	else {
	  tempNumBrother = tree.getNode(tempNum)->getFather()->getSon(0)->getId();
	}
	probaBranch+=log(exp(-duplicationExpectations[tempNumBrother])*(lossProbabilities[tempNumBrother])); 
	 branchNumbers[tempNumBrother]++; 
	tempNum = tree.getNode(tempNum)->getFather()->getId();
	//No duplication on branch leading to node now named tempNum, and no loss either
	probaBranch+=log(exp(-duplicationExpectations[tempNum])*(1-lossProbabilities[tempNum]));
	branchNumbers[tempNum]++;
	std::cout << "tempNum "<<tempNum<<std::endl;
      }
    }
    noeud->setBranchProperty(EVENTSPROBA, Number<double>(probaBranch));  
    //std::cout << "num : "<<num<<"Probabranch : "<<probaBranch<<std::endl;
    double lowlik = probaBranch+(dynamic_cast<const Number<double> *>(son0->getNodeProperty(LOWLIK))->getValue())+(dynamic_cast<const Number<double> *>(son1->getNodeProperty(LOWLIK))->getValue());
    noeud->setNodeProperty(LOWLIK, Number<double>(lowlik));
  }
*/
    



















/*
void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::vector<int> &branchNumbers, std::map <int,int> &geneNodeIdToDuplications, std::map <int, std::vector <int> > &geneNodeIdToLosses, std::map <int, std::vector <int> > &geneNodeIdToSpeciations) {
  if (noeud->isLeaf()) {
    //std::cout <<"IS LEAF !!!!!!!!"<<std::endl;
    const int temp = tree.getLeafId(seqSp[noeud->getName()]); 
    branchNumbers[temp]=branchNumbers[temp]+1; 
    //    Number<int> temp2 = new Number(temp);
    noeud->setNodeProperty(SPECIESID, Number<int>(temp));
    noeud->setBranchProperty(EVENT, std::string("S"));
    geneNodeIdToSpeciations[noeud->getId()].push_back(temp);
  }
  else{
    //std::cout << "Root degree : "<<noeud->degree()<<" Number of sons : "<<noeud->getNumberOfSons()<< std::endl;
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);  
    if ((dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue())==-1)
      {
	reconcile (tree, geneTree, son0, seqSp, lossNumbers, duplicationNumbers, branchNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    if ((dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue())==-1)
      {
	reconcile (tree, geneTree, son1, seqSp, lossNumbers, duplicationNumbers, branchNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    int a = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue()); 
    int b = (dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue());
    int a0=a;
    int b0=b;
    int olda=a;
    int oldb=b;
    //Within this loop, it is also possible to get gene losses, and where they occured !
    while (a!=b) { 
      if (a>b) {
	olda=a;
	a = tree.getNode(a)->getFather()->getId();
	//Recording gene losses
	std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, tree.getNode(a)->getSon(0)->getId());
	nodesIds0.push_back(tree.getNode(a)->getSon(0)->getId());
	std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, tree.getNode(a)->getSon(1)->getId());
	nodesIds1.push_back(tree.getNode(a)->getSon(1)->getId());
	if ((tree.getNode(a)->getSon(0)->getId()==olda)&&(!(VectorTools::contains(nodesIds1, b)))&&(b!=a))
	  {
	    int lostNodeId=tree.getNode(a)->getSon(1)->getId();
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
	    lossNumbers[lostNodeId]=lossNumber;
	    branchNumbers[a]=branchNumbers[a]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(a);
	    //std::cout <<"\tA One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
	else if ((tree.getNode(a)->getSon(1)->getId()==olda)&&(!(VectorTools::contains(nodesIds0, b)))&&(b!=a))
	  {
	    int lostNodeId=tree.getNode(a)->getSon(0)->getId();
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
	    lossNumbers[lostNodeId]=lossNumber;
	    branchNumbers[a]=branchNumbers[a]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(a);
	    //std::cout <<"\tB One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
      }
      else { //b>a
	oldb=b;
	b = tree.getNode(b)->getFather()->getId();
	//Recording gene losses
	std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, tree.getNode(b)->getSon(0)->getId());
	nodesIds0.push_back(tree.getNode(b)->getSon(0)->getId());
	std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, tree.getNode(b)->getSon(1)->getId());
	nodesIds1.push_back(tree.getNode(b)->getSon(1)->getId());
	if ((tree.getNode(b)->getSon(0)->getId()==oldb)&&(!(VectorTools::contains(nodesIds1, a)))&&(b!=a))
	  {
	    int lostNodeId=tree.getNode(b)->getSon(1)->getId();
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	    lossNumbers[lostNodeId]=lossNumber;
	    branchNumbers[b]=branchNumbers[b]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(b);
	    //std::cout <<"\tC One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
	else if ((tree.getNode(b)->getSon(1)->getId()==oldb)&&(!(VectorTools::contains(nodesIds0, a)))&&(b!=a))
	  {
	    int lostNodeId=tree.getNode(b)->getSon(0)->getId();
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber)); 
	    lossNumbers[lostNodeId]=lossNumber; 
	    branchNumbers[b]=branchNumbers[b]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(b);
	    //std::cout <<"\tD One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
      }
    }
    //We know a==b
    noeud->setNodeProperty(SPECIESID, Number<int>(a));
    if ((a==a0) || (b==b0)) //There has been a duplication !!
      {
	int dupNumber = (dynamic_cast<const Number<int> *>(tree.getNode(a)->getBranchProperty(DUPLICATIONS))->getValue())+1; 
	tree.getNode(a)->deleteBranchProperty(DUPLICATIONS);
	tree.getNode(a)->setBranchProperty(DUPLICATIONS, Number<int>(dupNumber)); 
	duplicationNumbers[a]=dupNumber;
	//std::cout <<"One duplication on node "<<a<<std::endl;
	geneNodeIdToDuplications[noeud->getId()] = a;
	noeud->setBranchProperty(EVENT, std::string("D"));
	//	geneNodeIdToDuplications[noeud->getId()] = a;
	//std::cout << "\tOne duplication on branch leading to node " << a <<std::endl;
	//The same species node number should be on the two sides of the duplication event, had there been no loss
	branchNumbers[a]=branchNumbers[a]+2;
	if ((a==a0) && (b==b0)) {//there has been no loss, here
	}
	else if (b==b0) { //The loss has occured before a0
	  //We need to place the loss event in the right lineage
	  if (olda==a0) { // only one loss !
	    int lostNodeId;
	    if (tree.getNode(olda)->getFather()->getSon(0)->getId()==olda) {
	      lostNodeId=tree.getNode(olda)->getFather()->getSon(1)->getId();
	    }
	    else {
	      lostNodeId=tree.getNode(olda)->getFather()->getSon(0)->getId();
	    }
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	    lossNumbers[lostNodeId]=lossNumber;
	    branchNumbers[olda]=branchNumbers[olda]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
	    //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
	    //std::cout <<"\t 1 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
	  else {// several losses
	    //get the list of nodes present in the subtree defined by olda
	  
	    if (tree.getNode(a)->getSon(0)->getId()==olda) {
	      int lostNodeId = tree.getNode(a)->getSon(1)->getId();
	      int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	      tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	      tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	      lossNumbers[lostNodeId]=lossNumber;
	      branchNumbers[olda]=branchNumbers[olda]+1;
	      geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	      geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
	      //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
	      //std::cout <<"\t 2 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	    }
	    // else if (VectorTools::contains(nodesIds1, a0)) {
	    else if(tree.getNode(a)->getSon(1)->getId()==olda) {
	      int lostNodeId = tree.getNode(a)->getSon(0)->getId();
	      int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	      tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	      tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	      lossNumbers[lostNodeId]=lossNumber;
	      branchNumbers[olda]=branchNumbers[olda]+1;
	      geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	      geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
	      //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
	      //std::cout <<"\t 3 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	    } 
	    else {std::cout <<"Problem reconcile !!"<<std::endl;
	      exit(-1);
	    }
	  }
	}
	else { 
	  //We need to place the loss event in the right lineage 
	  if (oldb==b0) { // only one loss !
	    int lostNodeId;
	    if (tree.getNode(oldb)->getFather()->getSon(0)->getId()==oldb) {
	      lostNodeId=tree.getNode(oldb)->getFather()->getSon(1)->getId();
	    }
	    else {
	      lostNodeId=tree.getNode(oldb)->getFather()->getSon(0)->getId();
	    }
	    int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	    tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	    tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	    lossNumbers[lostNodeId]=lossNumber;
	    branchNumbers[oldb]=branchNumbers[oldb]+1;
	    geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	    geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
	    //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
	    //std::cout <<"\t 4 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	  }
	  else {// several losses
	    //get the list of nodes present in the subtree defined by oldb
	   
	    if (tree.getNode(a)->getSon(0)->getId()==oldb) {
	      int lostNodeId = tree.getNode(b)->getSon(1)->getId();
	      int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	      tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	      tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	      lossNumbers[lostNodeId]=lossNumber; 
	      branchNumbers[oldb]=branchNumbers[oldb]+1; 
	      geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	      geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
	      //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
// std::cout <<"\t 5 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	    } 
	    // else if (VectorTools::contains(nodesIds1, b0)) {
	    else if(tree.getNode(a)->getSon(1)->getId()==oldb) {
	      int lostNodeId = tree.getNode(b)->getSon(0)->getId();
	      int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
	      tree.getNode(lostNodeId)->deleteBranchProperty(LOSSES);
	      tree.getNode(lostNodeId)->setBranchProperty(LOSSES, Number<int>(lossNumber));
	      lossNumbers[lostNodeId]=lossNumber;
	      branchNumbers[oldb]=branchNumbers[oldb]+1;  
	      geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
	      geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
	      //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
	      //std::cout <<"\t 6 One loss on the branch leading to node " <<lostNodeId<<std::endl;
	    }
	    else {std::cout <<"Problem reconcile !!"<<std::endl;
	      exit(-1);
	    }
	  }
	}

      }
    else 
      {
	branchNumbers[a]+=1;
	noeud->setBranchProperty(EVENT, std::string("S"));
	geneNodeIdToSpeciations[noeud->getId()].push_back(a);
      }
  }
}
*/






	/*
	  for (int i =0; i<branchNumbers.size() ; i++ ) {
	  std::cout <<"branch Number#"<< i<<" du Num: "<< duplicationNumbers[i]<<" du Prob: "<< duplicationProbabilities[i]<<" loss Num: "<< lossNumbers[i]<<" loss Prob: "<< lossProbabilities[i]<<" branch Num: "<< branchNumbers[i]<<std::endl;
	  }
	*/


















/**************************************************************************
 * This function reconciles a gene tree with a species tree with a post-order tree-traversal, using Zmasek and Eddy algorithm (2001). 
 * An added feature is that we also get loss events. 
 * We also fill std::maps containing links between nodes of the gene tree and events on the species tree. 
 * It is important to note that the number of events associated to a given node comes from the values in the "tree" object.
 * Special care must then be used to ascertain these numbers in the "tree" object are correct.
 **************************************************************************/
/*
void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::vector<int> &branchNumbers, std::map <int,int> &geneNodeIdToDuplications, std::map <int, std::vector <int> > &geneNodeIdToLosses, std::map <int, std::vector <int> > &geneNodeIdToSpeciations) {
  if (noeud->isLeaf()) {
    const int temp = tree.getLeafId(seqSp[noeud->getName()]); 
    branchNumbers[temp]=branchNumbers[temp]+1; 
    changeNodeProperty(*noeud, SPECIESID, Number<int>(temp));
    changeBranchProperty(*noeud, EVENT, BppString("S"));
    geneNodeIdToSpeciations[noeud->getId()].push_back(temp);
  }
  else{
    //std::cout << "Root degree : "<<noeud->degree()<<" Number of sons : "<<noeud->getNumberOfSons()<< std::endl;
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);  
    if ((dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue())==-1)
      {
        reconcile (tree, geneTree, son0, seqSp, lossNumbers, duplicationNumbers, branchNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    if ((dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue())==-1)
      {
        reconcile (tree, geneTree, son1, seqSp, lossNumbers, duplicationNumbers, branchNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    int a = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue()); 
    int b = (dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue());
    int a0=a;
    int b0=b;
    int olda=a;
    int oldb=b;
    //Within this loop, it is also possible to get gene losses, and where they occured !
    while (a!=b) { 
      if (a>b) {
        olda=a;
        a = tree.getNode(a)->getFather()->getId();
        //Recording gene losses
        std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, tree.getNode(a)->getSon(0)->getId());
        nodesIds0.push_back(tree.getNode(a)->getSon(0)->getId());
        std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, tree.getNode(a)->getSon(1)->getId());
        nodesIds1.push_back(tree.getNode(a)->getSon(1)->getId());
        if ((tree.getNode(a)->getSon(0)->getId()==olda)&&(!(VectorTools::contains(nodesIds1, b)))&&(b!=a))
          {
            int lostNodeId=tree.getNode(a)->getSon(1)->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(tree.getNode(lostNodeId)), LOSSES,  Number<int>(lossNumber));
            //tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            branchNumbers[a]=branchNumbers[a]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(a);
            //std::cout <<"\tA One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
        else if ((tree.getNode(a)->getSon(1)->getId()==olda)&&(!(VectorTools::contains(nodesIds0, b)))&&(b!=a))
          {
            int lostNodeId=tree.getNode(a)->getSon(0)->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(tree.getNode(lostNodeId)), LOSSES,  Number<int>(lossNumber));
            //tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            branchNumbers[a]=branchNumbers[a]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(a);
            //std::cout <<"\tB One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
      }
      else { //b>a
        oldb=b;
        b = tree.getNode(b)->getFather()->getId();
        //Recording gene losses
        std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, tree.getNode(b)->getSon(0)->getId());
        nodesIds0.push_back(tree.getNode(b)->getSon(0)->getId());
        std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, tree.getNode(b)->getSon(1)->getId());
        nodesIds1.push_back(tree.getNode(b)->getSon(1)->getId());
        if ((tree.getNode(b)->getSon(0)->getId()==oldb)&&(!(VectorTools::contains(nodesIds1, a)))&&(b!=a))
          {
            int lostNodeId=tree.getNode(b)->getSon(1)->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(tree.getNode(lostNodeId)), LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            branchNumbers[b]=branchNumbers[b]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(b);
            //std::cout <<"\tC One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
        else if ((tree.getNode(b)->getSon(1)->getId()==oldb)&&(!(VectorTools::contains(nodesIds0, a)))&&(b!=a))
          {
            int lostNodeId=tree.getNode(b)->getSon(0)->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*( tree.getNode(lostNodeId)), LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber; 
            branchNumbers[b]=branchNumbers[b]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(b);
            //std::cout <<"\tD One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
      }
    }
    
    changeNodeProperty(*(noeud), SPECIESID, Number<int>(a));
    //noeud->setNodeProperty(SPECIESID, Number<int>(a));
    if ((a==a0) || (b==b0)) //There has been a duplication !!
      {
        changeBranchProperty(*(noeud), EVENT, BppString("D"));
        //noeud->setBranchProperty(EVENT, std::string("D"));
        //	geneNodeIdToDuplications[noeud->getId()] = a;
        //std::cout << "\tOne duplication on branch leading to node " << a <<std::endl;
        //The same species node number should be on the two sides of the duplication event, had there been no loss
        if ((a==a0) && (b==b0)) {}//there has been no loss, here
        else if (b==b0) { //The loss has occured before a0
          //We need to place the loss event in the right lineage
          if (olda==a0) { // only one loss !
            int lostNodeId;
            if (tree.getNode(olda)->getFather()->getSon(0)->getId()==olda) {
              lostNodeId=tree.getNode(olda)->getFather()->getSon(1)->getId();
            }
            else {
              lostNodeId=tree.getNode(olda)->getFather()->getSon(0)->getId();
            }
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(tree.getNode(lostNodeId)), LOSSES, Number<int>(lossNumber));
             lossNumbers[lostNodeId]=lossNumber;
            branchNumbers[olda]=branchNumbers[olda]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
            //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
            //std::cout <<"\t 1 One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
          else {// several losses
            //get the list of nodes present in the subtree defined by olda
            if (tree.getNode(a)->getSon(0)->getId()==olda) {
              int lostNodeId = tree.getNode(a)->getSon(1)->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(tree.getNode(lostNodeId)),LOSSES, Number<int>(lossNumber)); 
              lossNumbers[lostNodeId]=lossNumber;
              branchNumbers[olda]=branchNumbers[olda]+1;
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
              //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
              //std::cout <<"\t 2 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            }
            // else if (VectorTools::contains(nodesIds1, a0)) {
            else if(tree.getNode(a)->getSon(1)->getId()==olda) {
              int lostNodeId = tree.getNode(a)->getSon(0)->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(tree.getNode(lostNodeId)),LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber;
              branchNumbers[olda]=branchNumbers[olda]+1;
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
              //  branchNumbers[tree.getNode(olda)->getFather()->getId()]=branchNumbers[tree.getNode(olda)->getFather()->getId()]+1;
              //std::cout <<"\t 3 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            } 
            else {std::cout <<"Problem reconcile !!"<<std::endl;
              exit(-1);
            }
          }
        }
        else { 
          //We need to place the loss event in the right lineage 
          if (oldb==b0) { // only one loss !
            int lostNodeId;
            if (tree.getNode(oldb)->getFather()->getSon(0)->getId()==oldb) {
              lostNodeId=tree.getNode(oldb)->getFather()->getSon(1)->getId();
            }
            else {
              lostNodeId=tree.getNode(oldb)->getFather()->getSon(0)->getId();
            }
            int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(tree.getNode(lostNodeId)),LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            branchNumbers[oldb]=branchNumbers[oldb]+1;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
            //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
            //std::cout <<"\t 4 One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
          else {// several losses
            //get the list of nodes present in the subtree defined by oldb
            if (tree.getNode(a)->getSon(0)->getId()==oldb) {
              int lostNodeId = tree.getNode(b)->getSon(1)->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*( tree.getNode(lostNodeId)), LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber; 
              branchNumbers[oldb]=branchNumbers[oldb]+1; 
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
              //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
              // std::cout <<"\t 5 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            } 
            // else if (VectorTools::contains(nodesIds1, b0)) {
            else if(tree.getNode(a)->getSon(1)->getId()==oldb) {
              int lostNodeId = tree.getNode(b)->getSon(0)->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(tree.getNode(lostNodeId)->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(tree.getNode(lostNodeId)), LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber;
              branchNumbers[oldb]=branchNumbers[oldb]+1;  
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
              //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
              //std::cout <<"\t 6 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            }
            else {std::cout <<"Problem reconcile !!"<<std::endl;
              exit(-1);
            }
          }
        }
        int spId= (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue());
        int dupNumber = (dynamic_cast<const Number<int> *>(tree.getNode(spId)->getBranchProperty(DUPLICATIONS))->getValue())+1; 
        changeBranchProperty(*(tree.getNode(spId)), DUPLICATIONS, Number<int>(dupNumber)); 
        duplicationNumbers[spId]=dupNumber;
        //std::cout <<"One duplication on node "<<spId<<std::endl;
        geneNodeIdToDuplications[noeud->getId()] = spId;
      }
	  else 
      {
        branchNumbers[a]+=1;
        changeBranchProperty(*(noeud), EVENT, BppString("S"));
        //noeud->setBranchProperty(EVENT, std::string("S"));
        geneNodeIdToSpeciations[noeud->getId()].push_back(a);
      }
  }
}
*/


/******************************************************************************************
 * In this version of the function, we do not fill the branch numbers std::vector.
 *****************************************************************************************/
/*
void reconcile (TreeTemplate<Node> & tree, TreeTemplate<Node> & geneTree, Node * noeud, std::map<std::string, std::string > seqSp, std::vector<int >  & lossNumbers, std::vector<int > & duplicationNumbers, std::map <int,int> &geneNodeIdToDuplications, std::map <int, std::vector <int> > &geneNodeIdToLosses, std::map <int, std::vector <int> > &geneNodeIdToSpeciations) {
  if (noeud->isLeaf()) {
    const int temp = tree.getLeafId(seqSp[noeud->getName()]); 
    changeNodeProperty(*(noeud), SPECIESID, Number<int>(temp));
    changeBranchProperty(*(noeud),EVENT, BppString("S")); 
    geneNodeIdToSpeciations[noeud->getId()].push_back(temp);
  }
  else{
    Node * son0=noeud->getSon(0);
    Node * son1=noeud->getSon(1);  
    if ((dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue())==-1)
      {
        reconcile (tree, geneTree, son0, seqSp, lossNumbers, duplicationNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    if ((dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue())==-1)
      {
        reconcile (tree, geneTree, son1, seqSp, lossNumbers, duplicationNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);
      }
    int a = (dynamic_cast<const Number<int> *>(son0->getNodeProperty(SPECIESID))->getValue()); 
    int b = (dynamic_cast<const Number<int> *>(son1->getNodeProperty(SPECIESID))->getValue());
    int a0=a;
    int b0=b;
    int olda=a;
    int oldb=b;
    Node * nodeA = tree.getNode(a);
    Node * oldNodeA = nodeA;
    Node * nodeB;
    if (a!=b) {
      nodeB = tree.getNode(b);
    }
    else {
      nodeB = nodeA;
    }
    Node * oldNodeB = nodeB;
    //Within this loop, it is also possible to get gene losses, and where they occured !
    while (a!=b) { 
      if (a>b) {
        olda=a;
        oldNodeA = nodeA;
        nodeA = nodeA->getFather();
        a = nodeA->getId();
        //Recording gene losses
        std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, nodeA->getSon(0)->getId());
        nodesIds0.push_back(nodeA->getSon(0)->getId());
        std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, nodeA->getSon(1)->getId());
        nodesIds1.push_back(nodeA->getSon(1)->getId());
        if ((nodeA->getSon(0)->getId()==olda)&&(!(VectorTools::contains(nodesIds1, b)))&&(b!=a))
          {
            Node * lostNode = nodeA->getSon(1);
            int lostNodeId=lostNode->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES,  Number<int>(lossNumber));
            //tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(a);
            //std::cout <<"\tA One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
        else if ((nodeA->getSon(1)->getId()==olda)&&(!(VectorTools::contains(nodesIds0, b)))&&(b!=a))
          {
            Node * lostNode = nodeA->getSon(0);
            int lostNodeId=lostNode->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES,  Number<int>(lossNumber));
            // tree.getNode(lostNodeId)->setBranchProperty(LOSSES,  Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(a);
            //std::cout <<"\tB One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
      }
      else { //b>a
        oldb=b;
        oldNodeB = nodeB;
        nodeB = nodeB->getFather();
        b = nodeB->getId();
        //Recording gene losses
        std::vector <int> nodesIds0 = TreeTemplateTools::getNodesId(tree, nodeB->getSon(0)->getId());
        nodesIds0.push_back(nodeB->getSon(0)->getId());
        std::vector <int> nodesIds1 = TreeTemplateTools::getNodesId(tree, nodeB->getSon(1)->getId());
        nodesIds1.push_back(nodeB->getSon(1)->getId());
        if ((nodeB->getSon(0)->getId()==oldb)&&(!(VectorTools::contains(nodesIds1, a)))&&(b!=a))
          {
            Node * lostNode = nodeB->getSon(1);
            int lostNodeId=lostNode->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(b);
            //std::cout <<"\tC One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
        else if ((nodeB->getSon(1)->getId()==oldb)&&(!(VectorTools::contains(nodesIds0, a)))&&(b!=a))
          {
            Node * lostNode = nodeB->getSon(0);
            int lostNodeId = lostNode->getId();
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber)); 
            lossNumbers[lostNodeId]=lossNumber; 
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(b);
            //std::cout <<"\tD One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
      }
    }
    
    changeNodeProperty(*(noeud), SPECIESID, Number<int>(a));
    // noeud->setNodeProperty(SPECIESID, Number<int>(a));
    if ((a==a0) || (b==b0)) //There has been a duplication !!
      {
        //	std::cout << "\n\n\n\t\t\tDUPLICATION on node "<<a<<std::endl;
        //	std::cout << TreeTools::nodeToParenthesis(geneTree, noeud->getId())<<std::endl;
        changeBranchProperty(*(noeud), EVENT, BppString("D"));
        //noeud->setBranchProperty(EVENT, std::string("D"));
        //	geneNodeIdToDuplications[noeud->getId()] = a;
        //std::cout << "\tOne duplication on branch leading to node " << a <<std::endl;
        //The same species node number should be on the two sides of the duplication event, had there been no loss
        if ((a==a0) && (b==b0)) {}//there has been no loss, here
        else if (b==b0) { //The loss has occured before a0
          //We need to place the loss event in the right lineage
          if (olda==a0) { // only one loss !
            int lostNodeId;
            Node * lostNode;
            if (oldNodeA->getFather()->getSon(0)->getId()==olda) {
              lostNode = oldNodeA->getFather()->getSon(1);
              lostNodeId =lostNode->getId();
            }
            else {
              lostNode = oldNodeA->getFather()->getSon(0);
              lostNodeId =lostNode->getId();
            }
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
            //std::cout <<"\t 1 One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
          else {// several losses
            //get the list of nodes present in the subtree defined by olda
            if (nodeA->getSon(0)->getId()==olda) {
              Node * lostNode = nodeA->getSon(1);
              int lostNodeId = lostNode->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber;
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
              
              //std::cout <<"\t 2 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            }
            // else if (VectorTools::contains(nodesIds1, a0)) {
            else if(nodeA->getSon(1)->getId()==olda) {
              Node * lostNode = nodeA->getSon(0);
              int lostNodeId = lostNode->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
               lossNumbers[lostNodeId]=lossNumber;
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(olda);
              //std::cout <<"\t 3 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            } 
            else {std::cout <<"Problem reconcile !!"<<std::endl;
              exit(-1);
            }
          }
        }
        else { 
          //We need to place the loss event in the right lineage 
          if (oldb==b0) { // only one loss !
            Node * lostNode;
            int lostNodeId;
            if (oldNodeB->getFather()->getSon(0)->getId()==oldb) {
              lostNode = oldNodeB->getFather()->getSon(1);
              lostNodeId=lostNode->getId();
            }
            else {
              lostNode = oldNodeB->getFather()->getSon(0);
              lostNodeId=lostNode->getId();
            }
            int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
            changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
            lossNumbers[lostNodeId]=lossNumber;
            geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
            geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
            //std::cout <<"\t 4 One loss on the branch leading to node " <<lostNodeId<<std::endl;
          }
          else {// several losses
            //get the list of nodes present in the subtree defined by oldb
            if (nodeA->getSon(0)->getId()==oldb) {
              Node * lostNode = nodeB->getSon(1);
              int lostNodeId = lostNode->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber; 
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
              //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
              // std::cout <<"\t 5 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            } 
            // else if (VectorTools::contains(nodesIds1, b0)) {
            else if(nodeA->getSon(1)->getId()==oldb) {
              Node * lostNode = nodeB->getSon(0);
              int lostNodeId = lostNode->getId();
              int lossNumber = (dynamic_cast<const Number<int> *>(lostNode->getBranchProperty(LOSSES))->getValue())+1;
              changeBranchProperty(*(lostNode), LOSSES, Number<int>(lossNumber));
              lossNumbers[lostNodeId]=lossNumber;
              geneNodeIdToLosses[noeud->getId()].push_back(lostNodeId);
              geneNodeIdToSpeciations[noeud->getId()].push_back(oldb);
              //  branchNumbers[tree.getNode(oldb)->getFather()->getId()]=branchNumbers[tree.getNode(oldb)->getFather()->getId()]+1;
              //std::cout <<"\t 6 One loss on the branch leading to node " <<lostNodeId<<std::endl;
            }
            else {std::cout <<"Problem reconcile !!"<<std::endl;
              exit(-1);
            }
          }
        }
        int spId= (dynamic_cast<const Number<int> *>(noeud->getNodeProperty(SPECIESID))->getValue());
        Node * nodeSpId = tree.getNode(spId);
        int dupNumber = (dynamic_cast<const Number<int> *>(nodeSpId->getBranchProperty(DUPLICATIONS))->getValue())+1; 
        changeBranchProperty(*(nodeSpId), DUPLICATIONS, Number<int>(dupNumber)); 
        duplicationNumbers[spId]=dupNumber;
        //std::cout <<"One duplication on node "<<spId<<std::endl;
        geneNodeIdToDuplications[noeud->getId()] = spId;
      }
    else 
      {
        changeBranchProperty(*(noeud), EVENT, BppString("S"));
        //noeud->setBranchProperty(EVENT, std::string("S"));
        geneNodeIdToSpeciations[noeud->getId()].push_back(a);
      }
  }
}

*/






/**************************************************************************
 * This function computes a reconciliation for a given root node
 **************************************************************************/
/*
double makeReconciliationAtGivenRoot (TreeTemplate<Node> * tree, 
                                      TreeTemplate<Node> * geneTree, 
                                      std::map<std::string, std::string > seqSp, 
                                      std::vector< double> lossProbabilities, 
                                      std::vector < double> duplicationProbabilities, 
                                      int MLindex)
{
  std::vector <int> lossNumbers;
  std::vector <int> duplicationNumbers; 
  std::vector <int> branchNumbers; 
  std::vector <int> numOlineages; 
  std::vector <int> num1lineages;
  std::vector <int> num2lineages;
  
  
  std::map <int,int> geneNodeIdToDuplications;
  std::map <int, std::vector <int> > geneNodeIdToLosses;
  std::map <int, std::vector <int> > geneNodeIdToSpeciations;
  resetSpeciesIdsAndLiks (*geneTree);
  
  int numSpNodes = tree->getNumberOfNodes(); 
  std::vector <int> vec;
  for (int i=0; i<numSpNodes; i++) {
    lossNumbers.push_back(0);
    duplicationNumbers.push_back(0);
    branchNumbers.push_back(0);
    geneNodeIdToDuplications.insert(std::make_pair(i, -1));
    geneNodeIdToLosses.insert(make_pair(i, vec));
    geneNodeIdToSpeciations.insert(make_pair(i, vec));
  }
  
  double MLRooting = UNLIKELY;
  if (geneTree->isRooted()) {
    geneTree->unroot();
  }
  if (MLindex != geneTree->getRootNode()->getId()){
    geneTree->newOutGroup(MLindex);
    resetLossesAndDuplications(*tree, lossProbabilities, duplicationProbabilities);
    resetVector(branchNumbers); 
    //  std::cout << "reconcile 1!!"<<std::endl;
    reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations); 
    computeScenarioScore (*tree, *geneTree, geneTree->getRootNode(), branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, numOlineages, num1lineages, num2lineages);
    MLRooting = (dynamic_cast<const Number<double> *>(geneTree->getRootNode()->getNodeProperty(LOWLIK))->getValue());
  }
  return MLRooting;
}
*/


/*****************************************************************************
 * Various heuristics are implemented in this function. The choice between the various heuristics depends upon the value of heuristicsLevel.
 * 1 : fastest heuristics : only a few nodes are tried for the roots (the number of the nodes tried depends upon speciesIdLimitForRootPosition), and for each root tried, the events are re-computed only for a subset of the tree.
 * 2 : All roots are tried, and for each root tried, the events are re-computed only for a subset of the tree.
 * 3 : All roots are tried, and for each root tried, the events are re-computed for all nodes of the tree (which should be useless unless there is a bug in the selection of the subset of the nodes.
 ****************************************************************************/
/*
double findMLReconciliation (TreeTemplate<Node> * spTree, 
                             TreeTemplate<Node> * geneTreeSafe, 
                             std::map<std::string, std::string > seqSp, 
                             std::vector<int> & lossNumbers, 
                             std::vector< double> lossProbabilities, 
                             std::vector< int> & duplicationNumbers, 
                             std::vector < double> duplicationProbabilities, 
                             int & MLindex, 
                             std::vector<int> &branchNumbers, 
                             int speciesIdLimitForRootPosition, 
                             int heuristicsLevel, 
                             std::vector <int> &num0lineages, 
                             std::vector <int> &num1lineages, 
                             std::vector <int> &num2lineages, 
                             std::set <int> &nodesToTryInNNISearch)
{
	TreeTemplate<Node> * geneTree = geneTreeSafe->clone();
	TreeTemplate<Node> * tree = spTree->clone();
	if (!geneTree->isRooted()) {
		std::cout << TreeTemplateTools::treeToParenthesis (*geneTree, true)<<std::endl;
		std::cout <<"!!!!!!gene tree is not rooted in findMLReconciliation !!!!!!"<<std::endl;
		exit(-1);
	}
	int oldRoot = geneTree->getRootNode()->getId();
  
	std::vector <int> geneNodes; // contains potential roots to try
	std::vector <Node *> allNodes;
	allNodes = geneTree->getNodes();
	double MLRooting = UNLIKELY;
	double currentScore = UNLIKELY;
	std::vector <int> bestLossNumbers = lossNumbers;
	std::vector <int> bestDuplicationNumbers = duplicationNumbers;
	std::vector <int> bestBranchNumbers = branchNumbers;
	std::vector <int> bestnum0lineages = num0lineages;
	std::vector <int> bestnum1lineages = num1lineages;
	std::vector <int> bestnum2lineages = num2lineages;
  
  std::map <int,int> geneNodeIdToDuplications;
  std::map <int, std::vector <int> > geneNodeIdToLosses;
  std::map <int, std::vector <int> > geneNodeIdToSpeciations; 
	std::vector <int> vec;
	for (int i = 0 ; i< allNodes.size(); i++ ) {
		geneNodeIdToDuplications.insert(std::make_pair(i, -1));
		geneNodeIdToLosses.insert(std::make_pair(i, vec));
		geneNodeIdToSpeciations.insert(std::make_pair(i, vec));
	}
	resetSpeciesIdsAndLiks (*geneTree);
  std::map <int, int > NodeIdToSpId; //correspondence between ids in the gene tree and ids in the species tree
  
	resetLossesAndDuplications(*tree, lossProbabilities, duplicationProbabilities);
	resetVector(branchNumbers);
	resetVector(num0lineages);
	resetVector(num1lineages);
	resetVector(num2lineages);
  reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations); 
	int rootSon1;	
	int rootSon2;
	if (heuristicsLevel > 1) {
		//we try all possible roots
		std::cout <<"All possible roots for gene trees are tried"<<std::endl;
		geneNodes = geneTree->getNodesId();
	}
	else {
		//we only try a smart subset of all nodes
		for (int i = 0; i< allNodes.size() ; i++) {
			int temp = (dynamic_cast<const Number<int> *>(allNodes[i]->getNodeProperty(SPECIESID))->getValue());
      NodeIdToSpId.insert(std::pair <int,int> (allNodes[i]->getId(),temp));
		}
    
		//Now we choose the nodes that we want to try
    rootSon1=geneTree->getRootNode()->getSon(0)->getId();
		rootSon2=geneTree->getRootNode()->getSon(1)->getId();
		int limit = -1;
    std::map<int,int >::iterator iter;
		//This loop aims at finding the smallest species node index in the tree (not all gene trees will have node with species id = 0!)
		for( iter = NodeIdToSpId.begin(); iter != NodeIdToSpId.end(); iter++ ) {
			if ((iter->second<limit)||(limit==-1)) {
				limit = iter->second;
			}
		}
		limit += speciesIdLimitForRootPosition;
    
		for( iter = NodeIdToSpId.begin(); iter != NodeIdToSpId.end(); iter++ ) {
			if ((iter->second <= limit) &&(iter->first!=rootSon1)&&(iter->first!=rootSon2)) {
        geneNodes.push_back(iter->first);
			}
		}
	}
	//We backup this tree and the gene tree, which will constitute our starting point for all roots tried, to avoid computing reconciliation for all nodes. Similarly, we backup the various std::maps.
	TreeTemplate<Node> * backup = geneTree->clone();
	TreeTemplate<Node> * backupTree = tree->clone();
  std::map <int, int> backupGeneNodeIdToDuplications = geneNodeIdToDuplications;
  std::map <int, std::vector <int> > backupGeneNodeIdToLosses = geneNodeIdToLosses;
  std::map <int, std::vector <int> > backupGeneNodeIdToSpeciations = geneNodeIdToSpeciations; 
	std::vector <int> backupLossNumbers = lossNumbers;
	std::vector <int> backupDuplicationNumbers = duplicationNumbers;
	std::vector <int> backupBranchNumbers = branchNumbers;
	std::vector <int> backupnum0lineages = num0lineages;
	std::vector <int> backupnum1lineages = num1lineages;
	std::vector <int> backupnum2lineages = num2lineages;
  
	computeScenarioScore (*tree, *geneTree, geneTree->getRootNode(), branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages);
	currentScore = (dynamic_cast<const Number<double> *>(geneTree->getRootNode()->getNodeProperty(LOWLIK))->getValue());
	//std::cout <<"startingScore : "<<currentScore<<std::endl;	
  
	MLRooting = currentScore;
	MLindex = oldRoot;
	bestLossNumbers = lossNumbers;
	bestDuplicationNumbers =duplicationNumbers;
	bestBranchNumbers = branchNumbers;
	bestnum0lineages = num0lineages;
	bestnum1lineages = num1lineages;
	bestnum2lineages = num2lineages;
	double backupScore=currentScore;
  
  nodesToTryInNNISearch.clear();
  for (int i =0 ; i < allNodes.size(); i++) {
    if (geneNodeIdToDuplications[i]!=-1) {
      nodesToTryInNNISearch.insert(i);
    }
  }
  
	//In this loop we try all rootings, starting from the starting rooting in geneTreeSafe
	for (int i = 0; i< geneNodes.size(); i++) {
		if ((heuristicsLevel != 1) ||(geneNodes[i]!=oldRoot)) {
      delete geneTree;
			delete tree;
			geneTree = backup->clone();
			tree = backupTree->clone();
			//We want to only update nodes whose status has changed between the two roots
			std::vector <int > nodesToUpdate = TreeTools::getPathBetweenAnyTwoNodes(*geneTree, oldRoot, geneNodes[i], true);
			resetSpeciesIdsAndLiksForGivenNodes (*geneTree, nodesToUpdate);
      
			geneNodeIdToDuplications = backupGeneNodeIdToDuplications;
			geneNodeIdToLosses = backupGeneNodeIdToLosses;
			geneNodeIdToSpeciations = backupGeneNodeIdToSpeciations; 
			lossNumbers = backupLossNumbers;
			duplicationNumbers = backupDuplicationNumbers;
			branchNumbers = backupBranchNumbers;
			num0lineages = backupnum0lineages;
			num1lineages = backupnum1lineages;
			num2lineages = backupnum2lineages;
			geneTree->newOutGroup(geneNodes[i]);
			resetLossesDuplicationsSpeciationsForGivenNodes(*tree, lossNumbers, lossProbabilities, duplicationNumbers, duplicationProbabilities, branchNumbers, nodesToUpdate, geneNodeIdToLosses, geneNodeIdToDuplications, geneNodeIdToSpeciations);
      
			if (heuristicsLevel == 3) {
				// Total reset :
				resetSpeciesIdsAndLiks (*geneTree);
				resetLossesAndDuplications(*tree, lossProbabilities, duplicationProbabilities);
				resetVector(branchNumbers);
				resetVector(num0lineages);
				resetVector(num1lineages);
				resetVector(num2lineages);
			}
      
      
			// std::cout << "#######################BEFORE RECONCILE##########################"<<std::endl; 
			reconcile(*tree, *geneTree, geneTree->getRootNode(), seqSp, lossNumbers, duplicationNumbers, geneNodeIdToDuplications, geneNodeIdToLosses, geneNodeIdToSpeciations);  
      
			// std::cout << "#######################AFTER RECONCILE##########################"<<std::endl;
      computeScenarioScore (*tree, *geneTree, geneTree->getRootNode(), branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages);
			currentScore = (dynamic_cast<const Number<double> *>(geneTree->getRootNode()->getNodeProperty(LOWLIK))->getValue());
      
			if (currentScore>MLRooting) {
				MLRooting = currentScore;
				MLindex = geneNodes[i];
				bestLossNumbers = lossNumbers;
				bestDuplicationNumbers =duplicationNumbers;
				bestBranchNumbers = branchNumbers;
				bestnum0lineages = num0lineages;
				bestnum1lineages = num1lineages;
				bestnum2lineages = num2lineages;
        nodesToTryInNNISearch.clear();
        for (int i =0 ; i < allNodes.size(); i++) {
          if (geneNodeIdToDuplications[i]!=-1) {
            nodesToTryInNNISearch.insert(i);
          }
        }
			}
      
		}
  }
  
  // std::cout <<"END OF findMLREconciliation  : MLindex :"<<MLindex<<" LK :"<<MLRooting<<std::endl;
  
  lossNumbers = bestLossNumbers;
  duplicationNumbers = bestDuplicationNumbers;
  branchNumbers = bestBranchNumbers;
  num0lineages = bestnum0lineages;
  num1lineages = bestnum1lineages;
  num2lineages = bestnum2lineages;
  delete tree;
  delete geneTree;
  delete backup;
  delete backupTree;
  
  return MLRooting;
  
  
  
}

*/







