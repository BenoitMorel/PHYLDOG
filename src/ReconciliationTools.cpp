/*
 * Copyright or © or Copr. Centre National de la Recherche Scientifique
 * contributor : Bastien Boussau (2009-2013)
 *
 * bastien.boussau@univ-lyon1.fr
 *
 * This software is a computer program whose purpose is to simultaneously build
 * gene and species trees when gene families have undergone duplications and
 * losses. It can analyze thousands of gene families in dozens of genomes
 * simultaneously, and was presented in an article in Genome Research. Trees and
 * parameters are estimated in the maximum likelihood framework, by maximizing
 * theprobability of alignments given the species tree, the gene trees and the
 * parameters of duplication and loss.
 *
 * This software is governed by the CeCILL license under French law and
 * abiding by the rules of distribution of free software.  You can  use,
 * modify and/ or redistribute the software under the terms of the CeCILL
 * license as circulated by CEA, CNRS and INRIA at the following URL
 * "http://www.cecill.info".
 *
 * As a counterpart to the access to the source code and  rights to copy,
 * modify and redistribute granted by the license, users are provided only
 * with a limited warranty  and the software's author,  the holder of the
 * economic rights,  and the successive licensors  have only  limited
 * liability.
 *
 * In this respect, the user's attention is drawn to the risks associated
 * with loading,  using,  modifying and/or developing or reproducing the
 * software by the user in light of its specific status of free software,
 * that may mean  that it is complicated to manipulate,  and  that  also
 * therefore means  that it is reserved for developers  and  experienced
 * professionals having in-depth computer knowledge. Users are therefore
 * encouraged to load and test the software's suitability as regards their
 * requirements in conditions enabling the security of their systems and/or
 * data to be ensured and,  more generally, to use and operate it in the
 * same conditions as regards security.
 *
 * The fact that you are presently reading this means that you have had
 * knowledge of the CeCILL license and that you accept its terms.
 */
/* This file contains various functions useful for reconciliations, such as reconciliation computation, printing of trees with integer indexes, search of a root with reconciliation...*/

#include "ReconciliationTools.h"

//#include "mpi.h" //SHOULD BE CORRECTED 13062017




void writeReconciledGeneTree ( map<string, string >& params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string>& seqSp, bool temporary ) {
    WHEREAMI( __FILE__ , __LINE__ );
  string suffix = ApplicationTools::getStringParameter ( "output.file.suffix", params, "", "", false, false );
  string reconcTree = ApplicationTools::getStringParameter ( "output.reconciled.tree.file", params, "reconciled.tree", "", false, false );
  reconcTree = reconcTree + suffix;
  if ( temporary ) {
    //   string temp = "temp";
    reconcTree = reconcTree + "_temp";
  }
  writeReconciledGeneTreeToFile (params, geneTree, speciesTree, seqSp, reconcTree);

}



void writeReconciledGeneTreeToFile ( map<string, string >& params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string>& seqSp, string& reconcTree ) {
    WHEREAMI( __FILE__ , __LINE__ );
  Nhx *nhx = new Nhx();
  std::ofstream out;

  breadthFirstreNumber ( *speciesTree );
  std::map <std::string, int> spId = computeSpeciesNamesToIdsMap ( *speciesTree );

  annotateGeneTreeWithDuplicationEvents ( *speciesTree,
                                          *geneTree,
                                          geneTree->getRootNode(),
                                          seqSp,
                                          spId );
  out.open ( reconcTree.c_str(), std::ios::out );
  nhx->write ( *geneTree, out );
  out.close();
  return;
}




/**************************************************************************
 * Assign arbitrary branch lengths (0.1) to branches if they do not have branch lengths
 **************************************************************************/


void assignArbitraryBranchLengths ( TreeTemplate<Node> & tree ) {
  std::vector <Node *> ns = tree.getNodes();
  for ( std::vector <Node *>::iterator it = ns.begin() ; it!= ns.end(); it++ ) {
    if ( ( *it )->hasFather() ) {
      if ( ! ( *it )->hasDistanceToFather() ) {
        ( *it )->setDistanceToFather ( 1 );
      }
    }
  }


}

/**************************************************************************
 * This function re-numbers nodes in a binary tree with a pre-order tree traversal, so that the root has got index 0 and the underlying nodes have higher indexes.
 **************************************************************************/

void reNumber ( TreeTemplate<Node> & tree, Node * noeud, int & index ) {
  noeud->setId ( index );
  index=index+1;
  if ( ! noeud->isLeaf() ) {
    Node * son0=noeud->getSon ( 0 );
    Node * son1=noeud->getSon ( 1 );
    reNumber ( tree, son0, index );
    reNumber ( tree, son1, index );
  }
}



void reNumber ( TreeTemplate<Node> & tree ) {
  Node * root = tree.getRootNode();
  int index=0;
  root->setId ( index );
  index=index+1;
  Node * son0=root->getSon ( 0 );
  Node * son1=root->getSon ( 1 );
  reNumber ( tree, son0, index );
  reNumber ( tree, son1, index );
}

/**************************************************************************
 * This function re-numbers nodes with a breadth-first traversal.
 * 0 =white, 1 = grey, 2 = black
 * It also updates std::vectors of duplication and loss probabilities.
 * It should also return a std::map giving the correspondence between depth levels and node Ids
 **************************************************************************/
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities ) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <double> dupProba=duplicationProbabilities;
  std::vector <double> lossProba=lossProbabilities;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    color.insert ( std::pair <Node *,int> ( nodes[i],0 ) );
  }
  std::queue <Node *> toDo;
  toDo.push ( tree.getRootNode() );
  color[tree.getRootNode()] = 1;
  double dupValue=duplicationProbabilities[tree.getRootNode()->getId()];
  double lossValue=lossProbabilities[tree.getRootNode()->getId()];
  dupProba[index]=dupValue;
  lossProba[index]=lossValue;
  tree.getRootNode()->setId ( index );
  std::vector <int> v;
  DepthToIds.insert ( std::pair <int, std::vector<int> > ( 0,v ) );
  DepthToIds[0].push_back ( index );
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while ( !toDo.empty() ) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for ( unsigned int j = 0 ; j< u->getNumberOfSons() ; j++ ) {
      sons.push_back ( u->getSon ( j ) );
    }
    for ( unsigned int j = 0; j< sons.size() ; j++ ) {
      if ( color[sons[j]]==0 ) {
        color[sons[j]]=1;
        dupValue=duplicationProbabilities[sons[j]->getId()];
        lossValue=lossProbabilities[sons[j]->getId()];
        dupProba[index]=dupValue;
        lossProba[index]=lossValue;
        sons[j]->setId ( index );
        if ( DepthToIds.count ( fatherDepth+1 ) ==0 ) {
          DepthToIds.insert ( std::pair <int, std::vector<int> > ( fatherDepth+1,v ) );
        }
        DepthToIds[fatherDepth+1].push_back ( index );
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push ( sons[j] );
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
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector<double> & duplicationProbabilities, std::vector <double> & lossProbabilities, std::vector <double> & coalBl ) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <double> dupProba=duplicationProbabilities;
  std::vector <double> lossProba=lossProbabilities;
  std::vector <double> cBl=coalBl;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    color.insert ( std::pair <Node *,int> ( nodes[i],0 ) );
  }
  std::queue <Node *> toDo;
  toDo.push ( tree.getRootNode() );
  color[tree.getRootNode()] = 1;
  double dupValue=duplicationProbabilities[tree.getRootNode()->getId()];
  double lossValue=lossProbabilities[tree.getRootNode()->getId()];
  double cValue = coalBl[tree.getRootNode()->getId()];
  dupProba[index]=dupValue;
  lossProba[index]=lossValue;
  cBl[index] = cValue;
  tree.getRootNode()->setId ( index );
  std::vector <int> v;
  DepthToIds.insert ( std::pair <int, std::vector<int> > ( 0,v ) );
  DepthToIds[0].push_back ( index );
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while ( !toDo.empty() ) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for ( unsigned int j = 0 ; j< u->getNumberOfSons() ; j++ ) {
      sons.push_back ( u->getSon ( j ) );
    }
    for ( unsigned int j = 0; j< sons.size() ; j++ ) {
      if ( color[sons[j]]==0 ) {
        color[sons[j]]=1;
        dupValue=duplicationProbabilities[sons[j]->getId()];
        lossValue=lossProbabilities[sons[j]->getId()];
        cValue = coalBl[sons[j]->getId()];
        dupProba[index]=dupValue;
        lossProba[index]=lossValue;
        cBl[index] = cValue;
        sons[j]->setId ( index );
        if ( DepthToIds.count ( fatherDepth+1 ) ==0 ) {
          DepthToIds.insert ( std::pair <int, std::vector<int> > ( fatherDepth+1,v ) );
        }
        DepthToIds[fatherDepth+1].push_back ( index );
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push ( sons[j] );
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
std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree, std::vector <double> & coalBl ) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <double> cBl=coalBl;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    color.insert ( std::pair <Node *,int> ( nodes[i],0 ) );
  }
  std::queue <Node *> toDo;
  toDo.push ( tree.getRootNode() );
  color[tree.getRootNode()] = 1;
  double cValue = coalBl[tree.getRootNode()->getId()];
  cBl[index] = cValue;
  tree.getRootNode()->setId ( index );
  std::vector <int> v;
  DepthToIds.insert ( std::pair <int, std::vector<int> > ( 0,v ) );
  DepthToIds[0].push_back ( index );
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while ( !toDo.empty() ) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for ( unsigned int j = 0 ; j< u->getNumberOfSons() ; j++ ) {
      sons.push_back ( u->getSon ( j ) );
    }
    for ( unsigned int j = 0; j< sons.size() ; j++ ) {
      if ( color[sons[j]]==0 ) {
        color[sons[j]]=1;
        cValue = coalBl[sons[j]->getId()];
        cBl[index] = cValue;
        sons[j]->setId ( index );
        if ( DepthToIds.count ( fatherDepth+1 ) ==0 ) {
          DepthToIds.insert ( std::pair <int, std::vector<int> > ( fatherDepth+1,v ) );
        }
        DepthToIds[fatherDepth+1].push_back ( index );
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push ( sons[j] );
      }
    }
    color[u]=2;
  }
  coalBl = cBl;
  return DepthToIds;
}


//There we do not update probabilities

std::map <int, std::vector <int> > breadthFirstreNumber ( TreeTemplate<Node> & tree ) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    color.insert ( std::pair <Node *,int> ( nodes[i],0 ) );
  }
  std::queue <Node *> toDo;
  toDo.push ( tree.getRootNode() );
  color[tree.getRootNode()] = 1;

  tree.getRootNode()->setId ( index );

  std::vector <int> v;
  DepthToIds.insert ( std::pair <int, std::vector<int> > ( 0,v ) );
  DepthToIds[0].push_back ( index );
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while ( !toDo.empty() ) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for ( unsigned int j = 0 ; j< u->getNumberOfSons() ; j++ ) {
      sons.push_back ( u->getSon ( j ) );
    }
    for ( unsigned int j = 0; j< sons.size() ; j++ ) {
      if ( color[sons[j]]==0 ) {
        color[sons[j]]=1;
        sons[j]->setId ( index );
        if ( DepthToIds.count ( fatherDepth+1 ) ==0 ) {
          DepthToIds.insert ( std::pair <int, std::vector<int> > ( fatherDepth+1,v ) );
        }
        DepthToIds[fatherDepth+1].push_back ( index );
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push ( sons[j] );
      }
    }
    color[u]=2;
  }
  return DepthToIds;
}

//There we do not update probabilities, and we reset "toComp" node properties
std::map <int, std::vector <int> > breadthFirstreNumberAndResetProperties ( TreeTemplate<Node> & tree ) {
  int index = 0;
  std::map<Node *, int> color ;
  std::map <int, std::vector <int> > DepthToIds; //A std::map where we store the correspondence between the depth of a node (number of branches between the root and the node) and the node id.
  std::map <int, int > IdsToDepths;
  std::vector <Node * > nodes = tree.getNodes();
  //All nodes white
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    color.insert ( std::pair <Node *,int> ( nodes[i],0 ) );
    if ( nodes[i]->hasNodeProperty ( "toComp" ) ) {
      nodes[i]->setNodeProperty ( "toComp", BppString ( "y" ) );
    }
  }
  std::queue <Node *> toDo;
  toDo.push ( tree.getRootNode() );
  color[tree.getRootNode()] = 1;
  tree.getRootNode()->setId ( index );
  std::vector <int> v;
  DepthToIds.insert ( std::pair <int, std::vector<int> > ( 0,v ) );
  DepthToIds[0].push_back ( index );
  IdsToDepths[index] = 0;
  index++;
  Node * u;
  while ( !toDo.empty() ) {
    u = toDo.front();
    toDo.pop();
    int fatherDepth = IdsToDepths[u->getId()];
    std::vector <Node *> sons;
    for ( unsigned int j = 0 ; j< u->getNumberOfSons() ; j++ ) {
      sons.push_back ( u->getSon ( j ) );
    }
    for ( unsigned int j = 0; j< sons.size() ; j++ ) {
      if ( color[sons[j]]==0 ) {
        color[sons[j]]=1;
        sons[j]->setId ( index );
        if ( DepthToIds.count ( fatherDepth+1 ) ==0 ) {
          DepthToIds.insert ( std::pair <int, std::vector<int> > ( fatherDepth+1,v ) );
        }
        DepthToIds[fatherDepth+1].push_back ( index );
        IdsToDepths[index] = fatherDepth+1;
        index++;
        toDo.push ( sons[j] );
      }
    }
    color[u]=2;
  }
  return DepthToIds;
}




/**************************************************************************
 * These functions reset the LOSSES and DUPLICATIONS values on the species tree and in the associated std::vectors.
 **************************************************************************/

void changeNodeProperty ( Node & noeud, const std::string & name, const Clonable & property ) {
  if ( noeud.hasNodeProperty ( name ) ) {
    noeud.deleteNodeProperty ( name );
    noeud.setNodeProperty ( name, property );
  }
  else {
    noeud.setNodeProperty ( name, property );
  }
}

void changeBranchProperty ( Node & noeud, const std::string & name, const Clonable & property ) {
  if ( noeud.hasBranchProperty ( name ) ) {
    noeud.deleteBranchProperty ( name );
    noeud.setBranchProperty ( name, property );
  }
  else {
    noeud.setBranchProperty ( name, property );
  }
}



void resetLossesAndDuplications ( TreeTemplate<Node> & tree, /*std::vector <int> &lossNumbers, */std::vector <double> &lossProbabilities, /*std::vector <int> &duplicationNumbers, */std::vector <double> &duplicationProbabilities ) {
  std::vector< int > nodesIds = tree.getNodesId ();
  //  Number<int> zero = Number<int>(0);
  BppString zero = BppString ( TextTools::toString ( 0 ) );
  for ( unsigned int i=0; i<nodesIds.size(); i++ ) {
    changeBranchProperty ( * ( tree.getNode ( nodesIds[i] ) ),LOSSES, zero );
    changeBranchProperty ( * ( tree.getNode ( nodesIds[i] ) ),DUPLICATIONS, zero );
    /*lossNumbers[nodesIds[i]] = 0;
     *        duplicationNumbers[nodesIds[i]] = 0;*/
  }
}

void resetLossesDuplicationsSpeciations ( TreeTemplate<Node> & tree, std::vector <int> &lossNumbers, std::vector <double> &lossProbabilities, std::vector <int> &duplicationNumbers, std::vector <double> &duplicationProbabilities, std::vector <int> &branchNumbers ) {
  std::vector< int > nodesIds = tree.getNodesId ();
  // Number<int> zero = Number<int>(0);
  BppString zero = BppString ( TextTools::toString ( 0 ) );

  for ( unsigned int i=0; i<nodesIds.size(); i++ ) {
    changeBranchProperty ( * ( tree.getNode ( nodesIds[i] ) ),LOSSES, zero );
    changeBranchProperty ( * ( tree.getNode ( nodesIds[i] ) ),DUPLICATIONS, zero );
    lossNumbers[nodesIds[i]] = 0;
    duplicationNumbers[nodesIds[i]] = 0;
    branchNumbers[nodesIds[i]] = 0;
  }
}



void resetLossesDuplicationsSpeciationsForGivenNodes ( TreeTemplate<Node> & tree, std::vector <int> & lossNumbers, std::vector <double> & lossProbabilities, std::vector <int> & duplicationNumbers, std::vector <double> & duplicationProbabilities, std::vector <int> & branchNumbers, std::vector <int> nodesToUpdate, std::map <int, std::vector<int> > & geneNodeIdToLosses, std::map <int, int > & geneNodeIdToDuplications, std::map <int, std::vector<int> > & geneNodeIdToSpeciations ) {
  for ( unsigned int i=0; i<nodesToUpdate.size(); i++ ) {
    unsigned int size = geneNodeIdToLosses[nodesToUpdate[i]].size();
    for ( unsigned int j = 0 ; j< size ; j++ ) {
      lossNumbers[geneNodeIdToLosses[nodesToUpdate[i]][j]]--;
      changeBranchProperty ( * ( tree.getNode ( geneNodeIdToLosses[nodesToUpdate[i]][j] ) ),LOSSES, BppString ( TextTools::toString ( ( lossNumbers[geneNodeIdToLosses[nodesToUpdate[i]][j]] ) ) ) );
    }
    size = geneNodeIdToSpeciations[nodesToUpdate[i]].size();
    for ( unsigned int j = 0 ; j< size ; j++ ) {
      branchNumbers[geneNodeIdToSpeciations[nodesToUpdate[i]][j]]--;
    }
    if ( geneNodeIdToDuplications[nodesToUpdate[i]] != -1 ) {
      duplicationNumbers[geneNodeIdToDuplications[nodesToUpdate[i]]]--;
      changeBranchProperty ( * ( tree.getNode ( geneNodeIdToDuplications[nodesToUpdate[i]] ) ), DUPLICATIONS, BppString ( TextTools::toString ( ( duplicationNumbers[geneNodeIdToDuplications[nodesToUpdate[i]]] ) ) ) );
    }
  }


  //We do not reset the 3 std::maps, because it is currently not useful. It might be useful provided the program is changed in some way.
}

/**************************************************************************
 * This function sets all values in a std::vector<int> to 0.
 **************************************************************************/

void printVector ( std::vector<int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    std::cout <<  "i : "<<i<< " vi : "<< v[i]<<std::endl;
  }
}

void printVectorLine ( std::vector<int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    std::cout <<  "i : "<<i<< " vi : "<< v[i]<<" ";
  }
  std::cout <<std::endl;
}

void resetVector ( std::vector<unsigned int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0;
  }
}

void resetVector ( std::vector<int> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0;
  }
}

void resetVector ( std::vector<double> & v ) {
  unsigned int temp=v.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[i]=0.0;
  }
}


void resetVectorForGivenNodes ( std::vector<int> & v, std::vector <int> nodesToUpdate ) {
  unsigned int temp=nodesToUpdate.size();
  for ( unsigned int i=0; i<temp ; i++ ) {
    v[nodesToUpdate[i]]=0;
  }
}

/**************************************************************************
 * This function resets species Ids in a tree.
 **************************************************************************/
void resetSpeciesIds ( TreeTemplate<Node> & tree ) {
  std::vector <Node *> nodes = tree.getNodes();
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    changeNodeProperty ( * ( nodes[i] ), SPECIESID, Number<int> ( -1 ) );
  }
}

void resetSpeciesIdsForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds ) {
  for ( unsigned int i = 0; i< nodesToUpdate.size() ; i++ ) {
    removedNodeIds.push_back ( ( dynamic_cast<const Number<int> *> ( tree.getNode ( nodesToUpdate[i] )->getNodeProperty ( SPECIESID ) )->getValue() ) );
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), SPECIESID, Number<int> ( -1 ) );
  }
}

void resetSpeciesIdsForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate ) {
  for ( unsigned int i = 0; i< nodesToUpdate.size() ; i++ ) {
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), SPECIESID, Number<int> ( -1 ) );
  }
}

void resetSpeciesIdsAndLiks ( TreeTemplate<Node> & tree ) {
  std::vector <Node *> nodes = tree.getNodes();
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    changeNodeProperty ( * ( nodes[i] ), SPECIESID, Number<int> ( -1 ) );
    changeNodeProperty ( * ( nodes[i] ), LOWLIK, Number<double> ( 1.0 ) );
    changeBranchProperty ( * ( nodes[i] ),EVENTSPROBA, Number<double> ( 1.0 ) );
  }
}

void resetSpeciesIdsAndLiksForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate, std::vector <int> & removedNodeIds ) {
  for ( unsigned int i = 0; i< nodesToUpdate.size() ; i++ ) {
    removedNodeIds.push_back ( ( dynamic_cast<const Number<int> *> ( tree.getNode ( nodesToUpdate[i] )->getNodeProperty ( SPECIESID ) )->getValue() ) );
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), SPECIESID, Number<int> ( -1 ) );
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), LOWLIK, Number<double> ( 1.0 ) );
    changeBranchProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), EVENTSPROBA, Number<double> ( 1.0 ) );
  }
}

void resetSpeciesIdsAndLiksForGivenNodes ( TreeTemplate<Node> & tree, std::vector<int > nodesToUpdate ) {
  for ( unsigned int i = 0; i< nodesToUpdate.size() ; i++ ) {
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), SPECIESID, Number<int> ( -1 ) );
    changeNodeProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), LOWLIK, Number<double> ( 1.0 ) );
    changeBranchProperty ( * ( tree.getNode ( nodesToUpdate[i] ) ), EVENTSPROBA, Number<double> ( 1.0 ) );
  }
}


/**************************************************************************
 * This function uses a std::map giving the link between levels and nodes ids and sends the maximum id corresponding to a limit level.
 **************************************************************************/

int getLimitIdFromNumberOfLevels ( std::map<int, std::vector <int> > DepthToIds, int levelsToRootAt ) {
  return ( VectorTools::max ( DepthToIds[levelsToRootAt] ) );
}



/**************************************************************************
 * Computes the rates of duplication and loss of the birth-death process along a
 * particular branch given counts of times where, at the end of the branch,
 * there were 0 gene (i), 1 gene (j) or 2 genes (k).
 * This is therefore an approximate formula (normally,
 * one should consider all possible counts, from 0 to +infinity).
 **************************************************************************/
void computeDuplicationAndLossProbabilities (const int i, const int j, const int k,
                                             double & lossProbability,
                                             double & duplicationProbability )
{
  double id=double ( i );
  double jd=double ( j );
  double kd=double ( k );

  //Linear Birth-Death process, from Bartholomay 1958
  double sum = id +jd +kd;
  double X = id/sum;
  double Y = jd/sum;
  double Z = kd/sum;
  //num2Lineages (Z) corresponds to the sum of cases where there were 2 or
  //more than 2 lineages at the end of a branch.
  double ZoverYPlusZ = Z/ ( Y+Z );
  double ZoverXTimesYPlusZ = ZoverYPlusZ * 1/X;

  //expected number of losses over a branch of length 1.
  double mu = log ( ( X-1 ) / ( ZoverYPlusZ-1 ) ) / ( ZoverXTimesYPlusZ -1 );

  //expected number of duplications over a branch of lenth 1.
  double lambda = ZoverXTimesYPlusZ * mu;

  /* If we truncate
   *      double Y2Z = pow(Y, 2)/Z;
   *      double XYZ = X*Y/Z;
   *      double oneXYZ = 1/XYZ;
   *
   *      //expected number of losses over a branch of length 1.
   *      double mu = (log(XYZ - X + Y2Z) - log(1 - X + Y2Z) ) / (1- oneXYZ);
   *
   *      //expected number of duplications over a branch of lenth 1.
   *      double lambda = oneXYZ * mu;
   */
  if ( ! ( std::isnan ( lambda ) ||std::isinf ( lambda ) ) ) {
    if ( lambda<=0.000001 ) {
      duplicationProbability = 0.000001;//SMALLPROBA;
    }
    else {
      duplicationProbability = lambda;
    }
  }
  else if ( std::isnan ( lambda ) ) {
    duplicationProbability = 0.000001;
  }



  if ( ! ( std::isnan ( mu ) ||std::isinf ( mu ) ) ) {
    if ( mu<=0.0000011 ) {
      lossProbability = 0.0000011;//SMALLPROBA;
    }
    else {
      lossProbability = mu;
    }
  }
  else if ( std::isnan ( mu ) ) {
    std::cerr<<"A branch loss probability could not be properly estimated. Set to 0.000001 by default. X:"<<X<<" Y:"<<Y<<" Z:"<<Z<<std::endl;
    lossProbability = 0.0000011;
  }

  return;
}

/*************************************************
 * For the root branch, a special trick is used. At the root, we never count families that have losses, because these families do not have any gene from the species involved in the study. So we set the number of losses at the average number observed in other branches.
 *************************************************/



void computeDuplicationAndLossProbabilitiesForAllBranches ( const std::vector <int> &numOGenes, const std::vector <int> &num1Genes, const std::vector <int> &num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities ) {
  //The trick for the root, already done in alterLineages...:
  /*
   *     int totNum0=0, totNum12=0;
   *     for (int i =1 ; i< lossProbabilities.size() ; i++) {
   *       totNum0+=numOGenes[i];
   *       totNum12+=num1Genes[i]+num2Genes[i];
}
//At branch 0, by definition, we never count cases where there has been a loss, so we set it to an average value.
numOGenes[0] = (int)floor( ((double)totNum0/(double)totNum12)*((double)num1Genes[0]+(double)num2Genes[0]));
*/
  for ( unsigned int i =0 ; i< lossProbabilities.size() ; i++ ) {
    computeDuplicationAndLossProbabilities ( numOGenes[i], num1Genes[i], num2Genes[i], lossProbabilities[i], duplicationProbabilities[i] );
  }
  return;
}


void computeAverageDuplicationAndLossProbabilitiesForAllBranches ( const std::vector <int> &numOGenes, const std::vector <int> &num1Genes, const std::vector <int> &num2Genes, std::vector <double> & lossProbabilities, std::vector<double> & duplicationProbabilities ) {
  //The trick for the root:
  /* int totNum0=0, totNum12=0;
   *     for (int i =1 ; i< lossProbabilities.size() ; i++) {
   *       totNum0+=numOGenes[i];
   *       totNum12+=num1Genes[i]+num2Genes[i];
}
numOGenes[0] = (int)floor( ((double)totNum0/(double)totNum12)*((double)num1Genes[0]+(double)num2Genes[0]));
*/
  int sumOGene = VectorTools::sum ( numOGenes );
  int sum1Gene = VectorTools::sum ( num1Genes );
  int sum2Gene = VectorTools::sum ( num2Genes );
  double dupProba, lossProba;
  computeDuplicationAndLossProbabilities ( sumOGene, sum1Gene, sum2Gene, lossProba, dupProba );
  int size = lossProbabilities.size();
  if ( ! ( std::isnan ( lossProba ) ||std::isinf ( lossProba ) || ( lossProba<0.0 ) ) ) {
    lossProbabilities.assign ( size, lossProba );
  }
  if ( ! ( std::isnan ( dupProba ) ||std::isinf ( dupProba ) || ( dupProba<0.0 ) ) ) {
    duplicationProbabilities.assign ( size, dupProba );
  }
  return;
}



/**************************************************************************
 * Computes the probability of a given number of lineages numberOfLineages at the end of a branch given that there was only one at the beginning of the branch.
 **************************************************************************/
double computeBranchProbability ( const double & duplicationProbability, const double & lossProbability, const int numberOfLineages ) {
  //Using a linear birth-death process (Bartholomay, 1958)
  /* int min ;
   *     if (numberOfLineages == 0) {
   *       min = 0;
}
else {
  min = 1;
}
*/
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
double computeLogBranchProbability ( const double & duplicationProbability, const double & lossProbability, const int numberOfLineages ) {
  //UNDONE TEST 1409
  /*  if (numberOfLineages<=1) {
   *        return 0;
}
else {
  return(-numberOfLineages+1);*/
  return ( log ( computeBranchProbability ( duplicationProbability, lossProbability, numberOfLineages ) ) );
  // }
}


/**************************************************************************
 * Computes the probability of a given number of lineages numberOfLineages at the end of a branch given that there was 0 at the beginning of the branch.
 **************************************************************************/
double computeBranchProbabilityAtRoot ( const double &duplicationProbability, double &lossProbability, const int numberOfLineages ) {
  if ( duplicationProbability==lossProbability ) {
    lossProbability+=0.0000001;
  }
  double beta = ( 1-exp ( duplicationProbability-lossProbability ) ) / ( lossProbability-duplicationProbability*exp ( duplicationProbability-lossProbability ) );
  if ( std::isnan ( beta ) ) {
    std::cout <<"ISNAN beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
  }
  if ( numberOfLineages==0 ) {
    std::cerr <<"Error in computeBranchProbabilityAtRoot: cannot compute P(0 lineage)!"<<std::endl;
    exit ( -1 );
    //    return (lossProbability*beta);
  }
  else {
    double temp = ( 1-duplicationProbability*beta ) *pow ( ( duplicationProbability*beta ), ( numberOfLineages-1 ) );
    if ( std::isnan ( temp ) ) {
      std::cout <<"ISNAN2 beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
    }
    else if ( temp == 0 ) {
      // std::cerr <<"TEMP equals 0; beta :"<<beta<<" duplicationProbability :"<<duplicationProbability<<" lossProbability:"<<lossProbability<<" numberOfLineages :"<<numberOfLineages<<std::endl;
      temp = NumConstants::VERY_TINY() ;
    }
    return ( temp );
  }
}


/**************************************************************************
 * Computes the log of the above probability (if it had not been changed not to!):
 * Update as of January 28th, 2010.
 * Now we do not consider gene originations anymore.
 * Therefore events at the root are considered as on any other branch.
 **************************************************************************/
double computeLogBranchProbabilityAtRoot ( const double &duplicationProbability, const double &lossProbability, const int numberOfLineages ) {
  //UNDONE TEST 1409
  /* if (numberOfLineages<=1) {
   *       return 0;
}
else {
  return(-numberOfLineages+1);
}*/
  //Update as of January 28th, 2010.
  //Now we do not consider gene originations anymore.
  //Therefore, we do not make a special case when we are at the root.
  // return(log(computeBranchProbabilityAtRoot(duplicationProbability, lossProbability, numberOfLineages)));
  //instead, we use:
  return ( log ( computeBranchProbability ( duplicationProbability, lossProbability, numberOfLineages ) ) );
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

void computeScenarioScore ( TreeTemplate<Node> & tree,
                            TreeTemplate<Node> & geneTree,
                            Node * noeud,
                            std::vector<int> &branchNumbers,
                            std::map <int, std::vector <int> > &geneNodeIdToSpeciations,
                            const std::vector<double> &duplicationProbabilities,
                            const std::vector<double> &lossProbabilities,
                            std::vector <int> &num0lineages,
                            std::vector <int> &num1lineages,
                            std::vector <int> &num2lineages )
{
  int num = ( dynamic_cast<const Number<int> *> ( noeud->getNodeProperty ( SPECIESID ) )->getValue() );
  Node * spNode = tree.getNode ( num );
  if ( noeud->isLeaf() ) {
    int fatherNum = ( dynamic_cast<const Number<int> *> ( noeud->getFather()->getNodeProperty ( SPECIESID ) )->getValue() );
    double probaBranch ;
    if ( * ( dynamic_cast < const std::string* > ( noeud->getFather()->getBranchProperty ( EVENT ) ) ) =="D" ) {
      //do nothing, we already counted the probability of the duplication event on one of the upper branches !
      probaBranch=0.0;
      changeNodeProperty ( * ( noeud ), NUMGENES, Number<int> ( 1 ) );
      if ( num!=fatherNum ) {
        //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
        probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], 1 );
        num1lineages[num]++;
        branchNumbers[num]++;
      }
    }
    else {
      //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
      probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], 1 );
      changeNodeProperty ( * ( noeud ), NUMGENES, Number<int> ( 1 ) );
      num1lineages[num]++;
      branchNumbers[num]++;
    }
    //if not at the root, there may be a long story of losses from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.

    if ( spNode->hasFather() ) {
      int tempNum = num;
      Node * spTempNode = spNode;
      //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
      int tempNumBrother;
      while ( ( tempNum!=fatherNum ) && ( spTempNode->hasFather() ) && ( spTempNode->getFather()->getId() !=fatherNum ) ) {
        //There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication
        if ( spTempNode->getFather()->getSon ( 0 )->getId() ==tempNum ) {
          tempNumBrother = spTempNode->getFather()->getSon ( 1 )->getId();
        }
        else {
          tempNumBrother = spTempNode->getFather()->getSon ( 0 )->getId();
        }

        probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0 );
        num0lineages[tempNumBrother]++;
        branchNumbers[tempNumBrother]++;
        spTempNode = spTempNode->getFather();
        tempNum = spTempNode->getId();
        //No duplication on branch leading to node now named tempNum, and no loss either
        probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNum], lossProbabilities[tempNum], 1 );
        num1lineages[tempNum]++;
        branchNumbers[tempNum]++;
      }
      if ( ( * ( dynamic_cast < const std::string* > ( noeud->getFather()->getBranchProperty ( EVENT ) ) ) =="D" ) && ( num!=fatherNum ) ) {
        if ( spTempNode->getFather()->getSon ( 0 )->getId() ==tempNum ) {
          tempNumBrother = spTempNode->getFather()->getSon ( 1 )->getId();
        }
        else {
          tempNumBrother = spTempNode->getFather()->getSon ( 0 )->getId();
        }
        probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0 );
        num0lineages[tempNumBrother]++;
        branchNumbers[tempNumBrother]++;
      }
    }

    changeBranchProperty ( * ( noeud ), EVENTSPROBA, Number<double> ( probaBranch ) );
    changeNodeProperty ( * ( noeud ), LOWLIK, Number<double> ( probaBranch ) );
  }

  else { //noeud is not a leaf
    Node * son0=noeud->getSon ( 0 );
    Node * son1=noeud->getSon ( 1 );
    if ( ( dynamic_cast<const Number<double> *> ( son0->getNodeProperty ( LOWLIK ) )->getValue() ) ==1.0 ) {
      computeScenarioScore ( tree, geneTree, son0, branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages );
    }
    if ( ( dynamic_cast<const Number<double> *> ( son1->getNodeProperty ( LOWLIK ) )->getValue() ) ==1.0 ) {
      computeScenarioScore ( tree, geneTree, son1, branchNumbers, geneNodeIdToSpeciations, duplicationProbabilities, lossProbabilities, num0lineages, num1lineages, num2lineages );
    }
    int num0 = ( dynamic_cast<const Number<int> *> ( son0->getNodeProperty ( SPECIESID ) )->getValue() );
    int num1 = ( dynamic_cast<const Number<int> *> ( son1->getNodeProperty ( SPECIESID ) )->getValue() );
    double probaBranch ;
    if ( num==num0||num==num1 ) { //there is a duplication at this branch
      if ( noeud->hasFather() ) {
        int fatherNum = ( dynamic_cast<const Number<int> *> ( noeud->getFather()->getNodeProperty ( SPECIESID ) )->getValue() );
        /* int brotherNum;
         *                 if (noeud->getFather()->getSon(0)->getId()==noeud->getId()) {
         *                   brotherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getSon(1)->getNodeProperty(SPECIESID))->getValue());
      }
      else {
        brotherNum = (dynamic_cast<const Number<int> *>(noeud->getFather()->getSon(0)->getNodeProperty(SPECIESID))->getValue());
      }*/

        if ( fatherNum == num ) {
          //this is not the first duplication: we cannot compute the score yet
          int numGenes = ( dynamic_cast<const Number<int> *> ( son0->getNodeProperty ( NUMGENES ) )->getValue() ) + ( dynamic_cast<const Number<int> *> ( son1->getNodeProperty ( NUMGENES ) )->getValue() );
          changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( numGenes ) );
          probaBranch = 0.0;
        }
        else {
          //this is the first duplication, we compute the probability
          branchNumbers[num]++;
          int numGenes = ( dynamic_cast<const Number<int> *> ( son0->getNodeProperty ( NUMGENES ) )->getValue() ) + ( dynamic_cast<const Number<int> *> ( son1->getNodeProperty ( NUMGENES ) )->getValue() );
          probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], numGenes );
          if ( numGenes==2 ) {
            num2lineages[num]++;
          }
          changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( 1 ) );
        }
      }
      else {//A duplication at the root
        //this is the first duplication, we compute the probability
        branchNumbers[num]++;
        int numGenes = ( dynamic_cast<const Number<int> *> ( son0->getNodeProperty ( NUMGENES ) )->getValue() ) + ( dynamic_cast<const Number<int> *> ( son1->getNodeProperty ( NUMGENES ) )->getValue() );
        probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], numGenes );
        if ( numGenes==2 ) {
          num2lineages[num]++;
        }
        changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( 1 ) );
      }
    }
    else { //no duplication, no loss on the branch leading to noeud

      if ( noeud->hasFather() ) {
        int fatherNum = ( dynamic_cast<const Number<int> *> ( noeud->getFather()->getNodeProperty ( SPECIESID ) )->getValue() );
        if ( * ( dynamic_cast < const std::string* > ( noeud->getFather()->getBranchProperty ( EVENT ) ) ) =="D" ) {
          probaBranch = 0.0;
          changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( 1 ) );
          if ( num!=fatherNum ) {
            //Probability of events along that branch : this was a simple speciation, so no duplication, no loss
            probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], 1 );
            num1lineages[num]++;
            branchNumbers[num]++;
          }
        }
        else {
          branchNumbers[num]++;
          probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], 1 );
          num1lineages[num]++;
          changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( 1 ) );
        }
      }
      else {
        branchNumbers[num]++;
        probaBranch = computeLogBranchProbability ( duplicationProbabilities[num], lossProbabilities[num], 1 );
        num1lineages[num]++;
        changeNodeProperty ( * ( noeud ),NUMGENES, Number<int> ( 1 ) );
      }
    }

    //if not at the root, there may be a long story from noeud to its father
    //so we need to take care of the potential losses between noeud and its father.
    if ( noeud->hasFather() && tree.getNode ( num )->hasFather() ) {
      int fatherNum = ( dynamic_cast<const Number<int> *> ( noeud->getFather()->getNodeProperty ( SPECIESID ) )->getValue() );
      if ( num!=fatherNum ) {
        int tempNum = num;
        Node * spTempNode = spNode;
        //we climb back in the species tree until we find the node whose id corresponds to the SPECIESID of the father of node noeud in the gene tree. Meanwhile, we take into account the gene losses, and the non-events.
        int tempNumBrother;
        while ( ( tempNum!=fatherNum ) && ( spTempNode->hasFather() ) && ( spTempNode->getFather()->getId() !=fatherNum ) ) {
          //There was a loss on branch leading to node brother to tempNum (tempNumBrother), and no duplication

          if ( spTempNode->getFather()->getSon ( 0 )->getId() ==tempNum ) {
            tempNumBrother = spTempNode->getFather()->getSon ( 1 )->getId();
          }
          else {
            tempNumBrother = spTempNode->getFather()->getSon ( 0 )->getId();
          }
          probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0 );
          num0lineages[tempNumBrother]++;
          branchNumbers[tempNumBrother]++;
          spTempNode = spTempNode->getFather();
          tempNum = spTempNode->getId();
          //No duplication on branch leading to node now named tempNum, and no loss either
          probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNum], lossProbabilities[tempNum], 1 );
          num0lineages[tempNum]++;
          branchNumbers[tempNum]++;
        }
        if ( ( * ( dynamic_cast < const std::string* > ( noeud->getFather()->getBranchProperty ( EVENT ) ) ) =="D" ) && ( num!=fatherNum ) ) {
          if ( spTempNode->getFather()->getSon ( 0 )->getId() ==tempNum ) {
            tempNumBrother = spTempNode->getFather()->getSon ( 1 )->getId();
          }
          else {
            tempNumBrother = spTempNode->getFather()->getSon ( 0 )->getId();
          }
          probaBranch+=computeLogBranchProbability ( duplicationProbabilities[tempNumBrother], lossProbabilities[tempNumBrother], 0 );
          num0lineages[tempNumBrother]++;
          branchNumbers[tempNumBrother]++;
        }
      }
    }

    changeBranchProperty ( * ( noeud ), EVENTSPROBA, Number<double> ( probaBranch ) );
    double lowlik = probaBranch+ ( dynamic_cast<const Number<double> *> ( son0->getNodeProperty ( LOWLIK ) )->getValue() ) + ( dynamic_cast<const Number<double> *> ( son1->getNodeProperty ( LOWLIK ) )->getValue() );
    changeNodeProperty ( * ( noeud ), LOWLIK, Number<double> ( lowlik ) );
  }
}



/**************************************************************************
 * This function writes a tree with a property whose format is integer.
 **************************************************************************/

std::string nodeToParenthesisWithIntNodeValues ( const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName ) throw ( NodeNotFoundException )
{
  if ( !tree.hasNode ( nodeId ) ) throw NodeNotFoundException ( "nodeToParenthesisWithIntNodeValues", nodeId );
  std::ostringstream s;
  if ( tree.isLeaf ( nodeId ) ) {
    s << tree.getNodeName ( nodeId ) << " "<< ( dynamic_cast<const Number<int> *> ( tree.getBranchProperty ( nodeId, propertyName ) )->getValue() );
  }
  else {
    s << "(";
    std::vector<int> sonsId = tree.getSonsId ( nodeId );
    s << nodeToParenthesisWithIntNodeValues ( tree, sonsId[0], bootstrap, propertyName );
    for ( unsigned int i = 1; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithIntNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
    s << ")";

    if ( bootstrap ) {
      if ( tree.hasBranchProperty ( nodeId, TreeTools::BOOTSTRAP ) )
        s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( nodeId, TreeTools::BOOTSTRAP ) )->getValue() );
    }
    else {
      if ( tree.hasBranchProperty ( nodeId, propertyName ) )
        s << ( dynamic_cast<const Number<int> *> ( tree.getBranchProperty ( nodeId, propertyName ) )->getValue() );
    }
  }
  if ( tree.hasDistanceToFather ( nodeId ) ) s << ":" << tree.getDistanceToFather ( nodeId );
  return s.str();
}




std::string treeToParenthesisWithIntNodeValues ( const Tree & tree, bool bootstrap, const std::string & propertyName )
{
  std::ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  std::vector<int> sonsId = tree.getSonsId ( rootId );
  if ( tree.isLeaf ( rootId ) ) {
    s << tree.getNodeName ( rootId );
    for ( unsigned int i = 0; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithIntNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
  }
  else {
    s << nodeToParenthesisWithIntNodeValues ( tree, sonsId[0], bootstrap, propertyName );
    for ( unsigned int i = 1; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithIntNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
  }
  s << ")";
  if ( bootstrap ) {
    if ( tree.hasBranchProperty ( rootId, TreeTools::BOOTSTRAP ) )
      s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( rootId, TreeTools::BOOTSTRAP ) )->getValue() );
  }
  else {
    if ( tree.hasBranchProperty ( rootId, propertyName ) )
      s << ( dynamic_cast<const Number<int> *> ( tree.getBranchProperty ( rootId, propertyName ) )->getValue() );
  }
  s << ";" << std::endl;
  return s.str();
}

/**************************************************************************
 * This function writes a tree with a property whose format is double.
 **************************************************************************/

std::string nodeToParenthesisWithDoubleNodeValues ( const Tree & tree, int nodeId, bool bootstrap, const std::string & propertyName ) throw ( NodeNotFoundException )
{
  if ( !tree.hasNode ( nodeId ) ) throw NodeNotFoundException ( "nodeToParenthesisWithDoubleNodeValues", nodeId );
  std::ostringstream s;
  if ( tree.isLeaf ( nodeId ) ) {
    s << tree.getNodeName ( nodeId ) << " "<< ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( nodeId, propertyName ) )->getValue() );
  }
  else {
    s << "(";
    std::vector<int> sonsId = tree.getSonsId ( nodeId );
    s << nodeToParenthesisWithDoubleNodeValues ( tree, sonsId[0], bootstrap, propertyName );
    for ( unsigned int i = 1; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithDoubleNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
    s << ")";

    if ( bootstrap ) {
      if ( tree.hasBranchProperty ( nodeId, TreeTools::BOOTSTRAP ) )
        s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( nodeId, TreeTools::BOOTSTRAP ) )->getValue() );
    }
    else {
      if ( tree.hasBranchProperty ( nodeId, propertyName ) )
        s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( nodeId, propertyName ) )->getValue() );
    }
  }
  if ( tree.hasDistanceToFather ( nodeId ) ) s << ":" << tree.getDistanceToFather ( nodeId );
  return s.str();
}




std::string treeToParenthesisWithDoubleNodeValues ( const Tree & tree, bool bootstrap, const std::string & propertyName )
{
  std::ostringstream s;
  s << "(";
  int rootId = tree.getRootId();
  std::vector<int> sonsId = tree.getSonsId ( rootId );
  if ( tree.isLeaf ( rootId ) ) {
    s << tree.getNodeName ( rootId );
    for ( unsigned int i = 0; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithDoubleNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
  }
  else {
    s << nodeToParenthesisWithDoubleNodeValues ( tree, sonsId[0], bootstrap, propertyName );
    for ( unsigned int i = 1; i < sonsId.size(); i++ ) {
      s << "," << nodeToParenthesisWithDoubleNodeValues ( tree, sonsId[i], bootstrap, propertyName );
    }
  }
  s << ")";
  if ( bootstrap ) {
    if ( tree.hasBranchProperty ( rootId, TreeTools::BOOTSTRAP ) )
      s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( rootId, TreeTools::BOOTSTRAP ) )->getValue() );
  }
  else {
    if ( tree.hasBranchProperty ( rootId, propertyName ) )
      s << ( dynamic_cast<const Number<double> *> ( tree.getBranchProperty ( rootId, propertyName ) )->getValue() );
  }
  s << ";" << std::endl;
  return s.str();
}



void setLossesAndDuplications ( TreeTemplate<Node> & tree,
                                std::vector <int> &lossNumbers,
                                std::vector <int> &duplicationNumbers )
{
  std::vector< int > nodesIds = tree.getNodesId ();

  for ( unsigned int i=0; i<nodesIds.size(); i++ ) {
    if ( tree.hasBranchProperty ( nodesIds[i], LOSSES ) ) {
      tree.getNode ( nodesIds[i] )->deleteBranchProperty ( LOSSES );
    }
    tree.getNode ( nodesIds[i] )->setBranchProperty ( LOSSES, BppString ( TextTools::toString ( ( lossNumbers[nodesIds[i]] ) ) ) );
    if ( tree.hasBranchProperty ( nodesIds[i], DUPLICATIONS ) ) {
      tree.getNode ( nodesIds[i] )->deleteBranchProperty ( DUPLICATIONS );
    }
    tree.getNode ( nodesIds[i] )->setBranchProperty ( DUPLICATIONS, BppString ( TextTools::toString ( ( duplicationNumbers[nodesIds[i]] ) ) ) );
  }
}

/**************************************************************************
 * This function sets the NUM0LINEAGES, NUM1LINEAGES and NUM2LINEAGES values on the species tree, given std::vectors
 **************************************************************************/
void assignNumLineagesOnSpeciesTree ( TreeTemplate<Node> & tree,
                                      std::vector <int> &num0Lineages,
                                      std::vector <int> &num1Lineages,
                                      std::vector <int> &num2Lineages )
{
  std::vector< int > nodesIds = tree.getNodesId ();
  std::string nums = "";
  for ( unsigned int i=0; i<nodesIds.size(); i++ ) {
    if ( tree.hasBranchProperty ( nodesIds[i], NUMLINEAGES ) ) {
      tree.getNode ( nodesIds[i] )->deleteBranchProperty ( NUMLINEAGES );
    }
    if ( tree.hasNodeProperty ( nodesIds[i], NUMLINEAGES ) ) {
      tree.getNode ( nodesIds[i] )->deleteNodeProperty ( NUMLINEAGES );
    }
    nums = TextTools::toString ( num0Lineages[nodesIds[i]] ) +"_"+TextTools::toString ( num1Lineages[nodesIds[i]] ) +"_"+TextTools::toString ( num2Lineages[nodesIds[i]] );
    tree.getNode ( nodesIds[i] )->setBranchProperty ( NUMLINEAGES, BppString ( nums ) );
    if ( tree.getNode ( nodesIds[i] )->isLeaf() ) {
      tree.getNode ( nodesIds[i] )->setName ( tree.getNode ( nodesIds[i] )->getName() + "_"+nums ) ;
    }
  }
}


/*****************************************************************************
 * This function returns the speciesID assigned to a leaf.
 *
 ****************************************************************************/


int assignSpeciesIdToLeaf ( Node * node,  const std::map<std::string, std::string > & seqSp,
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

void recoverLosses(Node *& node, int & a, const int & b, int & olda, int & a0,
                   const TreeTemplate<Node> & tree,
                   double & likelihoodCell,
                   const std::vector< double> & lossRates,
                   const std::vector< double> & duplicationRates,
                   int & dupData)
{

  // const Node* const node = tree.getNode(a);

  //std::cout <<"a "<<a <<std::endl;

//std::cout << TreeTemplateTools::treeToParenthesis(tree) <<std::endl;

  olda=a;

  Node* nodeA;

  if (node->hasFather()) {

    nodeA = node->getFather();

  }

  else {

    std::cout <<"Problem in recoverLosses, nodeA has no father"<<std::endl;

  }

  //a = node->getFather()->getId();


  a = nodeA->getId();

  // std::cout <<"a "<<a <<std::endl;

  node = nodeA;

  //  std::cout <<"a "<<a <<std::endl;

  //std::cout <<"here 10"<<std::endl;

  /* std::cout <<"nodeagetson0getId: "<<nodeA->getSon(0)->getId()<<std::endl;
   *
   s td::cout <<"nodeagetson1getId: "<<nodeA->getSon(1)->getId()<<std::endl;                                                                                                                                                                                    ***  *


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
   *
   l ikelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], dupData);                                                                                                                                                            ****

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

void recoverLossesWithDuplication ( const Node * nodeA,
                                    const int &a,
                                    const int &olda,
                                    const TreeTemplate<Node> & tree,
                                    double & likelihoodCell,
                                    const std::vector< double> & lossRates,
                                    const std::vector< double> & duplicationRates )
{
  //The loss has occured before a0
  // const Node * nodeA = tree.getNode(a);
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
  //UNDONE TEST 1009
  // if(!lostNode->isLeaf()) {
  likelihoodCell += computeLogBranchProbability ( duplicationRates[lostNode->getId()], lossRates[lostNode->getId()], 0 );
  // }
  // std::cout <<"recoverLossesWithDuplication loss on branch "<<lostNode->getId()<<std::endl;
  /* likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
   *     std::cout << "H 1 gene on branch "<<olda<<std::endl;*/
  //  std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
  return;
}

/*****************************************************************************
 * This function computes the lower conditional likelihood of a subtree and
 * assigns its summit node a species ID. Notations are influenced by
 * Zmasek and Eddy algorithm (2001).
 *
 ****************************************************************************/


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
                                                   bool atRoot )
{
  if ( rootLikelihood == 0.0 ) {
    int a, a0, olda;
    int b, b0, oldb;
    a = a0 = olda = son0SpId;
    b = b0 = oldb = son1SpId;
     /*std::cout <<"son0spid : "<<son0SpId<<std::endl;
       std::cout <<"son1spid : "<<son1SpId<<std::endl;
       std::cout <<"rootspid : "<<rootSpId<<std::endl;
     */

    /*  std::vector <int> ids = tree.getNodesId();
     *          for (int i =0; i<ids.size() ; i++) {
     *            std::cout <<"ids :"<<ids[i]<<std::endl;
  }
  */

    Node * temp0 = tree.getNode ( son0SpId );
    Node * temp1 = tree.getNode ( son1SpId );

    while ( a!=b ) { //There have been losses !
      if ( a>b ) {

        //  std::cout <<"before recoverLosses temp0id: "<<temp0.getId()<<std::endl;
        recoverLosses ( temp0, a, b, olda, a0, tree, rootLikelihood, lossRates, duplicationRates, son0DupData );
        /*  std::cout <<"HERE"<<std::endl;
         *                  std::cout <<"after recoverLosses temp0id: "<<temp0.getId()<<std::endl;*/
      }
      else {

        /*         std::cout <<"before recoverLosses temp1id: "<<temp1.getId()<<std::endl;
         *                        std::cout <<"HEHEHEEHEH\n\n"<<std::endl;*/
        recoverLosses ( temp1, b, a, oldb, b0, tree, rootLikelihood, lossRates, duplicationRates, son1DupData );
        // std::cout <<"after recoverLosses temp1id: "<<temp1.getId()<<std::endl;

      }
    }
    rootSpId = a;
    // std::cout <<"rootspid : "<<rootSpId<<std::endl;
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

        //  if (son0DupData>0) {
        /*  std::cout <<"C On branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
         *                  son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);*/
        // }
      }
      else { //The loss has occured before b0

        rootDupData += son0DupData+1;
        rootLikelihood-=computeLogBranchProbability ( duplicationRates[a0], lossRates[a0], son0DupData );
        recoverLossesWithDuplication ( temp1, b, oldb, tree, rootLikelihood, lossRates, duplicationRates );
        //  if (son1DupData>0) {
        /*  std::cout <<"D On branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
         *                  son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);*/
        //  }
      }
      //Counting the duplication(s) event(s)
      /*  if (a==a0) {
       *                rootDupData += son0DupData+1;
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
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbabilityAtRoot ( duplicationRates[a], lossRates[a], rootDupData );

      }
      else {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
        // std::cout<<"rootLikelihood5: "<< rootLikelihood<<std::endl;

      }
      //  }
    }
    else { //there was no duplication
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
      if ( atRoot ) {
        rootLikelihood += computeLogBranchProbabilityAtRoot ( duplicationRates[a], lossRates[a], rootDupData );

      }
      else {
        rootLikelihood += computeLogBranchProbability ( duplicationRates[a], lossRates[a], rootDupData );
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

double computeSubtreeLikelihoodPostorder ( TreeTemplate<Node> & spTree,
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
  int id=node->getId();
/*  std::cout <<  "computeSubtreeLikelihoodPostorder: id: "<< id << " ; " << TreeTemplateTools::treeToParenthesis(geneTree, true) <<std::endl;*/
  if ( node->isLeaf() ) {

    if ( likelihoodData[id][0]==0.0 ) {
      speciesIDs[id][0]=speciesIDs[id][1]=speciesIDs[id][2]=assignSpeciesIdToLeaf ( node, seqSp, spID );
      likelihoodData[id][0]=likelihoodData[id][1]=likelihoodData[id][2]=computeLogBranchProbability ( duplicationRates[speciesIDs[id][0]], lossRates[speciesIDs[id][0]], 1 );
      dupData[id][0] = dupData[id][1] = dupData[id][2] = 1;

      //  std::cout <<"leafLk id "<< id << " lk: "<<likelihoodData[id][0]<<std::endl;

    }
    /*  std::cout <<"at leaf "<<node->getName()<<std::endl;
     *          std::cout <<"lk "<<likelihoodData[id][0]<<std::endl;*/

    return ( likelihoodData[id][0] );
  }
  else {

    std::vector <Node *> sons = node->getSons();

    for ( unsigned int i = 0; i< sons.size(); i++ ) {
      computeSubtreeLikelihoodPostorder ( spTree, geneTree,
                                          sons[i], seqSp,
                                          spID, likelihoodData,
                                          lossRates, duplicationRates,
                                          speciesIDs, dupData );

    }

    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0, directionSon1;
    std::vector <Node *> neighbors = sons[0]->getNeighbors();
    for ( unsigned int i=0; i<neighbors.size(); i++ ) {
      if ( neighbors[i]==node ) {
        directionSon0 = i;
      }
    }
    neighbors = sons[1]->getNeighbors();
    for ( unsigned int i=0; i<neighbors.size(); i++ ) {
      if ( neighbors[i]==node ) {
        directionSon1 = i;
      }
    }

    /*  neighbors = node->getNeighbors();
     *          for (unsigned int i=0; i<neighbors.size(); i++) {
     *            if (neighbors[i]==node) {
     *              directionSon1 = i;
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
     *         std::cout << "son 0 lk "<<likelihoodData[idSon0][directionSon0]<< " directionSon0 "<<  directionSon0<<std::endl;
     *          std::cout << "son 1 lk "<<likelihoodData[idSon1][directionSon1]<<" directionSon1 "<<  directionSon1<<std::endl;
     *           std::cout <<"node ID "<<id<<"isRoot? "<<TreeTemplateTools::isRoot(*node)<<std::endl;*/

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
                                                TreeTemplateTools::isRoot ( *node ) );

    // std::cout <<"father lk "<< likelihoodData[id][0]<<std::endl;

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
void computeRootingLikelihood ( TreeTemplate<Node> & spTree,
                                Node * node,
                                std::vector <std::vector<double> > & likelihoodData,
                                const std::vector< double> & lossRates,
                                const std::vector < double> & duplicationRates,
                                std::vector <std::vector<int> > & speciesIDs,
                                std::vector <std::vector<int> > & dupData,
                                int sonNumber,
                                std::map <double, Node*> & LksToNodes )
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

  /*  std::cout <<"computeRootingLikelihood : nodeID "<<geneNodeId<<"; sonID "<<idNode1<<std::endl;
   *
   *      std::cout <<"likelihoodData[Father] "<< likelihoodData[idNode0][directionNode0]<<std::endl;
   *      std::cout <<"likelihoodData[Son] "<< likelihoodData[idNode1][directionNode1]<<std::endl;
   */
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
                                              dupData[idNode1][directionNode1], false );
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
                                              dupData[idSon1][directionSon1], true );
  // std::cout <<"LK FOUND "<<rootLikelihood<<std::endl;
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
                                        std::map <double, Node*> & LksToNodes )
{
  computeRootingLikelihood ( spTree, node,
                             likelihoodData, lossRates,
                             duplicationRates, speciesIDs,
                             dupData, sonNumber, LksToNodes );
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
  //  for (int i = 0; i< sons.size(); i++){
  {
    for ( unsigned int j =0; j<son->getNumberOfSons(); j++ ) {
      computeSubtreeLikelihoodPreorder ( spTree, geneTree,
                                         son, seqSp, spID,
                                         likelihoodData,
                                         lossRates, duplicationRates,
                                         speciesIDs, dupData, j, LksToNodes );
    }
  }
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




void recoverLossesAndLineages ( Node *& node, int & a, const int & b, int & olda, int & a0,
                                const TreeTemplate<Node> & tree,
                                int & dupData, std::vector<int> &num0lineages, std::vector<int> &num1lineages )
{
  // const Node* const node = tree.getNode(a);

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



void recoverLossesAndLineagesWithDuplication ( const Node * nodeA,
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
  // std::cout <<"G loss on branch "<<lostNode->getId()<<std::endl;
  /* likelihoodCell += computeLogBranchProbability(duplicationRates[olda], lossRates[olda], 1);
   *     std::cout << "H 1 gene on branch "<<olda<<std::endl;*/
  //  std::cout << "likelihoodCell "<< likelihoodCell <<std::endl;
  return;
}





/****************************************************************************/




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
      // if (son0DupData>0) {
      /*  std::cout <<"C On branch "<<a0<< "number of genes :"<<son0DupData<<std::endl;
       *             son0Likelihood += computeLogBranchProbability(duplicationRates[a0], lossRates[a0], son0DupData);*/
      //  }
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
      //  if (son1DupData>0) {
      /*  std::cout <<"D On branch "<<b0<< "number of genes :"<<son1DupData<<std::endl;
       *             son1Likelihood += computeLogBranchProbability(duplicationRates[b0], lossRates[b0], son1DupData);*/
      //  }
    }
    //Counting the duplication(s) event(s)
    /*  if (a==a0) {
     *         rootDupData += son0DupData+1;
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
    if ( rootDupData==1 ) {
      num1lineages[rootSpId]+=1;
    }
    else if ( rootDupData>=2 ) { //All branches with 2 or more lineages
      num2lineages[rootSpId]+=1;
    }
    //  }
  }
  else { //there was no duplication
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
        computeNumbersOfLineagesFromRoot ( spTree, geneTree, sons[i],
                                           seqSp, spID, num0lineages,
                                           num1lineages, num2lineages,
                                           speciesIDs, dupData,
                                           branchesWithDuplications );
      }
    }
    int idSon0 = sons[0]->getId();
    int idSon1 = sons[1]->getId();
    unsigned int directionSon0, directionSon1;
    std::vector <Node *> neighbors = sons[0]->getNeighbors();
    for ( unsigned int i=0; i<neighbors.size(); i++ ) {
      if ( neighbors[i]==node ) {
        directionSon0 = i;
      }
    }
    neighbors = sons[1]->getNeighbors();
    for ( unsigned int i=0; i<neighbors.size(); i++ ) {
      if ( neighbors[i]==node ) {
        directionSon1 = i;
      }
    }

    computeNumbersOfLineagesInASubtree ( *spTree, sons,
                                         speciesIDs[id][0],
                                         speciesIDs[idSon0][directionSon0],
                                         speciesIDs[idSon1][directionSon1],
                                         dupData[id][0], dupData[idSon0][directionSon0],
                                         dupData[idSon1][directionSon1],
                                         TreeTemplateTools::isRoot ( *node ),
                                         num0lineages, num1lineages,
                                         num2lineages, branchesWithDuplications );
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
                                const bool fillTables )
{

/*  std::cout<< "findMLReconciliationDR "<<std::endl;
  VectorTools::print(lossRates);
    VectorTools::print(duplicationRates);
std::cout << TreeTemplateTools::treeToParenthesis(*spTree, true)<<std::endl;
std::cout << TreeTemplateTools::treeToParenthesis(*geneTree, true)<<std::endl;*/

  if ( !geneTree->isRooted() )
  {
    std::cout << TreeTemplateTools::treeToParenthesis ( *geneTree, true ) <<std::endl;
    std::cout <<"!!!!!!gene tree is not rooted in findMLReconciliationDR !!!!!!"<<std::endl;
    exit ( -1 );
  }
  if ( !spTree->isRooted() )
  {
    std::cout << TreeTemplateTools::treeToParenthesis ( *spTree, true ) <<std::endl;
    std::cout <<"!!!!!!species tree is not rooted in findMLReconciliationDR !!!!!!"<<std::endl;
    exit ( -1 );
  }
  std::vector <double> nodeData ( 3, 0.0 );
  std::vector <std::vector<double> > likelihoodData ( geneTree->getNumberOfNodes(), nodeData );

  std::vector <int> nodeSpId ( 3, 0 );
  std::vector <std::vector<int> > speciesIDs ( geneTree->getNumberOfNodes(), nodeSpId );
  std::vector <std::vector<int> > dupData = speciesIDs;

  double initialLikelihood;
  //This std::map keeps rootings likelihoods. The key is the likelihood value, and the value is the node to put as outgroup.
  std::map <double, Node*> LksToNodes;

  Node * geneRoot = geneTree->getRootNode();


  {
    {
      //std::cout << "findMLReconciliationDR geneTree: " << TreeTemplateTools::treeToParenthesis ( *geneTree, true ) <<std::endl;

      initialLikelihood = computeSubtreeLikelihoodPostorder ( *spTree, *geneTree,
                                                              geneRoot, seqSp, spID,
                                                              likelihoodData, lossRates,
                                                              duplicationRates, speciesIDs, dupData );
    }
  }
 // std::cout << "computeSubtreeLikelihoodPostorder took: "<< "OMP DEACTIVATED"<<std::endl;


  std::vector <Node *> sons = geneRoot->getSons();

  if ( sons.size() !=2 ) {
    std::cerr <<"Error: "<<sons.size() << "sons at the root!"<<std::endl;
  }
  LksToNodes[initialLikelihood]=sons[0];
  //We fill the likelihood and species ID data for the root node.
  //We use "directions" 1 and 2 and leave "direction" 0 empty for coherence
  //with other nodes.
  likelihoodData[geneRoot->getId()][1] = likelihoodData[geneRoot->getSon ( 1 )->getId()][0];
  likelihoodData[geneRoot->getId()][2] = likelihoodData[geneRoot->getSon ( 0 )->getId()][0];
  speciesIDs[geneRoot->getId()][1] = speciesIDs[geneRoot->getSon ( 1 )->getId()][0];
  speciesIDs[geneRoot->getId()][2] = speciesIDs[geneRoot->getSon ( 0 )->getId()][0];
  dupData[geneRoot->getId()][1] = dupData[geneRoot->getSon ( 1 )->getId()][0];
  dupData[geneRoot->getId()][2] = dupData[geneRoot->getSon ( 0 )->getId()][0];


  {
    {
      for ( unsigned int i = 0; i< sons.size(); i++ ) {
        for ( unsigned int j =0; j<sons[i]->getNumberOfSons(); j++ ) {
          computeSubtreeLikelihoodPreorder ( *spTree, *geneTree,
                                             sons[i], seqSp, spID,
                                             likelihoodData,
                                             lossRates, duplicationRates,
                                             speciesIDs, dupData, j, LksToNodes );
        }
      }
    }
  }
 // std::cout << "computeSubtreeLikelihoodPreorder took: "<< "omp deactivated" <<std::endl;




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


    {
      {
        computeNumbersOfLineagesFromRoot ( spTree, tree,
                                           tree->getRootNode(),
                                           seqSp, spID,
                                           num0lineages, num1lineages,
                                           num2lineages, speciesIDs,
                                           dupData, nodesToTryInNNISearch );
      }
    }
  //  std::cout << "computeNumbersOfLineagesFromRoot took: "<< "omp deactivated" <<std::endl;

    delete tree;

  }

  //We return the best likelihood
  MLindex = LksToNodes.rbegin()->second->getId();

  return LksToNodes.rbegin()->first;



}


/**************************************************************************/
//We compute the average loss proportion observed among all species tree branches.
double computeAverageLossProportion ( std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages )
{
  double totNum0 = ( double ) VectorTools::sum ( num0lineages );
  double totNum1 = ( double ) VectorTools::sum ( num1lineages );
  double totNum2= ( double ) VectorTools::sum ( num2lineages );
  if ( totNum0+totNum1+totNum2==0 ) { //This happens at the beginning of the program
    resetLineageCounts ( num0lineages, num1lineages, num2lineages );
  }
  totNum0 = ( double ) VectorTools::sum ( num0lineages );
  totNum1 = ( double ) VectorTools::sum ( num1lineages );
  totNum2= ( double ) VectorTools::sum ( num2lineages );
  double prop = totNum0/ ( totNum0+totNum1+totNum2 );
  return prop;
}


void extractSubVectorsWithInternalLineages ( std::vector <int> & num0lineages,
                                             std::vector <int> & num1lineages,
                                             std::vector <int> & num2lineages,
                                             std::map <std::string, int> & genomeMissing,
                                             TreeTemplate<Node> & tree,
                                             std::vector <int> & num0,
                                             std::vector <int> & num1,
                                             std::vector <int> & num2 )
{
  std::vector <int> branchesToDiscard = tree.getLeavesId();
  //We also discard the root branch, which by definition does not have any loss
  branchesToDiscard.push_back ( 0 );

  sort ( branchesToDiscard.begin(), branchesToDiscard.end() );

  num0 = num0lineages;
  num1 = num1lineages;
  num2 = num2lineages;


  for ( std::vector<int>::reverse_iterator it = branchesToDiscard.rbegin(); it != branchesToDiscard.rend(); it++ ) {
    num0.erase ( num0.begin() +*it );
    num1.erase ( num1.begin() +*it );
    num2.erase ( num2.begin() +*it );
  }
  return;

}




//We compute the average loss proportion observed among species tree branches that
//are not found in "genomeMissing" to have missing data.

void extractSubVectorsWithCompletelySequencedLineages ( const std::vector <int> & num0lineages,
                                                        const std::vector <int> & num1lineages,
                                                        const std::vector <int> & num2lineages,
                                                        const std::map <std::string, int> & genomeMissing,
                                                        const TreeTemplate<Node> & tree,
                                                        std::vector <int> & num0,
                                                        std::vector <int> & num1,
                                                        std::vector <int> & num2 )
{
  std::vector <int> branchesToDiscard;
  for ( std::map<std::string, int >::const_iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++ ) {
    if ( it->second != 0 ) {
      branchesToDiscard.push_back ( tree.getLeafId ( it->first ) );
    }
  }
  //We also discard the root branch, which by definition does not have any loss
  branchesToDiscard.push_back ( 0 );

  sort ( branchesToDiscard.begin(), branchesToDiscard.end() );


 /* std::cout << "num0lineages.size(): " << num0lineages.size() << std::endl;
  std::cout <<  branchesToDiscard.size() << branchesToDiscard.size() << std::endl;*/
  num0 = num0lineages;
  num1 = num1lineages;
  num2 = num2lineages;

  for ( std::vector<int>::reverse_iterator it = branchesToDiscard.rbegin(); it != branchesToDiscard.rend(); it++ ) {
    num0.erase ( num0.begin() +*it );
    num1.erase ( num1.begin() +*it );
    num2.erase ( num2.begin() +*it );
  }
  return;

}




double computeAverageLossProportionOnCompletelySequencedLineages ( const std::vector <int> & num0lineages,
                                                                   const std::vector <int> & num1lineages,
                                                                   const std::vector <int> & num2lineages,
                                                                   const std::map <std::string, int> & genomeMissing,
                                                                   const TreeTemplate<Node> & tree )
{

  std::vector  <int> num0;
  std::vector  <int> num1;
  std::vector  <int> num2;

  extractSubVectorsWithCompletelySequencedLineages ( num0lineages,
                                                     num1lineages,
                                                     num2lineages,
                                                     genomeMissing,
                                                     tree,
                                                     num0,
                                                     num1,
                                                     num2 );

  double totNum0 = ( double ) VectorTools::sum ( num0 );
  double totNum1 = ( double ) VectorTools::sum ( num1 );
  double totNum2= ( double ) VectorTools::sum ( num2 );
  if ( totNum0+totNum1+totNum2==0 ) { //This happens at the beginning of the program
    totNum0 = 267;
    totNum1 = 723;
    totNum2 = 10;
  }
  double prop = totNum0/ ( totNum0+totNum1+totNum2 );
  return prop;
}


/**************************************************************************/
//We increase num0lineages in some external branches to account for the low sequence coverage of some genomes.
//We also set num0Lineages in the root branch.
void alterLineageCountsWithCoverages ( std::vector <int> & num0lineages,
                                       std::vector <int> & num1lineages,
                                       std::vector <int> & num2lineages,
                                       const std::map <std::string, int> & genomeMissing,
                                       const TreeTemplate<Node> & tree, bool average )
{
  std::vector  <int> num0;
  std::vector  <int> num1;
  std::vector  <int> num2;
  extractSubVectorsWithCompletelySequencedLineages ( num0lineages,
                                                     num1lineages,
                                                     num2lineages,
                                                     genomeMissing,
                                                     tree,
                                                     num0, num1, num2 );

  double avg0d = VectorTools::mean<int, double> ( num0 );
  double avg1d = VectorTools::mean<int, double> ( num1 );
  double avg2d = VectorTools::mean<int, double> ( num2 );

  if ( avg0d+avg1d+avg2d==0 ) { //This shouldn't happen but if it does...
    std::cerr<<"WARNING: No event on average on all branches. Setting default values."<<std::endl;
    avg0d = 267;
    avg1d = 723;
    avg2d = 10;
  }

  double propLoss = avg0d/ ( avg0d+avg1d+avg2d );
  if ( propLoss > 0.99 ) {
    propLoss = 0.99;
  }

  int avg0 = ( int ) avg0d;
  int avg1 = ( int ) avg1d;
  int avg2 = ( int ) avg2d;

  if ( average ) {
    for ( unsigned int i =0; i<num0lineages.size(); i++ ) {
      num0lineages[i]=avg0;
      num1lineages[i]=avg1;
      num2lineages[i]=avg2;
    }
  }
  else {
    //At branch 0 and branches surrounding the root, by definition, we never count cases where there has been a loss, so we set it to an average value.
    //We also put all num0, num1, and num2 values at the average values computed above.
    for ( unsigned int i =0; i<num0lineages.size(); i++ ) {
      if ( num0lineages[i] == 0 ) {
        num0lineages[i]=avg0;
        num1lineages[i]=avg1;
        num2lineages[i]=avg2;
      }
    }
  }

  //  Now we apply corrections for poorly sequenced genomes
  for ( std::map<std::string, int >::const_iterator it = genomeMissing.begin(); it != genomeMissing.end(); it++ ) {
    int id = tree.getLeafId ( it->first );
    int percent = it->second;
    if ( percent>0 ) {
      double percentd = ( double ) percent/100.0;
      double totLoss = propLoss +percentd;
      if ( totLoss>0.99 ) {
        totLoss=0.99;
      }

      while ( ( num1lineages[id] < 10 ) || ( num2lineages[id] < 10 ) ) {
        num1lineages[id] = num1lineages[id]*100;
        num2lineages[id] = num2lineages[id]*100;
        if ( num1lineages[id] < 1 ) {
          num1lineages[id] =1;
        }
        if ( num2lineages[id] < 1 ) {
          num2lineages[id] =1;
        }
      }
      num0lineages[id] = ( int ) ( totLoss * ( ( double ) num1lineages[id]+ ( double ) num2lineages[id] ) / ( 1.0-totLoss ) );
    }
  }
}

/**************************************************************************/
//We set all duplication and loss values to values obtained on a former dataset (loss rate 0.313, duplication rate 0.0159, which corresponds to num0lineages=267, num1lineages=723, num2lineages=10), and increase num0lineages in some external branches to account for the low sequence coverage of some genomes.
void alterLineageCountsWithCoveragesInitially ( std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages, const std::map <std::string, int> & genomeMissing, const TreeTemplate<Node> & tree ) {
  //Setting initial values for the std::vectors
  resetLineageCounts ( num0lineages, num1lineages, num2lineages );

  alterLineageCountsWithCoverages ( num0lineages, num1lineages, num2lineages, genomeMissing, tree, false );
}


void resetLineageCounts ( std::vector <int> & num0lineages, std::vector <int> & num1lineages, std::vector <int> & num2lineages ) {
  for ( unsigned int i =0; i<num0lineages.size(); i++ ) {
    num0lineages[i]=267;
    num1lineages[i]=723;
    num2lineages[i]=1;
    /* num0lineages[i]=267*(i+1);
     *         num1lineages[i]=723;
     *         num2lineages[i]=10;*/
  }
}

/**************************************************************************/
// This function removes properties associated to a tree.



void deleteSubtreeProperties ( Node &node ) {
  for ( unsigned int i = 0; i < node.getNumberOfSons(); i++ ) {
    Node * son = node.getSon ( i );
    deleteSubtreeProperties ( * son );
    son->deleteBranchProperties();
    son->deleteNodeProperties();
  }
}



void deleteTreeProperties ( TreeTemplate<Node> & tree ) {
  deleteSubtreeProperties ( * ( tree.getRootNode() ) );
}




/**************************************************************************/


std::map <std::string, int> computeSpeciesNamesToIdsMap ( TreeTemplate<Node> & tree ) {
  std::map <std::string, int> spId;
  std::vector <Node *> nodes = tree.getNodes();
  for ( unsigned int i = 0; i< nodes.size() ; i++ ) {
    if ( nodes[i]->isLeaf() ) {
      spId[nodes[i]->getName()] = nodes[i]->getId();
    }
  }
  return spId;
}

/**************************************************************************/



void computeDuplicationAndLossRatesForTheSpeciesTree ( std::string &branchProbaOptimization,
                                                       std::vector <int> & num0Lineages,
                                                       std::vector <int> & num1Lineages,
                                                       std::vector <int> & num2Lineages,
                                                       std::vector<double> & lossExpectedNumbers,
                                                       std::vector<double> & duplicationExpectedNumbers,
                                                       std::map <std::string, int> & genomeMissing,
                                                       TreeTemplate<Node> & tree )
{
  //outputting trees with branch lengths in numbers of events, before correction.
  computeDuplicationAndLossProbabilitiesForAllBranches ( num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers );
  std::cout << "Species tree with expected numbers of duplications as branch lengths:"<<std::endl;
  for ( unsigned int i =0; i<num0Lineages.size() ; i++ ) {
    tree.getNode ( i )->setBranchProperty ( "DUPLICATIONS", Number<double> ( duplicationExpectedNumbers[i] ) );
    if ( tree.getNode ( i )->hasFather() ) {
      tree.getNode ( i )->setDistanceToFather ( duplicationExpectedNumbers[i] );
    }
  }
  std::cout << treeToParenthesisWithDoubleNodeValues ( tree, false, "DUPLICATIONS" ) <<std::endl;
  std::cout << "Species tree with expected numbers of losses as branch lengths:"<<std::endl;
  for ( unsigned int i =0; i<num0Lineages.size() ; i++ ) {
    tree.getNode ( i )->setBranchProperty ( "LOSSES", Number<double> ( lossExpectedNumbers[i] ) );
    if ( tree.getNode ( i )->hasFather() ) {
      tree.getNode ( i )->setDistanceToFather ( lossExpectedNumbers[i] );
    }
  }
  std::cout << treeToParenthesisWithDoubleNodeValues ( tree, false, "LOSSES" ) <<std::endl;

  //Doing the correction:
  if ( branchProbaOptimization=="average" ) {
    std::cout << "Averaging branch-wise estimates."<<std::endl;
    alterLineageCountsWithCoverages ( num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree, true );
  }
  else {
    alterLineageCountsWithCoverages ( num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree, false );
  }

  computeDuplicationAndLossProbabilitiesForAllBranches ( num0Lineages, num1Lineages, num2Lineages, lossExpectedNumbers, duplicationExpectedNumbers );

}


/**************************************************************************/
void computeDuplicationAndLossRatesForTheSpeciesTreeInitially ( std::string &branchProbaOptimization,
                                                                std::vector <int> & num0Lineages,
                                                                std::vector <int> & num1Lineages,
                                                                std::vector <int> & num2Lineages,
                                                                std::vector<double> & lossProbabilities,
                                                                std::vector<double> & duplicationProbabilities,
                                                                const std::map <std::string, int> & genomeMissing,
                                                                const TreeTemplate<Node> & tree )
{
//  std::cout <<"Computing Initial DL Expected Numbers"<<std::endl;
  alterLineageCountsWithCoveragesInitially ( num0Lineages, num1Lineages, num2Lineages, genomeMissing, tree );
  computeDuplicationAndLossProbabilitiesForAllBranches ( num0Lineages, num1Lineages, num2Lineages, lossProbabilities, duplicationProbabilities );
}


/**************************************************************************/
/*Removes a leaf in a tree*/
void removeLeaf ( TreeTemplate<Node> & tree, std::string toRemove ) {
  // std::cout <<"in removeLeaf 1"<<std::endl;
  Node * NToRemove  = tree.getNode ( toRemove );
  if ( !NToRemove->hasFather() ) {
    // std::cout <<"node is root !!!!!"<<std::endl;
    Node * father = NToRemove->getSon ( 0 );
    Node * son1 = father->getSon ( 1 );
    father->removeFather();
    tree.newOutGroup ( son1->getId() );
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
  // std::cout <<"in removeLeaf 2"<<std::endl;
  Node * father = NToRemove->getFather();
  // std::cout <<"in removeLeaf 3"<<std::endl;
  Node * brother;
  for ( unsigned int i=0; i<father->getNumberOfSons(); i++ )
    if ( father->getSon ( i ) !=NToRemove ) {
      brother=father->getSon ( i );
      break;
    }
    // std::cout <<"in removeLeaf 4"<<std::endl;
    double distBro;
  try {
    distBro = brother->getDistanceToFather();
  }
  catch ( std::exception ) {
    distBro = 0.000001;
  }
  // std::cout <<"in removeLeaf 5"<<std::endl;
  if ( !father->hasFather() ) {
    // std::cout <<"father is root !!!!!"<<std::endl;
    // brother->removeFather();
    tree.rootAt ( brother->getId() );
    // std::cout <<"After rootAt"<<std::endl;
    //tree.newOutGroup(brother->getSon(0)->getId());
    for ( unsigned int i = 0; i<brother->getNumberOfSons(); i++ ) {
      if ( brother->getSon ( i ) ==NToRemove ) {
        brother->removeSon ( i );
        break;
      }
    }
    tree.newOutGroup ( brother->getSon ( 0 )->getId() );
    delete NToRemove;
    tree.resetNodesId();
    return;
  }
  double distFa;
  try {
    distFa = father->getDistanceToFather();
  }
  catch ( std::exception ) {
    distFa = 0.000001;
  }
  // std::cout <<"in removeLeaf 6"<<std::endl;
  Node * grandFather = father->getFather();
  grandFather->addSon ( brother );
  brother->setDistanceToFather ( distBro+distFa );
  for ( unsigned int i=0; i<grandFather->getNumberOfSons(); i++ )
    if ( grandFather->getSon ( i ) ==father ) {
      grandFather->removeSon ( i );
      break;
    }
    /*if (!grandFather->hasFather()) {
     *      tree.newOutGroup(grandFather->getSon(0)->getId());
}*/
    tree.resetNodesId();
    /*  delete NToRemove;
     *      delete father;*/
}

/**************************************************************************/

/******************************************************************************/
//This function cleans a vector of options read from a file.
/******************************************************************************/


void cleanVectorOfOptions ( std::vector<std::string> & listOptions, bool sizeConstraint ) {
  int maxStrSize=0;
  std::vector <int> toRemove;
  int i=0;
  for ( std::vector<std::string >::iterator it = listOptions.begin(); it != listOptions.end(); it++ ) {
    std::string arg = removeComments ( *it, std::string ( "#" ), std::string ( "\n" ) ); //remove shell comments.
    arg = removeComments ( arg, std::string ( "//" ), std::string ( "\n" ) ); //remove C simple comments.
    arg = removeComments ( arg, std::string ( "/*" ), std::string ( "*/" ) ); //remove C multiple comments.
    arg = TextTools::removeWhiteSpaces ( arg );
    *it = arg;
    int temp=it->length();
    if ( temp==0 ) {
      toRemove.push_back ( i );
    }
    else {
      if ( temp>maxStrSize ) {
        maxStrSize=temp;
      };
    }
    i++;
  }
  for ( std::vector<int>::reverse_iterator it = toRemove.rbegin(); it != toRemove.rend(); it++ ) {
    listOptions.erase ( listOptions.begin() +*it );
  }
  if ( sizeConstraint ) {
    if ( maxStrSize>MAXFILENAMESIZE ) {
      std::cout << "\nNames are too long, please abbreviate ! File names (including the path) need to be < "<< MAXFILENAMESIZE <<" letters long.\n"<<std::endl;
      exit ( -1 );
    }
  }
  return;
}

/**************************************************************************/


/**************************************************************************
 * This function distributes gene families so as to homogenize computational load between clients.
 **************************************************************************/

bool sortMaxFunction ( std::pair <std::string, double> i, std::pair <std::string, double> j ) {
  if ( i.second > j.second ) {
    return ( true );
  }
  else {
    return ( false );
  }
}

bool sortMinFunction ( std::pair <std::vector<std::string>, double> i, std::pair <std::vector<std::string>, double> j ) {
  if ( i.second < j.second ) {
    return ( true );
  }
  else {
    return ( false );
  }
}




void generateListOfOptionsPerClient ( std::vector <std::string> listOptions, int size, std::vector <std::vector<std::string> > &listOfOptionsPerClient, std::vector <unsigned int> &numbersOfGenesPerClient ) {
  //Here, two alternatives: either we do have information regarding the gene family sizes, or we don't.
  //Is there size information? == Is there a ":" in the first line?
  if ( TextTools::hasSubstring ( listOptions[0],":" ) ) {
    std::vector <std::pair <std::string, double> > elements;
    for ( unsigned int i = 0; i<listOptions.size() ; i++ ) {
      StringTokenizer st1 ( listOptions[i], ":", true );
      elements.push_back ( std::pair <std::string, double> ( st1.getToken ( 0 ), TextTools::toDouble ( st1.getToken ( 1 ) ) ) );
    }

    //Now sort the gene families by their size, in descending order
    sort ( elements.begin(), elements.end(), sortMaxFunction );
    //Now we assign gene families to nodes.
    //We start with big families, and then go in decreasing order.
    //First we assign the first families
    std::vector <std::pair <std::vector<std::string>, double> > listOfOptionsAndTotSizePerClient;

    std::vector<std::string> temp2;
    std::pair <std::vector<std::string>, double> temp3;
    unsigned int j = 0;

    for ( int i = 0; i<size-1 ; i++ ) {
      listOfOptionsAndTotSizePerClient.push_back ( temp3 );
      listOfOptionsAndTotSizePerClient[i].first.push_back ( elements[j].first );
      listOfOptionsAndTotSizePerClient[i].second = elements[j].second;
      j = j+1;
    }

    while ( j<listOptions.size() ) {
      //We sort listOfOptionsAndTotSizePerClient in increasing function
      sort ( listOfOptionsAndTotSizePerClient.begin(), listOfOptionsAndTotSizePerClient.end(), sortMinFunction );
      listOfOptionsAndTotSizePerClient[0].first.push_back ( elements[j].first );
      listOfOptionsAndTotSizePerClient[0].second = listOfOptionsAndTotSizePerClient[0].second + elements[j].second;
      j= j+1;
    }

    //Now all gene families must have been assigned to nodes.
    // std::vector <std::vector<std::string> > listOfOptionsPerClient;
    listOfOptionsPerClient.push_back ( temp2 );
    listOfOptionsPerClient[0].push_back ( std::string ( "####" ) ); //For the root node
    numbersOfGenesPerClient.push_back ( 0 );

    //We print the result of the assignment and fill listOfOptionsPerClient:
    for ( int i = 0; i<size-1 ; i++ ) {
      std::cout <<"Client "<<i<<" is in charge of "<< listOfOptionsAndTotSizePerClient[i].first.size() <<" gene families; Total Weight : "<<listOfOptionsAndTotSizePerClient[i].second<<std::endl;
      listOfOptionsPerClient.push_back ( listOfOptionsAndTotSizePerClient[i].first );
      numbersOfGenesPerClient.push_back ( listOfOptionsAndTotSizePerClient[i].first.size() );
    }
    return ;
    //generateListOfOptionsPerClient(listOptions, size, listOfOptionsPerClient, numbersOfGenesPerClient);
  }
  else {
    int numberOfGenesPerClient = ( int ) ( listOptions.size() ) / ( size -1 );
    int reste = ( listOptions.size() ) % ( size -1 );
    std::cout <<"Number of genes per client : "<<numberOfGenesPerClient<< " and extra "<<reste<<std::endl;
    for ( int i = 0 ; i< size-1 ; i++ ) {
      if ( i<reste ) {
        numbersOfGenesPerClient.push_back ( numberOfGenesPerClient+1 );
      }
      else {
        numbersOfGenesPerClient.push_back ( numberOfGenesPerClient );
      }
    }
    std::vector<unsigned int>::iterator it2 = numbersOfGenesPerClient.begin();
    numbersOfGenesPerClient.insert ( it2, int ( 0 ) ); //For the server, we insert a "dumb" option file at the beginning of the std::vector, so only clients compute the reconciliation
    int currentFile = 0;
    std::vector<std::string> temp2;
    for ( int i = 0 ; i< size ; i++ ) {
      listOfOptionsPerClient.push_back ( temp2 );
      if ( i == 0 ) {
        listOfOptionsPerClient[i].push_back ( std::string ( "####" ) );
      }
      else {
        for ( unsigned int j = 0 ; j < numbersOfGenesPerClient[i] ; j++ ) {
          listOfOptionsPerClient[i].push_back ( listOptions[currentFile] );
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

std::string removeComments (
  const std::string & s,
  const std::string & begin,
  const std::string & end )
{
  std::string r = s;
  std::string::size_type last = 0;
  do {
    std::string::size_type first = r.find ( begin, last );
    if ( first == std::string::npos ) return r; //No shell comment.
    //else:
    last = r.find ( end, first );
    if ( last == std::string::npos ) {
      r.erase ( r.begin() + first, r.end() );
    }
    else {
      r.erase ( r.begin() + first, r.begin() + last );
    }
  }
  while ( last != std::string::npos );
  return r;
}




/**************************************************************************
 * Annotate gene tree with duplication events. Adds D=Y to node properties where
 * there was a duplication, and D=N where there was no duplication.
 * Also adds species node id for each node of the gene tree, using :S=spId.
 *************************************************************************/

void annotateGeneTreeWithDuplicationEvents ( TreeTemplate<Node> & spTree,
                                             TreeTemplate<Node> & geneTree,
                                             Node * node,
                                             const std::map<std::string, std::string > seqSp,
                                             const std::map<std::string, int > spID )
{
  if ( node->isLeaf() ) {
    node->setNodeProperty ( "S", BppString ( TextTools::toString ( assignSpeciesIdToLeaf ( node, seqSp, spID ) ) ) );
    node->setBranchProperty ( "Ev", BppString ( "S" ) );
    return;
  }
  else {
    std::vector <Node *> sons = node->getSons();
    for ( unsigned int i = 0; i< sons.size(); i++ ) {
      annotateGeneTreeWithDuplicationEvents ( spTree, geneTree, sons[i], seqSp, spID );
    }

    int a = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[0]->getNodeProperty ( "S" ) ) )->toSTL() );
    int b = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[1]->getNodeProperty ( "S" ) ) )->toSTL() );

    int aold = a;
    int bold = b;
    int lossA = 0;
    int lossB = 0;

    while ( a!=b ) {
      if ( a>b ) {
        a = spTree.getNode ( a )->getFather()->getId();
        lossA = lossA +1;
      }
      else {
        b = spTree.getNode ( b )->getFather()->getId();
        lossB = lossB + 1;
      }
    }
    sons[0]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossA ) ) );
    sons[1]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossB ) ) );
    node->setNodeProperty ( "S", BppString ( TextTools::toString ( a ) ) );
    if ( ( a == aold ) || ( a == bold ) ) {
      node->setBranchProperty ( "Ev", BppString ( "D" ) );
    }
    else {
      node->setBranchProperty ( "Ev", BppString ( "S" ) );
    }

    return;
  }
}


/**************************************************************************
 * Annotate gene tree with duplication events. Adds D=Y to node properties where
 * there was a duplication, and D=N where there was no duplication.
 * Also adds species node id for each node of the gene tree, using :S=spId.
 *************************************************************************/

void annotateGeneTreeWithScoredDuplicationEvents ( TreeTemplate<Node> & spTree,
                                                   TreeTemplate<Node> & geneTree,
                                                   Node * node,
                                                   const std::map<std::string, std::string > seqSp,
                                                   const std::map<std::string, int > spID )
{
  //    std::cout << "annotateGeneTreeWithScoredDuplicationEvents "<<std::endl;
  if ( node->isLeaf() ) {
    //     std::cout << "annotateGeneTreeWithScoredDuplicationEvents Leaf"<<std::endl;

    int sp = assignSpeciesIdToLeaf ( node, seqSp, spID );
    node->setNodeProperty ( "S", BppString ( TextTools::toString ( sp ) ) );
    node->setNodeProperty ( "spPresentInSubtree", BppVector<int> ( 1, sp ) );
    node->setBranchProperty ( "Ev", BppString ( "S" ) );
    //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents Leaf 2"<<std::endl;

    return;
  }
  else {
    //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf "<<std::endl;

    std::vector <Node *> sons = node->getSons();
    for ( unsigned int i = 0; i< sons.size(); i++ ) {
      annotateGeneTreeWithScoredDuplicationEvents ( spTree, geneTree, sons[i], seqSp, spID );
    }

    //std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 2"<<std::endl;

    int a = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[0]->getNodeProperty ( "S" ) ) )->toSTL() );
    int b = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[1]->getNodeProperty ( "S" ) ) )->toSTL() );
    //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 3"<<std::endl;

    int aold = a;
    int bold = b;
    int lossA = 0;
    int lossB = 0;

    while ( a!=b ) {
      if ( a>b ) {
        a = spTree.getNode ( a )->getFather()->getId();
        lossA = lossA +1;
      }
      else {
        b = spTree.getNode ( b )->getFather()->getId();
        lossB = lossB + 1;
      }
    }
    //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 4"<<std::endl;

    sons[0]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossA ) ) );
    sons[1]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossB ) ) );
    node->setNodeProperty ( "S", BppString ( TextTools::toString ( a ) ) );
    //VectorTools::print( (dynamic_cast<const BppVector<int> *>(sons[0]->getNodeProperty("spPresentInSubtree")))->toSTL() );
    //VectorTools::print( (dynamic_cast<const BppVector<int> *>(sons[1]->getNodeProperty("spPresentInSubtree")))->toSTL() );

    std::vector<int> vec = VectorTools::vectorUnion ( ( dynamic_cast<const BppVector<int> *> ( sons[0]->getNodeProperty ( "spPresentInSubtree" ) ) )->toSTL(), ( dynamic_cast<const BppVector<int> *> ( sons[1]->getNodeProperty ( "spPresentInSubtree" ) ) )->toSTL() );
    node->setNodeProperty ( "spPresentInSubtree", BppVector<int> ( vec.begin(), vec.end() ) );
    //    std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 5"<<std::endl;

    if ( ( a == aold ) || ( a == bold ) ) {
      node->setBranchProperty ( "Ev", BppString ( "D" ) );
      node->setNodeProperty ( "Ev", BppString ( "D" ) );
      //Need to compute the intersection (and its size) and divide by the size of the union (computed above)
      //to compute the score of the duplication event.
      std::vector<int> vec2 = VectorTools::vectorIntersection ( ( dynamic_cast<const BppVector<int> *> ( sons[0]->getNodeProperty ( "spPresentInSubtree" ) ) )->toSTL(), ( dynamic_cast<const BppVector<int> *> ( sons[1]->getNodeProperty ( "spPresentInSubtree" ) ) )->toSTL() );
      node->setNodeProperty ( "Score", BppString ( TextTools::toString ( vec2.size() / vec.size() ) ) );
      node->setBranchProperty ( "Score", BppString ( TextTools::toString ( vec2.size() / vec.size() ) ) ); //TEMPORARY, JUST FOR DEBUGGING
    }
    else {
      node->setBranchProperty ( "Ev", BppString ( "S" ) );
    }
    //  std::cout << "annotateGeneTreeWithScoredDuplicationEvents No leaf 6"<<std::endl;

    return;
  }
}

//To sort in descending order
bool cmp ( int a, int b ) {
  return a > b;
}

//To sort in ascending order
bool anticmp ( int a, int b ) {
  return a < b;
}


/**************************************************************************/
VectorSiteContainer * getSequencesFromOptions ( map <string, string>  params, Alphabet* alphabet, bool& cont )
{

  VectorSiteContainer * sites;

  //Sequences and model of evolution
  std::string seqFile = ApplicationTools::getStringParameter ( "input.sequence.file",params,"none" );
  if ( !FileTools::fileExists ( seqFile ) ) {
    std::cerr << "Error: Sequence file "<< seqFile <<" not found." << std::endl;
    cont = false;
    /*          MPI::COMM_WORLD.Abort(1);
     *                  exit(-1);*/
  }
  else {
    VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer ( alphabet, params );
    //Removing white spaces in the names
    std::vector<std::string> names = allSites->getSequencesNames();
    for (int i = 0 ; i < names.size() ; ++i) {
      names[i] = TextTools:: removeSurroundingWhiteSpaces(names[i]);
    }
    allSites->setSequencesNames(names);

    unsigned int numSites = allSites->getNumberOfSites();
    ApplicationTools::displayResult ( "Number of sequences", TextTools::toString ( allSites->getNumberOfSequences() ) );
    ApplicationTools::displayResult ( "Number of sites", TextTools::toString ( numSites ) );

    std::vector <std::string> seqsToRemove;

    if ( numSites == 0 ) {
      std::cout<<"WARNING: Discarding a family whose alignment is 0 site long: "<< seqFile <<std::endl;
      cont = false;
      delete allSites;
    }
    else {
      unsigned int minPercentSequence = ApplicationTools::getIntParameter ( "sequence.removal.threshold",params,0 );
      unsigned int threshold = ( int ) ( ( double ) minPercentSequence * ( double ) numSites / 100 );

      if ( minPercentSequence > 0 ) {
        for ( int j = allSites->getNumberOfSequences()-1 ; j >= 0 ; j-- ) {
          if ( SequenceTools::getNumberOfCompleteSites ( allSites->getSequence ( j ) ) < threshold ) {
            ApplicationTools::displayResult ( "Removing a short sequence:", allSites->getSequence ( j ).getName() );
            // allSites->deleteSequence(i);
            seqsToRemove.push_back ( allSites->getSequence ( j ).getName() );
          }
        }
      }

      for ( unsigned int j =0 ; j<seqsToRemove.size(); j++ ) {
        std::vector <std::string> seqNames = allSites->getSequencesNames();
        if ( VectorTools::contains ( seqNames, seqsToRemove[j] ) ) {
          allSites->deleteSequence ( seqsToRemove[j] );
        }
        else
          std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene alignment."<<std::endl;
      }


      ApplicationTools::displayResult ( "# sequences post size-based removal:", TextTools::toString ( allSites->getNumberOfSequences() ) );
      seqsToRemove.clear();

      if ( allSites->getNumberOfSequences() <= 1 ) {
        std::cout << "Only one sequence left: discarding gene family "<<std::endl;
        cont = false;
      }
      else {
        sites = SequenceApplicationTools::getSitesToAnalyse ( *allSites, params );
        delete allSites;

        cont = true;
      }
    }
  }
  return sites;

}


/**************************************************************************/
SubstitutionModel*   getModelFromOptions ( map <string,string> params, Alphabet *alphabet, VectorSiteContainer *sites, bool& cont )
{
  /****************************************************************************
   *    //Then we need to get the substitution model.
   *****************************************************************************/
  SubstitutionModel* model = PhylogeneticsApplicationTools::getSubstitutionModel ( alphabet, 00, sites, params );
  cont = true;

  return model;
}


/**************************************************************************/
DiscreteDistribution* getRateDistributionFromOptions ( map <string,string> params, SubstitutionModel* model, bool& cont )
{
  DiscreteDistribution* rDist    = 0;
  if ( model->getNumberOfStates() > model->getAlphabet()->getSize() ) {
    // Markov-modulated Markov model!
    rDist = new ConstantDistribution ( 1. );
  }
  else {
    rDist = PhylogeneticsApplicationTools::getRateDistribution ( params );
  }
  cont = true;
  return rDist;
}


/**************************************************************************/
TreeTemplate<Node> * getTreeFromOptions ( map <string,string> params, Alphabet *alphabet, VectorSiteContainer * sites, SubstitutionModel* model, DiscreteDistribution* rDist, bool& cont )
{
  string file = ApplicationTools::getStringParameter ( "input.sequence.file",params,"none" );

  if (sites->getSequencesNames().size() < 3 || sites->getSequencesNames()[0].size() < 3) {
    std::cout << "getTreeFromOptions todobenoit hacking matrix < 3 case exception" << std::endl;
    return 0;
  }
  TreeTemplate<Node> *  rootedTree = 00;
  // Get the initial gene tree
  string initTree = ApplicationTools::getStringParameter ( "init.gene.tree", params, "user", "", false, false );
  ApplicationTools::displayResult ( "Input gene tree", initTree );
  TreeTemplate<Node> * unrootedGeneTree = 0;
  if ( initTree == "user" ) {
    std::string geneTree_File =ApplicationTools::getStringParameter ( "gene.tree.file",params,"none" );
    if ( geneTree_File=="none" ) {
      std::cout << "\n\nNo Gene tree was provided. The option init.gene.tree is set to user (by default), which means that the option gene.tree.file must be filled with the path of a valid tree file. \nIf you do not have a gene tree file, the program can start from a random tree, if you set init.gene.tree at random, or can build a gene tree with BioNJ or a PhyML-like algorithm with options bionj or phyml.\n\n" << std::endl;
      //MPI::COMM_WORLD.Abort ( 1 ); //SHOULD BE CORRECTED 13062017
      exit ( -1 );
    }
    Newick newick ( true );
    if ( !FileTools::fileExists ( geneTree_File ) ) {
      std::cerr << "Error: geneTree_File "<< geneTree_File <<" not found." << std::endl;
      std::cerr << "Building a bionj tree instead for gene " << geneTree_File << std::endl;
      unrootedGeneTree = buildBioNJTree ( params, sites, model, rDist, alphabet );
      if ( rootedTree ) {
        delete rootedTree;
        rootedTree =0;
      }
      rootedTree = unrootedGeneTree->clone();
      rootedTree->newOutGroup ( 0 );
      //exit(-1);
    }
    else {
      if ( rootedTree ) {
        delete rootedTree;
        rootedTree =0;
      }
      rootedTree = dynamic_cast < TreeTemplate < Node > * > ( newick.read ( geneTree_File ) );
    }
    if ( !rootedTree->isRooted() ) {
      if ( unrootedGeneTree ) {
        delete unrootedGeneTree;
        unrootedGeneTree = 0;
      }
      unrootedGeneTree = rootedTree->clone();
      //std::cout << "The gene tree is not rooted ; the root will be searched."<<std::endl;
      rootedTree->newOutGroup ( 0 );
    }
    else {
      if ( unrootedGeneTree ) {
        delete unrootedGeneTree;
        unrootedGeneTree = 0;
      }
      unrootedGeneTree = rootedTree->clone();
      unrootedGeneTree->unroot();
    }
    ApplicationTools::displayResult ( "Gene Tree file", geneTree_File );
    ApplicationTools::displayResult ( "Number of leaves", TextTools::toString ( rootedTree->getNumberOfLeaves() ) );
  }
  else if ( ( initTree == "bionj" ) || ( initTree == "phyml" ) ) { //build a BioNJ starting tree, and possibly refine it using PhyML algorithm
    unrootedGeneTree = buildBioNJTree ( params, sites, model, rDist, alphabet );

    if ( initTree == "phyml" ) { //refine the tree using PhyML algorithm (2003)
      refineGeneTreeUsingSequenceLikelihoodOnly ( params, unrootedGeneTree, sites, model, rDist, file, alphabet );
    }
    if ( rootedTree ) {
      delete rootedTree;
      rootedTree = 0;
    }
    rootedTree = unrootedGeneTree->clone();
    breadthFirstreNumber ( *rootedTree );
    //  std::cout << " Problem tree? : "<<TreeTemplateTools::treeToParenthesis(*rootedTree, true) << std::endl;

    rootedTree->newOutGroup ( rootedTree->getLeavesId() [0] );
  }
  else throw Exception ( "Unknown init gene tree method. init.gene.tree should be 'user', 'bionj', or 'phyml'." );
  //   }
  cont = true;
  return rootedTree;
}



/**************************************************************************/
void getCorrespondanceSequenceSpeciesFromOptions ( map< string, string > params, bool& cont, map <string, string> &seqSp, map<string, deque<string> > &spSeq )
{

  /****************************************************************************
   *    //Then we need to get the file containing links between sequences and species.
   *****************************************************************************/
  std::string taxaseqFile = ApplicationTools::getStringParameter ( "taxaseq.file",params,"none" );
  if ( taxaseqFile=="none" ) {
    std::cout << "\n\nNo taxaseqfile was provided. Cannot compute a reconciliation between a species tree and a gene tree using sequences if the relation between the sequences and the species is not explicit !\n" << std::endl;
    std::cout << "phyldog species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
    std::cerr << "\n\nNo taxaseqfile was provided. Cannot compute a reconciliation between a species tree and a gene tree using sequences if the relation between the sequences and the species is not explicit !\n" << std::endl;
    std::cerr << "phyldog species.tree.file=bigtree taxaseq.file=taxaseqlist gene.tree.file= genetree sequence.file=sequences.fa output.tree.file=outputtree\n"<<std::endl;
    cont = false;
    /* MPI::COMM_WORLD.Abort(1);
     *            exit(-1);*/
  }
  if ( !FileTools::fileExists ( taxaseqFile ) ) {
    std::cerr << "Error: taxaseqfile "<< taxaseqFile <<" not found." << std::endl;
    std::cout << "Error: taxaseqfile "<< taxaseqFile <<" not found." << std::endl;
    cont = false;
    /*  MPI::COMM_WORLD.Abort(1);
     *          exit(-1);*/
  }

  if ( cont ) {
    //Getting the relations between species and sequence names
    //In this file, the format is expected to be as follows :
    /*
     *         SpeciesA:sequence1;sequence2
     *         SpeciesB:sequence5
     *         SpeciesC:sequence3;sequence4;sequence6
     *         ...
     */
    //We use a std::map to record the links between species names and sequence names
    //For one species name, we can have several sequence names

    std::ifstream inSpSeq ( taxaseqFile.c_str() );
    std::string line;
    while ( getline ( inSpSeq,line ) ) {
      //We divide the line in 2 : first, the species name, second the sequence names
      StringTokenizer st1 ( line, ":", true );
      //Then we divide the sequence names
      if ( st1.numberOfRemainingTokens () >1 ) {
        StringTokenizer st2 ( st1.getToken ( 1 ), ";", true );
        if ( spSeq.find ( st1.getToken ( 0 ) ) == spSeq.end() ) {
          deque<string> temp = st2.getTokens();
          for ( unsigned int j = 0 ; j < temp.size() ; j++ ) {
            temp[j] = TextTools::removeSurroundingWhiteSpaces( temp[j] );
          }
          spSeq.insert ( make_pair ( TextTools::removeSurroundingWhiteSpaces( st1.getToken ( 0 ) ), temp ));
        }
        else {
          for ( unsigned int j = 0 ; j < ( st2.getTokens() ).size() ; j++ )
            spSeq.find ( TextTools::removeSurroundingWhiteSpaces( st1.getToken ( 0 ) ) )->second.push_back ( TextTools::removeSurroundingWhiteSpaces( st2.getTokens() [j] ) );
        }
      }
    }
    //Building seqSp
    for ( std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++ ) {
      for ( std::deque<std::string >::iterator it2 = ( it->second ).begin(); it2 != ( it->second ).end(); it2++ ) {
        seqSp.insert ( make_pair ( *it2, it->first ) );
      }
    }
  }
}


void removeUselessSequencesFromAlignment ( const TreeTemplate<Node>* spTree, bpp::VectorSiteContainer * sites, bool& cont, std::map<std::string, std::deque<std::string> > &spSeq, std::string file) {

  std::vector <std::string> seqsToRemove;

  //At the same time, we gather sequences we will have to remove from the
  //alignment and from the gene tree
  std::vector <std::string> spNamesToTake = spTree->getLeavesNames();
  for ( std::map<std::string, std::deque<std::string> >::iterator it = spSeq.begin(); it != spSeq.end(); it++ ) {
    if ( !VectorTools::contains ( spNamesToTake,it->first ) ) {
      for ( std::deque<std::string >::iterator it2 = ( it->second ).begin(); it2 != ( it->second ).end(); it2++ ) {
        std::cout<<"Removing sequence of species not considered: " << it->first <<std::endl;
        seqsToRemove.push_back ( *it2 );
      }
    }
  }
  //If we need to remove all sequences or all sequences except one,
  //better remove the gene family
  if ( seqsToRemove.size() >=sites->getNumberOfSequences()-1 ) {
    cont=false;
    std::cout <<"All or almost all sequences have been removed: avoiding family "<< file <<std::endl;
  }

  if ( cont ) {
    //We need to prune the alignment so that they contain
    //only sequences from the species under study.
    for ( unsigned int j =0 ; j<seqsToRemove.size(); j++ ) {
      std::vector <std::string> seqNames = sites->getSequencesNames();
      if ( VectorTools::contains ( seqNames, seqsToRemove[j] ) ) {
        sites->deleteSequence ( seqsToRemove[j] );
      }
      else
        std::cout<<"Sequence "<<seqsToRemove[j] <<"is not present in the gene alignment."<<std::endl;
    }
  }
  return;

}


void qualityControlGeneTree ( TreeTemplate<Node>* geneTree, bpp::VectorSiteContainer * sites, bool& cont, std::string file) {
  std::vector <std::string> seqsToRemove;
  //Going through the gene tree to see if leaves have branches that are too long.
  std::vector <Node*> leaves = geneTree->getLeaves();
  // std::cout << "leaves.size(): "<<leaves.size() <<std::endl;
  for ( unsigned int j = 0 ; j < leaves.size() ; j++ ) {
    if ( leaves[j] -> hasFather() && leaves[j]->hasDistanceToFather() && leaves[j]->getDistanceToFather() >= 2.0 ) {
      std::cout << "WARNING: Removing sequence "<< leaves[j]->getName() <<" from family "<<file<< " because its branch is unreasonably long (>=2.0)."<<std::endl;
      seqsToRemove.push_back ( leaves[j]->getName() );
      //removing the corresponding sequence, if present
      if ( sites->hasSequence ( leaves[j]->getName() ) )
        sites->deleteSequence ( leaves[j]->getName() );
    }

    // Removing sequences without leaves in the tree.
    if ( ! sites->hasSequence ( leaves[j]->getName() )) {
        seqsToRemove.push_back ( leaves[j]->getName() );
    }

  }

  if ( sites->getNumberOfSequences() >1 ) {
    //Pruning sequences from the gene tree
    for ( unsigned int j =0 ; j<seqsToRemove.size(); j++ ) {
      std::vector <std::string> leafNames = geneTree->getLeavesNames();
      if ( VectorTools::contains ( leafNames, seqsToRemove[j] ) ) {
        removeLeaf ( *geneTree, seqsToRemove[j] );
        /* if (unrootedGeneTree) {
         *                     delete unrootedGeneTree;
         *                     unrootedGeneTree = 0;
      }
      unrootedGeneTree = geneTree->clone();
      if (!geneTree->isRooted()) {
        std::cout <<"gene tree is not rooted!!! "<< taxaseqFile<<std::endl;
      }
      unrootedGeneTree->unroot();*/
      }
      else
        std::cout<<"Sequence "<<seqsToRemove[j] <<" is not present in the gene tree."<<std::endl;
    }
  }
  //If we have only one sequence in the end, we do not make a tree
  else {
    cont=false;
    std::cout <<"All or almost all sequences have been removed: avoiding family "<< file <<std::endl;
  }
  return;

}


bpp::Alphabet* getAlphabetFromOptions ( std::map <std::string, std::string>  params, bool& cont)
{
  Alphabet *alphabet = bpp::SequenceApplicationTools::getAlphabet ( params, "", false );
  cont = true;
  return alphabet;
}




/******************************************************************************/
// a is the species in which there was a duplication
// just below a, there was a loss
int recoverLossClosestToDuplication(TreeTemplate<Node> * spTree, int a, int aold) {
    Node* nodea = spTree->getNode ( a );
    std::vector < Node *> sonsA = nodea->getSons();
    std::vector < int > underlyingNodes = TreeTools::getNodesId(*spTree, sonsA[0]->getId());
    int lostBranch;
    if (VectorTools::contains<int>(underlyingNodes, aold)) {
      lostBranch = sonsA[1]->getId();
    }else {
      lostBranch = sonsA[0]->getId();
    }
    return lostBranch;
}

/******************************************************************************/

void setNDProperty(TreeTemplate<Node> * tree) {
  std::vector<Node* > nodes = tree->getNodes();
  //for (auto n = nodes.begin(); n!= nodes.end(); ++n) {
  for ( unsigned int i = 0 ; i < nodes.size() ; i++ ) {
      nodes[i]->setNodeProperty ( "ND", BppString ( TextTools::toString ( nodes[i]->getId() ) ) );
  }
  return;
}

/******************************************************************************/
// Recovers the events of duplication and loss for a given gene tree wrt a species tree.



      //vector on branches, map on sites, last vector on substitution types

vector < std::string > recoverDuplicationsAndLosses(
                    TreeTemplate<Node> * cTree,
                    TreeTemplate<Node> * spTree,
                    const string &familyName)
{
  setNDProperty(cTree);
  vector<int> ids = cTree->getNodesId();
  ids.pop_back(); //remove root id.

  vector < std::string > outputMatrix;
  string line = "";
  size_t numNodes = cTree->getNumberOfNodes();
  for (size_t i = 0; i < numNodes; ++i) {
    Node * node = cTree->getNode(ids[i]);
     if ( node->hasNodeProperty("S") ) {
        if (node->getNumberOfSons()>0 ) { // node has children, could be a duplication, and could contain losses
          int fatherSpID = TextTools::toInt ( ( dynamic_cast<const BppString*> ( node->getNodeProperty ( "S" ) ) )->toSTL() );
          std::vector <Node *> sons = node->getSons();
          int a = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[0]->getNodeProperty ( "S" ) ) )->toSTL() );
          int b = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[1]->getNodeProperty ( "S" ) ) )->toSTL() );
         // std::cout << "fatherSpID: " << fatherSpID<< "; a: "<< a << "; b: "<< b << std::endl;;
          int aold = a;
          int bold = b;
         // std::cout << "aold: "<< aold << " bold: "<< bold << std::endl;
          while ( a!=b ) {
              if ( a>b ) {
                int olda = a;
                Node* nodea = spTree->getNode ( a );
                a = nodea->getFather()->getId();
                if (a == fatherSpID) {  }
                else {
                  std::vector <Node *> sonsA = nodea->getFather()->getSons();
                  int lostBranch;
                  if (sonsA[0]->getId() == olda) {
                    lostBranch = sonsA[1]->getId();
                  }
                  else {
                    lostBranch = sonsA[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
                outputMatrix.push_back(line);
                      // lossA = lossA +1;
                }
              }
              else {
                int oldb = b;
                Node* nodeb = spTree->getNode ( b );
                b = nodeb->getFather()->getId();
                if (b == fatherSpID) {   }
                else {
                  std::vector <Node *> sonsb = nodeb->getFather()->getSons();
                  int lostBranch;
                  if (sonsb[0]->getId() == oldb) {
                    lostBranch = sonsb[1]->getId();
                  }
                  else {
                    lostBranch = sonsb[0]->getId();
                  }
                  line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
              outputMatrix.push_back(line);
                }
          //                lossB = lossB + 1;
              }
            }
        /*  sons[0]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossA ) ) );
          sons[1]->setBranchProperty ( "L", BppString ( TextTools::toString ( lossB ) ) );
          node->setNodeProperty ( "S", BppString ( TextTools::toString ( a ) ) );*/
          if ( ( a == aold ) || ( a == bold ) ) { //There was a duplication
            node->setBranchProperty ( "Ev", BppString ( "D" ) );
            int dupBranch = a;
            line = "event(" +  TextTools::toString(dupBranch) + ",\"" + familyName + "\"," + "duplication" + ")" ;
          outputMatrix.push_back(line);

            //We also need to check whether there have been losses
            if (aold > bold) { //loss in the lineage leading to a
                       /* std::cout << "BIS aold: "<< aold << " bold: "<< bold << std::endl;
                        std::cout << "BIS a: "<< a << " b: "<< b << std::endl;*/

              int lostBranch = recoverLossClosestToDuplication(spTree, a, aold);
              line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
              outputMatrix.push_back(line);

            }
            else if (bold>aold) { //loss in the lineage leading to b
                        //std::cout << "TER aold: "<< aold << " bold: "<< bold << std::endl;
                int lostBranch = recoverLossClosestToDuplication(spTree, b, bold);
                line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
              outputMatrix.push_back(line);
            }
          }
        }
      }
      else {
        std::cout << "No S Node property"<<std::endl;
      }
  }
  return outputMatrix;
}

/******************************************************************************/
// Outputs to file the numbers of duplications and losses per species.


void outputNumbersOfEventsPerFamilyPerSpecies( map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, bool temporary ) {

  WHEREAMI( __FILE__ , __LINE__ );
string suffix = ApplicationTools::getStringParameter ( "output.file.suffix", params, "", "", false, false );
string evFile = ApplicationTools::getStringParameter ( "output.events.file", params, "events.txt", "", false, false );
evFile = evFile + suffix;
if ( temporary ) {
  //   string temp = "temp";
  evFile = evFile + "_temp";
}

outputNumbersOfEventsToFile( params, geneTree,  speciesTree, seqSp, familyName, evFile);


}



void outputNumbersOfEventsToFile( map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, std::string& evFile ) {

  WHEREAMI( __FILE__ , __LINE__ );
  std::ofstream out;

breadthFirstreNumber ( *speciesTree );
std::map <std::string, int> spId = computeSpeciesNamesToIdsMap ( *speciesTree );

annotateGeneTreeWithDuplicationEvents ( *speciesTree,
                                        *geneTree,
                                        geneTree->getRootNode(),
                                        seqSp,
                                        spId );

vector < std::string > lines = recoverDuplicationsAndLosses(geneTree, speciesTree, familyName);
out.open ( evFile.c_str(), std::ios::out );
for (unsigned int i = 0; i < lines.size(); ++i) {
  out << lines[i] <<std::endl;
}
out.close();
return;


}



/******************************************************************************/
// Recovers the events of duplication for a given gene tree wrt a species tree,
// so that we can tell which genes are orthologous and which genes are paralogous.
// The output we want is as follows:
/*
DUPLICATIONS : 4        LOSSES : 20
    ORTHOLOGY RELATIONSHIP: ENSOGAP00000018864_F05772@Otolemur_garnettii, ENSSSCP00000021040_F05772@Sus_scrofa, ENSP00000325562_F05772@Homo_sapiens, ENSSTOP00000023581_F05772@Ictidomys_tridecemlineatus, ENSMLUP00000012455_F05772@Myotis_lucifugus, ENSCJAP00000012390_F05772@Callithrix_jacchus, APAMM00000012449_F05772@Mus_musculus, ENSMLUP00000019915_F05772@Myotis_lucifugus, ENSMPUP00000016053_F05772@Mustela_putorius_furo, ENSFCAP00000002823_F05772@Felis_catus, ENSOARP00000015537_F05772@Ovis_aries, ENSCPOP00000007201_F05772@Cavia_porcellus    <====>    ENSDNOP00000033050_F05772@Dasypus_novemcinctus, ENSDNOP00000021811_F05772@Dasypus_novemcinctus, ENSLAFP00000013597_F05772@Loxodonta_africana
    PARALOGY RELATIONSHIP: ENSSTOP00000023581_F05772@Ictidomys_tridecemlineatus, ENSMPUP00000016053_F05772@Mustela_putorius_furo, ENSOGAP00000018864_F05772@Otolemur_garnettii    <====>    ENSSSCP00000021040_F05772@Sus_scrofa, ENSP00000325562_F05772@Homo_sapiens, ENSMLUP00000012455_F05772@Myotis_lucifugus, ENSCJAP00000012390_F05772@Callithrix_jacchus, APAMM00000012449_F05772@Mus_musculus, ENSMLUP00000019915_F05772@Myotis_lucifugus, ENSFCAP00000002823_F05772@Felis_catus, ENSOARP00000015537_F05772@Ovis_aries, ENSCPOP00000007201_F05772@Cavia_porcellus
...
 */
vector < std::string > recoverOrthologsAndParalogs(
                    TreeTemplate<Node> * cTree,
                    TreeTemplate<Node> * spTree,
                    const string &familyName)
{
  setNDProperty(cTree);
  vector<int> ids = cTree->getNodesId();
  ids.pop_back(); //remove root id.

  vector < std::string > outputMatrix;
  string line = "";
  size_t numDups = 0;
  size_t numLosses = 0;
  size_t numNodes = cTree->getNumberOfNodes();
  for (size_t i = 0; i < numNodes; ++i) {
    Node * node = cTree->getNode(ids[i]);
     if ( node->hasNodeProperty("S") ) {
        if (node->getNumberOfSons()>0 ) { // node has children, could be a duplication, and could contain losses
          int fatherSpID = TextTools::toInt ( ( dynamic_cast<const BppString*> ( node->getNodeProperty ( "S" ) ) )->toSTL() );
          std::vector <Node *> sons = node->getSons();
          int a = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[0]->getNodeProperty ( "S" ) ) )->toSTL() );
          int b = TextTools::toInt ( ( dynamic_cast<const BppString*> ( sons[1]->getNodeProperty ( "S" ) ) )->toSTL() );
         // std::cout << "fatherSpID: " << fatherSpID<< "; a: "<< a << "; b: "<< b << std::endl;;
          int aold = a;
          int bold = b;
         // std::cout << "aold: "<< aold << " bold: "<< bold << std::endl;
          while ( a!=b ) {
              if ( a>b ) {
                int olda = a;
                Node* nodea = spTree->getNode ( a );
                a = nodea->getFather()->getId();
                if (a == fatherSpID) {  }
                else {
                  std::vector <Node *> sonsA = nodea->getFather()->getSons();
                  int lostBranch;
                  if (sonsA[0]->getId() == olda) {
                    lostBranch = sonsA[1]->getId();
                  }
                  else {
                    lostBranch = sonsA[0]->getId();
                  }
		  //line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
		numLosses += 1;
                }
              }
              else {
                int oldb = b;
                Node* nodeb = spTree->getNode ( b );
                b = nodeb->getFather()->getId();
                if (b == fatherSpID) {   }
                else {
                  std::vector <Node *> sonsb = nodeb->getFather()->getSons();
                  int lostBranch;
                  if (sonsb[0]->getId() == oldb) {
                    lostBranch = sonsb[1]->getId();
                  }
                  else {
                    lostBranch = sonsb[0]->getId();
                  }
                  //line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
		  numLosses += 1;
                }
              }
            }
          if ( ( a == aold ) || ( a == bold ) ) { //There was a duplication
            node->setBranchProperty ( "Ev", BppString ( "D" ) );
            int dupBranch = a;
	    numDups += 1;
	    line = familyName + "\tPARALOGY RELATIONSHIP: ";
	    // Assuming binary tree
	    std::vector<std::string> leafNames = TreeTemplateTools::getLeavesNames(*(node->getSon(0)));
	    for (size_t j = 0; j < leafNames.size()-1; ++j) {
	      line = line + leafNames[j] +", ";
	    }
	    line = line + leafNames[leafNames.size()-1] + " <===> ";
	    leafNames = TreeTemplateTools::getLeavesNames(*(node->getSon(1)));
	    for (size_t j = 0; j < leafNames.size()-1; ++j) {
	      line = line + leafNames[j] +", ";
	    }
	    line = line + leafNames[leafNames.size()-1];
	    outputMatrix.push_back(line);
            //line = "event(" +  TextTools::toString(dupBranch) + ",\"" + familyName + "\"," + "duplication" + ")" ;

            //We also need to check whether there have been losses
            if (aold > bold) { //loss in the lineage leading to a
                       /* std::cout << "BIS aold: "<< aold << " bold: "<< bold << std::endl;
                        std::cout << "BIS a: "<< a << " b: "<< b << std::endl;*/

              int lostBranch = recoverLossClosestToDuplication(spTree, a, aold);
              //line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
	      numLosses += 1;
            }
            else if (bold>aold) { //loss in the lineage leading to b
                        //std::cout << "TER aold: "<< aold << " bold: "<< bold << std::endl;
                int lostBranch = recoverLossClosestToDuplication(spTree, b, bold);
                //line = "event(" +  TextTools::toString(lostBranch) + ",\"" + familyName + "\"," + "loss" + ")" ;
		numLosses += 1;
            }
          }
	  else { // There was no duplication, hence we have a speciation node
	    	    line = familyName + "\tORTHOLOGY RELATIONSHIP: ";
	    // Assuming binary tree
		    std::vector<std::string> leafNames = TreeTemplateTools::getLeavesNames(*(node->getSon(0)));
		    for (size_t j = 0; j < leafNames.size()-1; ++j) {
	      line = line + leafNames[j] +", ";
	    }
		    line = line + leafNames[leafNames.size()-1] + " <===> ";
	    leafNames = TreeTemplateTools::getLeavesNames(*(node->getSon(1)));
	    for (size_t j = 0; j < leafNames.size()-1; ++j) {
	      line = line + leafNames[j] +", ";
	    }
	    line = line + leafNames[leafNames.size()-1];
	    outputMatrix.push_back(line);
	  }
        }
      }
      else {
        std::cout << "No S Node property"<<std::endl;
      }
  }
  line = familyName + "\tDUPLICATIONS : "+TextTools::toString(numDups)+"\tLOSSES : "+TextTools::toString(numLosses);
  outputMatrix.push_back(line);
  return outputMatrix;
}



/******************************************************************************/
// Outputs to file orthologous and paralogous genes
void outputOrthologousAndParalogousGenes(map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, bool temporary ) {

  WHEREAMI( __FILE__ , __LINE__ );
  string suffix = ApplicationTools::getStringParameter ( "output.file.suffix", params, "", "", false, false );
  string orFile = ApplicationTools::getStringParameter ( "output.orthologs.file", params, "orthologs.txt", "", false, false );
  orFile = orFile + suffix;
  if ( temporary ) {
    //   string temp = "temp";
    orFile = orFile + "_temp";
  }

  outputOrthologousAndParalogousGenesToFile(params, geneTree, speciesTree, seqSp, familyName, orFile);

}


void outputOrthologousAndParalogousGenesToFile(map<string, string > params, TreeTemplate<Node> *geneTree,  TreeTemplate<Node> *speciesTree, const std::map <std::string, std::string> seqSp, std::string& familyName, std::string& orFile ) {

  WHEREAMI( __FILE__ , __LINE__ );
  std::ofstream out;

  breadthFirstreNumber ( *speciesTree );
  std::map <std::string, int> spId = computeSpeciesNamesToIdsMap ( *speciesTree );

  annotateGeneTreeWithDuplicationEvents ( *speciesTree,
    *geneTree,
    geneTree->getRootNode(),
    seqSp,
    spId );

    vector < std::string > lines = recoverOrthologsAndParalogs(geneTree, speciesTree, familyName);

    out.open ( orFile.c_str(), std::ios::out );
    out << lines[lines.size()-1] <<std::endl;
    for (unsigned int i = 0; i < lines.size()-1; ++i) {
      out << lines[i] <<std::endl;
    }
    out.close();
    return;

}
