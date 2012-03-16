//
// File: DLGeneTreeLikelihood.cpp
// Created by: Bastien Boussau
// Created on: Tue Oct 04 18:14:51 2011
//

/*
 Copyright or � or Copr. CNRS, (November 16, 2004)
 
 This software is a computer program whose purpose is to provide classes
 for phylogenetic data analysis.
 
 This software is governed by the CeCILL  license under French law and
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


#include "DLGeneTreeLikelihood.h"
#include "GenericTreeExplorationAlgorithms.h"

// From Utils:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>

// From NumCalc:
#include <Bpp/Numeric/AutoParameter.h>

// From the STL:
#include <iostream>
//using namespace std;
using namespace bpp;

/******************************************************************************/
DLGeneTreeLikelihood::DLGeneTreeLikelihood(
                                                           const Tree & tree,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
                                                           TreeTemplate<Node> & geneTreeWithSpNames,
                                                           const std::map <std::string, std::string> seqSp,
                                                           std::map <std::string,int> spId,
                                                           //std::vector <int> & lossNumbers, 
                                                           std::vector <double> & lossProbabilities, 
                                                           //std::vector <int> & duplicationNumbers, 
                                                           std::vector <double> & duplicationProbabilities,
                                                           //std::vector <int> & branchNumbers,
                                                           std::vector <int> & num0Lineages,
                                                           std::vector <int> & num1Lineages,
                                                           std::vector <int> & num2Lineages, 
                                                           int speciesIdLimitForRootPosition,
                                                           int heuristicsLevel,
                                                           int & MLindex, 
                                                           bool checkRooted,
                                                           bool verbose, 
                                                           bool rootOptimization, 
                                                           bool considerSequenceLikelihood, 
                                                           bool DLStartingGeneTree, 
                                                           unsigned int sprLimit)
throw (Exception):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp(seqSp), _spId(spId)
{
    
    nniLk_ = new NNIHomogeneousTreeLikelihood(tree, model, rDist, checkRooted, verbose); 
    _spTree = spTree.clone();
    _rootedTree = rootedTree.clone();
    _geneTreeWithSpNames = geneTreeWithSpNames.clone();
    //_lossNumbers = lossNumbers;
    _lossProbabilities = lossProbabilities;
    //_duplicationNumbers = duplicationNumbers;
    _duplicationProbabilities = duplicationProbabilities;
    //_branchNumbers = branchNumbers;
    _num0Lineages=num0Lineages;
    _num1Lineages=num1Lineages;
    _num2Lineages=num2Lineages;
    _scenarioLikelihood = UNLIKELY;
   // _sequenceLikelihood = UNLIKELY;
    _MLindex = MLindex;
    _rootOptimization = rootOptimization; 
    /*_tentativeDuplicationNumbers = duplicationNumbers;
     _tentativeLossNumbers = lossNumbers; 
     _tentativeBranchNumbers = branchNumbers;*/
    _tentativeNum0Lineages =num0Lineages;
    _tentativeNum1Lineages =num1Lineages; 
    _tentativeNum2Lineages =num2Lineages;
    _tentativeMLindex = MLindex;
    _totalIterations = 0;
    _counter = 0;
    _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
    _heuristicsLevel = heuristicsLevel;
    _optimizeSequenceLikelihood = true;
    _optimizeReconciliationLikelihood = true;
    _considerSequenceLikelihood = considerSequenceLikelihood;
    _DLStartingGeneTree = DLStartingGeneTree;
    sprLimit_ = sprLimit;
    // _listOfPreviousRoots = new std::vector <int> ();
}

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(
                                                           const Tree & tree,
                                                           const SiteContainer & data,
                                                           SubstitutionModel * model,
                                                           DiscreteDistribution * rDist,
                                                           TreeTemplate<Node> & spTree,
                                                           TreeTemplate<Node> & rootedTree,
                                                           TreeTemplate<Node> & geneTreeWithSpNames,
                                                           const std::map <std::string, std::string> seqSp,
                                                           std::map <std::string,int> spId,
                                                           //std::vector <int> & lossNumbers, 
                                                           std::vector <double> & lossProbabilities, 
                                                           //std::vector <int> & duplicationNumbers, 
                                                           std::vector <double> & duplicationProbabilities, 
                                                           //std::vector <int> & branchNumbers, 
                                                           std::vector <int> & num0Lineages,
                                                           std::vector <int> & num1Lineages,
                                                           std::vector <int> & num2Lineages, 
                                                           int speciesIdLimitForRootPosition, 
                                                           int heuristicsLevel,
                                                           int & MLindex, 
                                                           bool checkRooted,
                                                           bool verbose,
                                                           bool rootOptimization, 
                                                           bool considerSequenceLikelihood, 
                                                           bool DLStartingGeneTree,
                                                           unsigned int sprLimit)
throw (Exception):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp (seqSp), _spId(spId)
{
    nniLk_ = new NNIHomogeneousTreeLikelihood(tree, data, model, rDist, checkRooted, verbose);
    _spTree = spTree.clone();
    _rootedTree = rootedTree.clone();
    _geneTreeWithSpNames = geneTreeWithSpNames.clone();
    _lossProbabilities = lossProbabilities;
    _duplicationProbabilities = duplicationProbabilities; 
    _num0Lineages=num0Lineages;
    _num1Lineages=num1Lineages;
    _num2Lineages=num2Lineages;
    _scenarioLikelihood = UNLIKELY;
  //  _sequenceLikelihood = UNLIKELY;
    _MLindex = MLindex;
    _rootOptimization = rootOptimization; 
    _tentativeNum0Lineages =num0Lineages;
    _tentativeNum1Lineages =num1Lineages; 
    _tentativeNum2Lineages =num2Lineages;
    _tentativeMLindex = MLindex;
    _totalIterations = 0; 
    _counter = 0;
    _speciesIdLimitForRootPosition = speciesIdLimitForRootPosition;
    _heuristicsLevel = heuristicsLevel;
    _optimizeSequenceLikelihood = true;
    _optimizeReconciliationLikelihood = true;
    _considerSequenceLikelihood = considerSequenceLikelihood;
    _DLStartingGeneTree = DLStartingGeneTree;
    sprLimit_ = sprLimit;
    
}

/******************************************************************************/

DLGeneTreeLikelihood::DLGeneTreeLikelihood(const DLGeneTreeLikelihood & lik):
nniLk_(0), _spTree(0), _rootedTree(0), _geneTreeWithSpNames(0), _seqSp (lik._seqSp), _spId(lik._spId)
{
    nniLk_ = lik.nniLk_->clone(); 
    _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone()) ;
    _rootedTree = dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone()) ;
    _geneTreeWithSpNames = dynamic_cast<TreeTemplate<Node> *> (lik._geneTreeWithSpNames->clone()) ;
    _lossProbabilities = lik._lossProbabilities;
    _duplicationProbabilities = lik._duplicationProbabilities; 
    _num0Lineages=lik._num0Lineages;
    _num1Lineages=lik._num1Lineages;
    _num2Lineages=lik._num2Lineages;
    _scenarioLikelihood = lik._scenarioLikelihood;
   // _sequenceLikelihood = lik._sequenceLikelihood;
    _MLindex = lik._MLindex;
    _rootOptimization = lik._rootOptimization; 
    _tentativeNum0Lineages =lik._tentativeNum0Lineages;
    _tentativeNum1Lineages =lik._tentativeNum1Lineages;
    _tentativeNum2Lineages =lik._tentativeNum2Lineages;
    _tentativeMLindex = lik._MLindex;
    _totalIterations = lik._totalIterations;
    _counter = lik._counter;
    _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
    _heuristicsLevel = lik._heuristicsLevel;
    _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
    _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
    _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
    _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;
    _considerSequenceLikelihood = lik._considerSequenceLikelihood;
    _DLStartingGeneTree = lik._DLStartingGeneTree;
    sprLimit_ = lik.sprLimit_;
}

/******************************************************************************/

DLGeneTreeLikelihood & DLGeneTreeLikelihood::operator=(const DLGeneTreeLikelihood & lik)
{
    if (nniLk_) delete nniLk_;
    nniLk_ = lik.nniLk_->clone(); 
    //NNIHomogeneousTreeLikelihood::operator=(lik);
    if (_spTree) delete _spTree;
    _spTree = dynamic_cast<TreeTemplate<Node> *> (lik._spTree->clone());
    if (_rootedTree) delete _rootedTree;
    _rootedTree= dynamic_cast<TreeTemplate<Node> *> (lik._rootedTree->clone());
    if (_geneTreeWithSpNames) delete _geneTreeWithSpNames;
    _geneTreeWithSpNames = dynamic_cast<TreeTemplate<Node> *> (lik._geneTreeWithSpNames->clone()) ;
    _spId = lik._spId;
    _lossProbabilities = lik._lossProbabilities;
    _duplicationProbabilities = lik._duplicationProbabilities;
    _num0Lineages=lik._num0Lineages;
    _num1Lineages=lik._num1Lineages;
    _num2Lineages=lik._num2Lineages;
    _scenarioLikelihood = lik._scenarioLikelihood;
   // _sequenceLikelihood = lik._sequenceLikelihood;
    _MLindex = lik._MLindex;
    _rootOptimization = lik._rootOptimization;
    _tentativeNum0Lineages =lik._tentativeNum0Lineages;
    _tentativeNum1Lineages =lik._tentativeNum1Lineages;
    _tentativeNum2Lineages =lik._tentativeNum2Lineages;
    _tentativeMLindex = lik._MLindex;
    _totalIterations = lik._totalIterations;
    _counter = lik._counter;
    _speciesIdLimitForRootPosition = lik._speciesIdLimitForRootPosition;
    _heuristicsLevel = lik._heuristicsLevel;
    _nodesToTryInNNISearch = lik._nodesToTryInNNISearch;
    _tentativeNodesToTryInNNISearch = lik._tentativeNodesToTryInNNISearch;
    _optimizeSequenceLikelihood = lik._optimizeSequenceLikelihood;
    _optimizeReconciliationLikelihood = lik._optimizeReconciliationLikelihood ;
    _considerSequenceLikelihood = lik._considerSequenceLikelihood;
    _DLStartingGeneTree = lik._DLStartingGeneTree;
    sprLimit_ = lik.sprLimit_;
    return *this;
}




/******************************************************************************/

DLGeneTreeLikelihood::~DLGeneTreeLikelihood()
{
    if (nniLk_) delete nniLk_;
    if (_spTree) delete _spTree;
    if (_rootedTree) delete _rootedTree;
    if (_geneTreeWithSpNames) delete _geneTreeWithSpNames;
}



/******************************************************************************/

void DLGeneTreeLikelihood::initParameters()
{
    // std::cout << "in initParameters"<<std::endl;
    if (_considerSequenceLikelihood) {
        nniLk_->initParameters();
    }
    if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
        exit(-1);
        //    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
    }
    else {
        _scenarioLikelihood = findMLReconciliationDR (_spTree, _rootedTree, 
                                                      _seqSp, _spId, _lossProbabilities, 
                                                      _duplicationProbabilities, _MLindex, 
                                                      _num0Lineages, _num1Lineages,
                                                      _num2Lineages, _nodesToTryInNNISearch); 
    }
    _MLindex = -1;
    // std::cout << "in ReconciliationTreeLikelihood::initParameters : _num0Lineages, _num1Lineages, _num2Lineages : "<< TextTools::toString(VectorTools::sum(_num0Lineages))<<" "<<TextTools::toString(VectorTools::sum(_num1Lineages))<<" "<< TextTools::toString(VectorTools::sum(_num2Lineages))<<std::endl;
    
//    std::cout <<"INITIAL _scenarioLikelihood "<<_scenarioLikelihood <<" nniLk_->getLogLikelihood(): "<< nniLk_->getLogLikelihood()<<std::endl;
}



/******************************************************************************/

void DLGeneTreeLikelihood::resetMLindex() {
    _MLindex = -1;
}




/******************************************************************************/
/* We need to introduce in the likelihood computation the scenario likelihood */
double DLGeneTreeLikelihood::getLogLikelihood() const
{
    double ll = 0;
    // _sequenceLikelihood=UNLIKELY;
    //  if (_sequenceLikelihood ==UNLIKELY) {
    
    //TEST
    /*
     if ((_sequenceLikelihood == UNLIKELY)||(_optimizeSequenceLikelihood==true)) {
     Vdouble * lik = & likelihoodData_->getRootRateSiteLikelihoodArray();
     const std::vector<unsigned int> * w = & likelihoodData_->getWeights();
     for(unsigned int i = 0; i < nbDistinctSites_; i++)
     {
     ll += (* w)[i] * log((* lik)[i]);
     //std::cout << i << "\t" << (* w)[i] << "\t" << log((* lik)[i]) << std::endl;
     }
     _sequenceLikelihood = ll;
     // std::cout <<"getLogLikelihood _sequenceLikelihood "<< _sequenceLikelihood<<std::endl;
     }
     */
    //END TEST
    
    // ll = getValue();
    
    /* else {
     std::cout <<"NOT RECOMPUTING IN GETLOGLIKELIHOOD : "<< _sequenceLikelihood<<std::endl; 
     }*/
    //Now we add the scenario likelihood
    //  std::cout << "COMPUTING LOGLIKELIHOOD :" << ll  << " SCENARIO LOGLIKELIHOOD :" << _scenarioLikelihood <<" TOTAL : "<< ll + _scenarioLikelihood <<std::endl;
    
    //TEST 16 02 2010
    // DRHomogeneousTreeLikelihood::getLogLikelihood();
    //  setMinuslogLikelihood_ (_sequenceLikelihood);
    if (_considerSequenceLikelihood) {
        ll = nniLk_->getLogLikelihood() + _scenarioLikelihood;
    }
    else {
        ll = _scenarioLikelihood;
    }
    // ll = _sequenceLikelihood ;
    return ll;
}
/******************************************************************************/
//returns -loglikelihood
//As a reminder: _sequenceLikelihood < 0, _scenarioLikelihood < 0, getValue >0
double DLGeneTreeLikelihood::getValue() const  
throw (Exception)
{
    {
        if(!nniLk_->isInitialized()) throw Exception("reconciliationTreeLikelihood::getValue(). Instance is not initialized.");
        // return (-getLogLikelihood());
        //TEST 16 02 2010
        // std::cout<<"\t\t\t_sequenceLikelihood: "<<_sequenceLikelihood<< " _scenarioLikelihood: "<<_scenarioLikelihood<<std::endl;
        if (_considerSequenceLikelihood) {
            return (- nniLk_->getLogLikelihood() - _scenarioLikelihood);
        }
        else {
            return (- _scenarioLikelihood);
        }
        //return (minusLogLik_ - _scenarioLikelihood);
        //return (-minusLogLik_);
    }
}



/******************************************************************************/


void DLGeneTreeLikelihood::fireParameterChanged(const ParameterList & params)
{
    if (_considerSequenceLikelihood) {
        nniLk_->applyParameters();
        nniLk_->matchParametersValues(params);
    }
    //If we need to update the reconciliation likelihood
    if (_optimizeReconciliationLikelihood) {
        computeReconciliationLikelihood();
    }
    
    
}


/******************************************************************************/


void DLGeneTreeLikelihood::computeSequenceLikelihood()
{
    if ( _considerSequenceLikelihood && (_optimizeSequenceLikelihood==true) ) {
        nniLk_->computeTreeLikelihood ();
        /*
        nniLk_->computeSubtreeLikelihoodPostfix(dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getRootNode());
        nniLk_->computeSubtreeLikelihoodPrefix(nniLk_->getTree().getRootNode());
        nniLk_->computeRootLikelihood();*/
    }
}




/******************************************************************************/

void DLGeneTreeLikelihood::computeReconciliationLikelihood()
{
    resetLossesAndDuplications(*_spTree, /*_lossNumbers, */_lossProbabilities, /*_duplicationNumbers, */_duplicationProbabilities);
    if (_heuristicsLevel>0) {
        std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
        MPI::COMM_WORLD.Abort(1);
        exit(-1);
        //    _scenarioLikelihood = findMLReconciliation (&_spTree, &_rootedTree, _seqSp, _lossNumbers, _lossProbabilities, _duplicationNumbers, _duplicationProbabilities, _MLindex, _branchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch); 
    }
    else {
        //    _scenarioLikelihood = findMLReconciliationDR (&_spTree, &_rootedTree, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _MLindex, _num0Lineages, _num1Lineages, _num2Lineages, _nodesToTryInNNISearch);
        _scenarioLikelihood = findMLReconciliationDR (_spTree, _rootedTree, 
                                                      _seqSp, _spId, 
                                                      _lossProbabilities, _duplicationProbabilities, 
                                                      _tentativeMLindex, 
                                                      _tentativeNum0Lineages, _tentativeNum1Lineages, 
                                                      _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
        _MLindex = _tentativeMLindex;
        _num0Lineages = _tentativeNum0Lineages;
        _num1Lineages = _tentativeNum1Lineages;
        _num2Lineages = _tentativeNum2Lineages;
        _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
    }
}



/******************************************************************************/


void DLGeneTreeLikelihood::computeTreeLikelihood()
{
    if (_considerSequenceLikelihood )
    {
        computeSequenceLikelihood();
    }
    computeReconciliationLikelihood();  
}

/******************************************************************************/
/*
 * This function tries a given NNI. It takes the rooted tree, makes an NNI on it, and computes the likelihood of the best scenario for this new topology. If this likelihood is better than the current scenario likelihood, the sequence likelihood is computed on the unrooted tree.
 
 */
double DLGeneTreeLikelihood::testNNI(int nodeId) const throw (NodeException)
{
    //int nodeId = son->getId();
    //  std::cout<<"IN TESTNNI, nodeId "<< nodeId << "_nodesToTryInNNISearch.size() "<< _nodesToTryInNNISearch.size()<<std::endl;
    // std::cout << "before "<<TreeTools::treeToParenthesis (*tree_, true)<<std::endl;
    //If the NNI is around a branch where a duplication was found, 
    //or if we just try all branches because the starting gene trees are parsimonious in
    //numbers of DL.
    if ((_nodesToTryInNNISearch.count(nodeId)==1) /*|| _DLStartingGeneTree*/) {
        TreeTemplate<Node> * treeForNNI = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
        
        _tentativeMLindex = _MLindex;
        /* _tentativeLossNumbers = _lossNumbers;
         _tentativeDuplicationNumbers = _duplicationNumbers;
         _tentativeBranchNumbers = _branchNumbers;*/
        _tentativeNum0Lineages = _num0Lineages;
        _tentativeNum1Lineages = _num1Lineages;
        _tentativeNum2Lineages =_num2Lineages;
        _tentativeNodesToTryInNNISearch.clear();
        
        //We first estimate the likelihood of the scenario: if not better than the current scenario, no need to estimate the branch length !
        //We use the same procedure as in doNNI !
        const Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getNode(nodeId);
        
        
        if(!son->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'son' must not be the root node.", nodeId);
        const Node * parent = son->getFather();
        
        if(!parent->hasFather()) throw NodeException("DLGeneTreeLikelihood::testNNI(). Node 'parent' must not be the root node.", parent->getId());
        const Node * grandFather = parent->getFather();
        
        //From here: Bifurcation assumed.
        //In case of multifurcation, an arbitrary uncle is chosen.
        //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
        unsigned int parentPosition = grandFather->getSonPosition(parent);
        
        Node * sonForNNI    = treeForNNI->getNode(nodeId);
        Node * parentForNNI = sonForNNI->getFather();
        Node * grandFatherForNNI = parentForNNI->getFather();
        Node * uncleForNNI = grandFatherForNNI->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
        parentForNNI->removeSon(sonForNNI);
        grandFatherForNNI->removeSon(uncleForNNI);
        parentForNNI->addSon(uncleForNNI);
        grandFatherForNNI->addSon(sonForNNI);
        
        //Now we root the tree sent to findMLReconciliation as in _rootedTree
        int id = treeForNNI->getRootNode()->getId();
        if(TreeTemplateTools::hasNodeWithId(*(_rootedTree->getRootNode()->getSon(0)),id)) {
            treeForNNI->newOutGroup(_rootedTree->getRootNode()->getSon(1)->getId());
        }
        else {
            treeForNNI->newOutGroup(_rootedTree->getRootNode()->getSon(0)->getId());
        }
        double ScenarioMLValue = 0;
        _totalIterations = _totalIterations+1;
        
        /* //If we want to optimize the root or if we are at the first try
         if (_rootOptimization){
         if (_heuristicsLevel>0) {
         std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
         exit(-1);
         // ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
         }
         else {
         ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
         }
         }
         else {
         if (_heuristicsLevel>0) {
         std::cout <<"Sorry, these heuristics are no longer available. Try option 0."<<std::endl;
         exit(-1);
         // ScenarioMLValue =  findMLReconciliation (&_spTree, treeForNNI, _seqSp, _tentativeLossNumbers, _lossProbabilities, _tentativeDuplicationNumbers, _duplicationProbabilities, _tentativeMLindex, _tentativeBranchNumbers, _speciesIdLimitForRootPosition, _heuristicsLevel, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
         }
         else {
         ScenarioMLValue =  findMLReconciliationDR (&_spTree, treeForNNI, _seqSp, _spId, _lossProbabilities, _duplicationProbabilities, _tentativeMLindex, _tentativeNum0Lineages, _tentativeNum1Lineages, _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch); 
         }
         
         
         
         
         }*/
        ScenarioMLValue =  findMLReconciliationDR (_spTree, treeForNNI/*&_rootedTree*/, 
                                                   _seqSp, _spId, 
                                                   _lossProbabilities, _duplicationProbabilities, 
                                                   _tentativeMLindex, 
                                                   _tentativeNum0Lineages, _tentativeNum1Lineages, 
                                                   _tentativeNum2Lineages, _tentativeNodesToTryInNNISearch, false); //false so that _tentativeNum*Lineages are not updated
        
        
        if (treeForNNI) delete treeForNNI;
       // std::cout<<"???WORTH computing the sequence likelihood "<< ScenarioMLValue<< " "<< _scenarioLikelihood<<std::endl;
        if (_considerSequenceLikelihood ) 
        {
            if  (ScenarioMLValue >  _scenarioLikelihood) 
            { //If it is worth computing the sequence likelihood
                //Retrieving arrays of interest:
               // std::cout << "before nniLk_->testNNI(nodeId);"<<std::endl;
                double newLkMinusOldLk = nniLk_->testNNI(nodeId);
               // std::cout << "after nniLk_->testNNI(nodeId);"<<std::endl;

                
                /*
                const DRASDRTreeLikelihoodNodeData * parentData = & nniLk_->getLikelihoodData()->getNodeData(parent->getId());
                const VVVdouble                    * sonArray   = & parentData->getLikelihoodArrayForNeighbor(son->getId());
                std::vector<const Node *> parentNeighbors = TreeTemplateTools::getRemainingNeighbors(parent, grandFather, son);
                unsigned int nbParentNeighbors = parentNeighbors.size();
                std::vector<const VVVdouble *> parentArrays(nbParentNeighbors);
                std::vector<const VVVdouble *> parentTProbs(nbParentNeighbors);
                for(unsigned int k = 0; k < nbParentNeighbors; k++)
                {
                    const Node * n = parentNeighbors[k]; // This neighbor
                    parentArrays[k] = & parentData->getLikelihoodArrayForNeighbor(n->getId()); 
                    parentTProbs[k] = & pxy_[n->getId()];
                }
                const DRASDRTreeLikelihoodNodeData * grandFatherData = & nniLk_->getLikelihoodData()->getNodeData(grandFather->getId());
                const VVVdouble                    * uncleArray      = & grandFatherData->getLikelihoodArrayForNeighbor(uncle->getId()); 
                std::vector<const Node *> grandFatherNeighbors = TreeTemplateTools::getRemainingNeighbors(grandFather, parent, uncle);
                unsigned int nbGrandFatherNeighbors = grandFatherNeighbors.size();
                
                std::vector<const VVVdouble *> grandFatherArrays;
                std::vector<const VVVdouble *> grandFatherTProbs;
                for(unsigned int k = 0; k < nbGrandFatherNeighbors; k++)
                {
                    const Node * n = grandFatherNeighbors[k]; // This neighbor
                    if(grandFather->getFather() == NULL || n != grandFather->getFather())
                    {
                        grandFatherArrays.push_back(& grandFatherData->getLikelihoodArrayForNeighbor(n->getId())); 
                        grandFatherTProbs.push_back(& pxy_[n->getId()]);
                    }
                }
                //Compute array 1: grand father array
                VVVdouble array1 = *sonArray;
                resetLikelihoodArray(array1);
                grandFatherArrays.push_back(sonArray);
                grandFatherTProbs.push_back(& pxy_[son->getId()]);
                if(grandFather->hasFather())
                {
                    computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, & grandFatherData->getLikelihoodArrayForNeighbor(grandFather->getFather()->getId()), & pxy_[grandFather->getId()], array1, nbGrandFatherNeighbors, nbDistinctSites_, nbClasses_, nbStates_, false); 
                    
                }  
                else
                {
                    computeLikelihoodFromArrays(grandFatherArrays, grandFatherTProbs, array1, nbGrandFatherNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false); 
                    
                    //This is the root node, we have to account for the ancestral frequencies:
                    for(unsigned int i = 0; i < nbDistinctSites_; i++)
                    {
                        for(unsigned int j = 0; j < nbClasses_; j++)
                        {
                            for(unsigned int x = 0; x < nbStates_; x++)
                                array1[i][j][x] *= rootFreqs_[x];
                        }
                    }
                }
                
                
                //Compute array 2: parent array
                VVVdouble array2 = *uncleArray;
                resetLikelihoodArray(array2);
                parentArrays.push_back(uncleArray);
                parentTProbs.push_back(& pxy_[uncle->getId()]);
                computeLikelihoodFromArrays(parentArrays, parentTProbs, array2, nbParentNeighbors + 1, nbDistinctSites_, nbClasses_, nbStates_, false); 
                //Initialize BranchLikelihood:
                brLikFunction_->initModel(nniLk->getModel(), nniLk->getRateDistribution());
                brLikFunction_->initLikelihoods(&array1, &array2);
                ParameterList parameters;	
                unsigned int pos = 0;
                while (pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
                if(pos == nodes_.size()) throw Exception("NNIHomogeneousTreeLikelihood::testNNI. Invalid node id.");
                Parameter brLen = getParameter("BrLen" + TextTools::toString(pos));
                brLen.setName("BrLen");
                parameters.addParameter(brLen);
                brLikFunction_->setParameters(parameters);
                //Re-estimate branch length:
                brentOptimizer_->setMessageHandler(NULL);
                brentOptimizer_->setFunction(brLikFunction_);
                brentOptimizer_->getStopCondition()->setTolerance(0.1);
                brentOptimizer_->setInitialInterval(brLen.getValue(), brLen.getValue()+0.01);
                brentOptimizer_->init(parameters); 
                brentOptimizer_->optimize(); 
                brLenNNIValues_[nodeId] = brentOptimizer_->getParameters().getParameter("BrLen").getValue();
                brLikFunction_->resetLikelihoods(); //Array1 and Array2 will be destroyed after this function call.
                //We should not keep pointers towards them...
                */
                
                //Return the resulting likelihood:
             //   double temp = getValue() ; 
                
                //std::cout<<"temp "<< temp<< "; brLikFunction_->getValue()"<< brLikFunction_->getValue()<<std::endl;
                
          //      double tot = brLikFunction_->getValue() - ScenarioMLValue - temp; // -newsequencelogLk - (newscenariologLk) - (-currenttotallogLk); if <0, worth doing 
           //     if (tot<0) 
                double tot = - ScenarioMLValue + _scenarioLikelihood + newLkMinusOldLk;
               // if (newLkMinusOldLk<0) 
                    if (tot < 0)
                {
                    _tentativeScenarioLikelihood=ScenarioMLValue;
                }
                return tot;
                // std::cout << "after "<<TreeTools::treeToParenthesis (tree_, true)<<std::endl;
                
            }
            else 
            {
                _tentativeMLindex = -1;
                return 1;
            }
        }
        else 
        {
            if  (ScenarioMLValue >  _scenarioLikelihood) 
            {
                _tentativeScenarioLikelihood=ScenarioMLValue;
                return (- ScenarioMLValue + _scenarioLikelihood);
            }
            else 
            {
                _tentativeMLindex = -1;
                return 1;
            }
        }
    }
    else {
        _tentativeMLindex = -1;
        return 1;
    }
}

/*******************************************************************************/

void DLGeneTreeLikelihood::doNNI(int nodeId) throw (NodeException)
{
    //std::cout<<"\t\t\tIN DONNI "<< std::endl;
    //std::cout << TreeTools::treeToParenthesis(*tree_, true) << std::endl;
    //Perform the topological move, the likelihood array will have to be recomputed...
    
    nniLk_->doNNI(nodeId);

    /*
    Node * son    = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->getNode(nodeId);
    if(!son->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'son' must not be the root node.", nodeId);
    Node * parent = son->getFather();
    if(!parent->hasFather()) throw NodeException("DRHomogeneousTreeLikelihood::doNNI(). Node 'parent' must not be the root node.", parent->getId());
    Node * grandFather = parent->getFather();
    //From here: Bifurcation assumed.
    //In case of multifurcation, an arbitrary uncle is chosen.
    //If we are at root node with a trifurcation, this does not matter, since 2 NNI are possible (see doc of the NNISearchable interface).
    unsigned int parentPosition = grandFather->getSonPosition(parent);
    Node * uncle = grandFather->getSon(parentPosition > 1 ? 0 : 1 - parentPosition);
    //Swap nodes:
    parent->removeSon(son);
    grandFather->removeSon(uncle);
    parent->addSon(uncle);
    grandFather->addSon(son);
    unsigned int pos = 0;
    while(pos < nodes_.size() && nodes_[pos]->getId() != parent->getId()) pos++;
    if(pos == nodes_.size()) throw Exception("ReconciliationTreeLikelihood::doNNI. Invalid node id.");
    
    //Julien Proposition :
    
    std::string name = "BrLen"+ TextTools::toString(pos);

    
    if (_considerSequenceLikelihood ) 
    {
        if(brLenNNIValues_.find(nodeId) != brLenNNIValues_.end())
        {
            double length = brLenNNIValues_[nodeId];
            brLenParameters_.setParameterValue(name, length);
            getParameter_(name).setValue(length);
            parent->setDistanceToFather(length);
        }
        else std::cerr << "ERROR, branch not found: " << nodeId << std::endl;
        try { brLenNNIParams_.addParameter(brLenParameters_.getParameter(name)); }
        catch(ParameterException & ex)
        {
            std::cerr << "DEBUG:" << std::endl;
            brLenNNIParams_.printParameters(std::cerr);
            std::cerr << "DEBUG:" << name << std::endl;
        }
        brLenNNIParams_[brLenNNIParams_.size()-1].removeConstraint();
    }
    else 
    {
        double length = 0.1;
        // brLenParameters_.setParameterValue(name, length);
        // getParameter_(name).setValue(length);
        parent->setDistanceToFather(length);
    }
    
    */
    //In case of copy of this object, we must remove the constraint associated to this stored parameter:
    //(It should be also possible to update the pointer in the copy constructor,
    //but we do not need the constraint info here...).
    
    _MLindex = _tentativeMLindex;
    _duplicationNumbers = _tentativeDuplicationNumbers;
    _lossNumbers = _tentativeLossNumbers;
    _branchNumbers = _tentativeBranchNumbers;
    
    _nodesToTryInNNISearch = _tentativeNodesToTryInNNISearch;
    
   // std::cout <<"_tentativeScenarioLikelihood: "<<_tentativeScenarioLikelihood<<std::endl;
    _scenarioLikelihood = _tentativeScenarioLikelihood;// + _brLikFunction->getValue();
    
    //Now we need to update _rootedTree
    TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
    //First we root this temporary tree as in _rootedTree (same lines as in testNNI)
    int id = tree->getRootNode()->getId();
    if(TreeTemplateTools::hasNodeWithId(*(_rootedTree->getRootNode()->getSon(0)),id)) {
        tree->newOutGroup(_rootedTree->getRootNode()->getSon(1)->getId());
    }
    else {
        tree->newOutGroup(_rootedTree->getRootNode()->getSon(0)->getId());
    }
    //Then we root this tree according to MLindex
    tree->newOutGroup(_MLindex);
    //We update _rootedTree
    if (_rootedTree) delete _rootedTree;
    _rootedTree = tree->clone();
    delete tree;
    //we need to update the sequence likelihood
    if (_considerSequenceLikelihood ) 
    {
        OptimizeSequenceLikelihood(true);
    }
    OptimizeReconciliationLikelihood(true);
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getDuplicationNumbers(){
    return _duplicationNumbers;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getLossNumbers(){
    return _lossNumbers;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::getBranchNumbers(){
    return _branchNumbers;
}


/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get0LineagesNumbers() const {
    return _num0Lineages;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get1LineagesNumbers() const {
    return _num1Lineages;
}

/*******************************************************************************/

std::vector <int> DLGeneTreeLikelihood::get2LineagesNumbers() const {
    return _num2Lineages;
}

/*******************************************************************************/

void DLGeneTreeLikelihood::setProbabilities(std::vector <double> duplicationProbabilities, std::vector <double> lossProbabilities){
    
    _lossProbabilities = lossProbabilities;
    _duplicationProbabilities = duplicationProbabilities; 
}

/*******************************************************************************/

int DLGeneTreeLikelihood::getRootNodeindex(){
    return _MLindex;
}

/*******************************************************************************/
/*
void DLGeneTreeLikelihood::resetSequenceLikelihood(){
    _sequenceLikelihood = UNLIKELY;
}
*/
/*******************************************************************************/

double DLGeneTreeLikelihood::getSequenceLikelihood() {
    return nniLk_->getValue(); 
}
/*******************************************************************************/

void DLGeneTreeLikelihood::initialize() {
    nniLk_->initialize();
    return;
}



/*******************************************************************************/
void DLGeneTreeLikelihood::print () const {
    std::cout << "Species tree:"<<std::endl;
    std::cout << TreeTools::treeToParenthesis (getSpTree(), true)<<std::endl;
    std::cout << "Gene family rooted tree:"<<std::endl;
    std::cout << TreeTools::treeToParenthesis (getRootedTree(), true)<<std::endl;
    std::cout << "Gene family tree:"<<std::endl;
    std::cout << TreeTools::treeToParenthesis (nniLk_->getTree(), true)<<std::endl;
    std::cout << "0 lineage numbers"<<std::endl;
    VectorTools::print(get0LineagesNumbers());
    std::cout << "1 lineage numbers"<<std::endl;
    VectorTools::print(get1LineagesNumbers());
    std::cout << "2 lineages numbers"<<std::endl;
    VectorTools::print(get2LineagesNumbers());
    std::cout << "Expected numbers of losses"<<std::endl;
    VectorTools::print(_lossProbabilities);
    std::cout << "Expected numbers of duplications"<<std::endl;
    VectorTools::print(_duplicationProbabilities);
    std::cout << "Root index"<<std::endl;
    std::cout << _MLindex <<std::endl;
    
    
}


/************************************************************************
 * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
 * and executes the ones with the highest likelihood. 
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRs(map<string, string> params) {
    /*  std::cout <<"\t\t\tStarting MLSearch : current tree : "<< std::endl;
     std::cout<< TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;*/
    if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
        //We don't do SPRs
        //std::cout << "WE DONT DO SPRS"<<std::endl;
        computeReconciliationLikelihood();
        return;
    }
   // std::cout << "WE DO SPRS"<<std::endl;

    std::vector <int> nodeIdsToRegraft;
    bool betterTree;
    TreeTemplate<Node> * treeForSPR = 0;
    TreeTemplate<Node> * bestTree = 0;
    // if (getLogLikelihood()==UNLIKELY) 
    computeReconciliationLikelihood();
    //    ReconciliationTreeLikelihood * bestTreeLogLk = this->clone();
    double logL = getLogLikelihood();
    //std::cout << "logL before mapping: " <<getSequenceLikelihood()<<std::endl;
    //initial perturbation of branch lengths, to use suboptimal ones for speed
    /*ParameterList bls = nniLk_->getBranchLengthsParameters () ;
     for (unsigned int i  = 0 ; i < bls.size()/3 ; i++) {
     bls[i].setValue(0.1);
     }
     nniLk_->matchParametersValues(bls);*/
    //  std::cout << "logL before mapping 2: " <<getSequenceLikelihood()<<std::endl;    
    //TEST 14 03 2012, if it works, it's inefficient

    NNIHomogeneousTreeLikelihood * drlk = 0;
    drlk = nniLk_->clone();
   // optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (nniLk_), 0.1, params);
    
    int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
    {
        params[ std::string("optimization.max_number_f_eval")] = 100;
    }
    PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
    params[ std::string("optimization.max_number_f_eval")] = backup;
    
    delete nniLk_;
    nniLk_ = drlk->clone();
    
    delete drlk;
    drlk = 0;
    
    /*  optimizeBLForSPRs(dynamic_cast<DRTreeLikelihood*> (nniLk_),
     0.1, params);*/
        
    
    
    // std::cout<< "Unrooted starting tree: "<<TreeTools::treeToParenthesis(nniLk_->getTree(), true)<< std::endl;
    /*
     bls = nniLk_->getBranchLengthsParameters () ;
     for (unsigned int i  = 0 ; i < bls.size() ; i++) {
     if (_rootedTree->getNode(i)->hasFather()) {
     _rootedTree->getNode(i)->setDistanceToFather(bls[i].getValue());
     }
     }*/
    
    logL = getLogLikelihood();
    //  std::cout << "logL after mapping: " <<getSequenceLikelihood()<<std::endl;
    
    double bestlogL = logL;
    double candidateScenarioLk ;
    double bestSequenceLogL = getSequenceLikelihood();
    double bestScenarioLk = getScenarioLikelihood();
    // std::cout << "LOGL: "<<logL << "ScenarioLK: "<< bestScenarioLk <<"; sequenceLK: "<<getSequenceLikelihood() << std::endl;
    int numIterationsWithoutImprovement = 0;
    NNIHomogeneousTreeLikelihood * bestDrlk = 0;
    breadthFirstreNumber (*_rootedTree);
    
    
    // DRHomogeneousTreeLikelihood drlk;
    //  std::cout<< "Starting tree: "<<TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    
    string parentDup;
    string nodeDup;
    string numLoss = "0";
    
    bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    
    
    while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    {
        
        annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                               *_rootedTree, 
                                               _rootedTree->getRootNode(), 
                                               _seqSp, _spId); 
        
        for (int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
        {
            Node * n = _rootedTree->getNode(nodeForSPR);
            /*            if (n->getFather()->hasBranchProperty("Ev")) {
             parentDup = (dynamic_cast<const BppString *>(n->getFather()->getBranchProperty("Ev")))->toSTL() ;
             }*/
            if (n->hasBranchProperty("L")) {
                numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
            }
            /* else {
             parentDup = "F" ;}*/
            //   std::cout << "parentDup : "<< parentDup <<std::endl;
            /*  if ( _rootedTree->getNode(nodeForSPR)->hasNodeProperty("D")) { 
             nodeDup = (dynamic_cast<const BppString *>(_rootedTree->getNode(nodeForSPR)->getNodeProperty("D")))->toSTL() ; }
             else {
             nodeDup = "F" ; }*/
            //if ( parentDup == "D"  /*|| nodeDup == "Y"*/) {
            if ( numLoss != "0"  /*|| nodeDup == "Y"*/) {
                
                //CHANGE12 11 2011 buildVectorOfRegraftingNodesLimitedDistance(*_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                buildVectorOfRegraftingNodesGeneTree(*_spTree, *_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                
                /*              
                 if (_rootedTree->getNode(nodeForSPR)->isLeaf()) {
                 std::cout << "nodeForSPR: "<< _rootedTree->getNode(nodeForSPR)->getName()<<" ; "<< nodeForSPR << "; Node ids to regraft: " <<std::endl;
                 }
                 else {
                 std::cout << "nodeForSPR: "<< nodeForSPR << "; Node ids to regraft: " <<std::endl;
                 }
                 VectorTools::print(nodeIdsToRegraft);
                 */
                betterTree = false;
                for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
                {
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    treeForSPR = _rootedTree->clone();
                    
                    makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false);
                                        
                    //breadthFirstreNumber (*treeForSPR);
                    
                    //Compute the DL likelihood
                    candidateScenarioLk =  findMLReconciliationDR (_spTree, treeForSPR, 
                                                                   _seqSp, _spId, 
                                                                   _lossProbabilities, 
                                                                   _duplicationProbabilities, 
                                                                   _tentativeMLindex, 
                                                                   _tentativeNum0Lineages, 
                                                                   _tentativeNum1Lineages, 
                                                                   _tentativeNum2Lineages, 
                                                                   _tentativeNodesToTryInNNISearch, false); 
                    // std::cout<< "candidate tree: "<< TreeTools::treeToParenthesis(*treeForSPR, true)<< std::endl;
                    
                    if (candidateScenarioLk > bestScenarioLk)// - 0.1) //We investigate the sequence likelihood if the DL likelihood is not bad
                    {
                        
                        if (computeSequenceLikelihoodForSPR) {
                            
                            // std::cout << "good candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                            if (drlk) {
                                delete drlk;
                                drlk = 0;
                            }
                           // std::cout << "COMPUTING SEQLK: " << TreeTools::treeToParenthesis(*treeForSPR, true) <<std::endl;
                            drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                            //  drlk  = new DRTreeLikelihood (*treeForSPR, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                            
                            
                            //drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                            /*          ParameterList bls = nniLk_->getBranchLengthsParameters () ;
                             for (unsigned int i  = 0 ; i < bls.size()/3 ; i++) {
                             bls[i].setValue(0.1);
                             }
                             drlk->matchParametersValues(bls);
                             */
                            drlk->initialize();
                            /*std::cout << "TEST getValue: "<< drlk->getValue() <<std::endl;
                            std::cout<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                           
                            //  std::cout << "Good SPR; Sequence lk before optimization: "<< -drlk->getValue() << std::endl;
                            //  _rootedTree = treeForSPR->clone();
                            // makeSPR(*tree_, nodeForSPR, nodeIdsToRegraft[i]);
                            // tree_->unroot();
                            //  tree_ = treeForSPR->clone();
                            
                            //Compute sequence likelihood and improve branch lengths on the gene tree.
                            /*  std::cout<< TreeTools::treeToParenthesis(*tree_, true)<< std::endl;
                             optimizeBLMapping(this,
                             0.1);
                             */
                            // tree_->unroot();
                            //fireParameterChanged(getParameters());
                            //Compute sequence likelihood and improve branch lengths on the gene tree.
                            //  std::cout<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;
                            
                            //  params[ std::string("optimization.topology")] = "false";
                                                        
                            /*TEMP18012012 */
                            /*optimizeBLMappingForSPRs does not work for some families...
                            optimizeBLMappingForSPRs(dynamic_cast<DRTreeLikelihood*> (drlk),
                                                     0.1, params);*/
                            /*                            optimizeBLForSPRs(dynamic_cast<DRTreeLikelihood*> (drlk),
                             0.1, params);
                             */  
                           /* std::cout << "before norm opt: " <<std::endl;
                            std::cout<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            
                            //TEST 13 03 2012
                            int backup = ApplicationTools::getIntParameter("optimization.max_number_f_eval", params, false, "", true, false);
                            {
                                params[ std::string("optimization.max_number_f_eval")] = 100;
                            }
                            
                            PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
                            params[ std::string("optimization.max_number_f_eval")] = backup;

                           /* std::cout << "TEST getValue2: "<< drlk->getValue() <<std::endl;
                            std::cout<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            
                            /*std::cout << "Good SPR; Sequence lk after optimization: "<< -drlk->getValue() << std::endl;
                             std::cout<<"Good SPR; optimized tree: "<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;*/
                            logL = candidateScenarioLk - drlk->getValue();
                        }
                        else {
                            logL = candidateScenarioLk - bestSequenceLogL;
                        }
                        
                    }
                    else { 
                        // std::cout << "bad candidateScenarioLk: "<< candidateScenarioLk<<"; bestScenarioLk: "<< bestScenarioLk<< std::endl;
                        logL =logL - 10;
                    }
                    //If the candidate tree has a DL + sequence Lk better than the current best
                    if (logL - 0.01 > bestlogL) 
                    {
                        betterTree = true;
                        bestlogL =logL;
                        bestScenarioLk = candidateScenarioLk;
                        if (computeSequenceLikelihoodForSPR) {
                            bestSequenceLogL = drlk->getValue();
                            
                            if (bestDrlk) {
                                delete bestDrlk;
                                bestDrlk = 0;
                            }
                            
                            bestDrlk = drlk->clone();
                        }
                        if (bestTree) {
                            delete bestTree;
                            bestTree = 0;
                        }
                        
                        bestTree = treeForSPR->clone();
                        /*  std::cout << "\t\t\tSPRs: Better candidate tree likelihood : "<<bestlogL<< std::endl;
                         std::cout << "\t\t\tTotal likelihood: "<<logL <<"; Reconciliation likelihood: "<< bestScenarioLk << ", Sequence likelihood: "<< bestSequenceLogL <<", for tree: "<< std::endl;*/
                        /*TreeTools::treeToParenthesis(bestTreeLogLk->getRootedTree()) << */
                    }
                    else {
                        //   copyContentsFrom(*bestTreeLogLk);
                        //  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
                    }
                }
                if (betterTree) //If, among all the SPRs tried, a better tree has been found 
                {
                    
                    logL = bestlogL; 
                    numIterationsWithoutImprovement = 0;
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    if (_rootedTree) 
                    {
                        delete _rootedTree;
                        _rootedTree = 0;
                    }
                    _rootedTree = bestTree->clone();
                    breadthFirstreNumber (*_rootedTree);
                    
                    if (bestTree) {
                        delete bestTree;
                        bestTree = 0;
                    }
                    
                    /* bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(bestDrlk->getTree()))->clone();
                     vector<int> drlkNodes = bestTree->getNodesId();
                     
                     for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
                     if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
                     _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
                     }*/
                    
                    /*    std::cout <<"\t\t\tSPRs: Improvement! : "<<numIterationsWithoutImprovement<< std::endl;
                     std::cout << "\t\t\tNew total Likelihood value "<<logL<< std::endl;
                     std::cout <<"\t\t\tNumber of gene trees tried : "<<index<< std::endl;*/
                    /*  if (bestTree) {
                     delete bestTree;
                     bestTree = 0;
                     }*/
                    if (computeSequenceLikelihoodForSPR) {                        
                        if (nniLk_) {
                            delete nniLk_;
                            nniLk_ = 0;
                        }
                        nniLk_ = bestDrlk->clone();                        
                    }
                    /* std::cout << "NOT ROOTED: "<<std::endl;
                     std::cout<< TreeTools::treeToParenthesis(drlk->getTree(), true)<< std::endl;
                     std::cout << "ROOTED: "<<std::endl;
                     std::cout<< TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
                     */
                }
                else 
                {
                    // logL = bestlogL;  
                    numIterationsWithoutImprovement++;
                    //  std::cout <<"\t\t\tSPRs: Number of iterations without improvement : "<<numIterationsWithoutImprovement << "; Total number  of iterations: "<< index << std::endl;
                }
                
                if (treeForSPR) 
                {
                    delete treeForSPR;
                    treeForSPR = 0;
                }
                if (bestTree) {
                    delete bestTree;
                    bestTree = 0;
                }
                if (bestDrlk) {
                    delete bestDrlk;
                    bestDrlk = 0;
                }
                if (drlk) {
                    delete drlk;
                    drlk = 0;
                }
            }
            else {
                numIterationsWithoutImprovement++;
            }
        }
    }
    
    _scenarioLikelihood = bestScenarioLk;
    //final branch lengths optimization
    
    _rootedTree->resetNodesId();
    if (drlk) {
        delete drlk;
        drlk = 0;
    }

    bestTree = _rootedTree->clone();

    drlk = new NNIHomogeneousTreeLikelihood (*bestTree, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
    drlk->initialize();
    
    PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
  
    if (nniLk_) {
        delete nniLk_;
        nniLk_ = 0;
    }
    
    nniLk_ = drlk->clone();
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
    
    vector<int> drlkNodes = bestTree->getNodesId();
    for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
        if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
            _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
    }
    
    //One more reconciliation, to update the "_num*Lineages" vectors.
    computeReconciliationLikelihood();
        
    Nhx *nhx = new Nhx();
    annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                           *_rootedTree, 
                                           _rootedTree->getRootNode(), 
                                           _seqSp, _spId); 
    cout << "Reconciled tree: "<<endl;
    nhx->write(*_rootedTree, cout);
    /*  cout << "unrooted tree: "<<endl;
     nhx->write(*bestTree, cout);*/
    //    std::cout<< TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    if (nhx) delete nhx;
    
}





/************************************************************************
 * Tries all possible (and potentially valuable) SPRs and computes the DL likelihood.
 * For the three best ones, if the DL lk is better than the current best, computes the sequence lk.
 * Repeat until no SPR can improve the gene tree. 
 * This algorithm is slower and less efficient than the other one.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeSPRs2(map<string, string> params) { 
    /*  std::cout <<"\t\t\tStarting MLSearch : current tree : "<< std::endl;
     std::cout<< TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;*/
    if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
        //We don't do SPRs
      //  std::cout << "WE DONT DO SPRS"<<std::endl;
        computeReconciliationLikelihood();
        return;
    }
  //  std::cout << "WE DO SPRS"<<std::endl;
    
    std::vector <int> nodeIdsToRegraft;
    bool betterTree = false;
    TreeTemplate<Node> * treeForSPR = 0;
    TreeTemplate<Node> * bestTree = 0;
    // if (getLogLikelihood()==UNLIKELY) 
    computeReconciliationLikelihood();
    //    ReconciliationTreeLikelihood * bestTreeLogLk = this->clone();
    double logL = getLogLikelihood();

    std::map <double, std::vector<int> > allDLLKsAndIndices;
    std::vector< TreeTemplate<Node>* > allTreesTested;
    unsigned int index = 0;
    
    optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (nniLk_),
                             0.1, params);
    logL = getLogLikelihood();
    
    double bestlogL = logL;
    double candidateScenarioLk ;
    double bestSequenceLogL = getSequenceLikelihood();
    double bestScenarioLk = getScenarioLikelihood();
    // std::cout << "bestLogL: "<< bestlogL << "bestScenarioLK: "<< bestScenarioLk <<"; bestSequenceLK: "<<bestSequenceLogL << std::endl;
    int numIterationsWithoutImprovement = 0;
    NNIHomogeneousTreeLikelihood * drlk = 0;
    NNIHomogeneousTreeLikelihood * bestDrlk = 0;
    breadthFirstreNumber (*_rootedTree);
    
    
    // DRHomogeneousTreeLikelihood drlk;
    //  std::cout<< "Starting tree: "<<TreeTools::treeToParenthesis(*_rootedTree, true)<< std::endl;
    
    string parentDup;
    string nodeDup;
    string numLoss = "0";
    
    bool computeSequenceLikelihoodForSPR = ApplicationTools::getBooleanParameter("compute.sequence.likelihood.in.sprs", params, true, "", false, false);
    
    
    while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    {
        betterTree=false;
        allDLLKsAndIndices.clear();
        allTreesTested.clear();        
        index = 0;
        annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                               *_rootedTree, 
                                               _rootedTree->getRootNode(), 
                                               _seqSp, _spId); 
        for (unsigned int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--) 
        {
            Node * n = _rootedTree->getNode(nodeForSPR);
            /*            if (n->getFather()->hasBranchProperty("Ev")) {
             parentDup = (dynamic_cast<const BppString *>(n->getFather()->getBranchProperty("Ev")))->toSTL() ;
             }*/
            if (n->hasBranchProperty("L")) {
                numLoss = (dynamic_cast<const BppString *>(n->getBranchProperty("L")))->toSTL() ;
            }
            if ( numLoss != "0" ) {
                buildVectorOfRegraftingNodesGeneTree(*_spTree, *_rootedTree, nodeForSPR, sprLimit_, nodeIdsToRegraft);
                for (unsigned int i =0 ; i<nodeIdsToRegraft.size() ; i++) 
                {
                    if (treeForSPR) 
                    {
                        delete treeForSPR;
                        treeForSPR = 0;
                    }
                    treeForSPR = _rootedTree->clone();
                    
                    makeSPR(*treeForSPR, nodeForSPR, nodeIdsToRegraft[i], false);
                    
                    //breadthFirstreNumber (*treeForSPR);
                    
                    //Compute the DL likelihood
                    candidateScenarioLk =  findMLReconciliationDR (_spTree, treeForSPR, 
                                                                   _seqSp, _spId, 
                                                                   _lossProbabilities, 
                                                                   _duplicationProbabilities, 
                                                                   _tentativeMLindex, 
                                                                   _tentativeNum0Lineages, 
                                                                   _tentativeNum1Lineages, 
                                                                   _tentativeNum2Lineages, 
                                                                   _tentativeNodesToTryInNNISearch, false); 
                    if (allDLLKsAndIndices.find(- candidateScenarioLk) != allDLLKsAndIndices.end() ) {
                        allDLLKsAndIndices[- candidateScenarioLk].push_back(index);
                    }
                    else {
                        allDLLKsAndIndices[- candidateScenarioLk] = vector<int> ();
                        allDLLKsAndIndices[- candidateScenarioLk].push_back(index);
                    }
                    allTreesTested.push_back (treeForSPR->clone());
                    index = index +1;
                   // logL =logL - 10;
                }
            } // End if ( numLoss != "0" )
            else {
                numIterationsWithoutImprovement++;
            }
        } //End for (int nodeForSPR=_rootedTree->getNumberOfNodes()-1 ; nodeForSPR >0; nodeForSPR--)
        
        int numberOfTopologiesToTry = 3;
        
        vector <int> treeIndicesToTry;
        vector <double> scenarioLksToTry;
        for(std::map<double, std::vector<int> >::iterator it = allDLLKsAndIndices.begin(); it != allDLLKsAndIndices.end(); it++){
           // std::cout << "Tested scenario Lk: "<< it->first <<" Found for "<<it->second.size()<<" topologies"<<std::endl;
            if (- it->first > bestScenarioLk) {
                for (unsigned int i = 0 ; i < it->second.size(); i++) {
                    scenarioLksToTry.push_back(- it->first);
                    treeIndicesToTry.push_back(it->second[i]);
                }
                if (treeIndicesToTry.size() >= numberOfTopologiesToTry) break;
            }
        }
       // std::cout << "treeIndicesToTry.size(): " << treeIndicesToTry.size() <<std::endl;
        /////For the top k candidate topologies, we compute the loglk
        for (unsigned int k = 0 ; k < treeIndicesToTry.size() ; k++)
        {
            if (treeForSPR) 
            {
                delete treeForSPR;
                treeForSPR = 0;
            }
            treeForSPR = allTreesTested[treeIndicesToTry[k]]->clone();
            candidateScenarioLk = scenarioLksToTry[k];
            
            
            if (computeSequenceLikelihoodForSPR) {
                if (drlk) {
                    delete drlk;
                    drlk = 0;
                }
               // std::cout << "COMPUTING SEQLK: " << TreeTools::treeToParenthesis(*treeForSPR, true) <<std::endl;
                drlk  = new NNIHomogeneousTreeLikelihood (*treeForSPR, *(nniLk_->getData()), nniLk_->getSubstitutionModel(), nniLk_->getRateDistribution(), true, false);
                drlk->initialize();
                optimizeBLMappingForSPRs(dynamic_cast<DRHomogeneousTreeLikelihood*> (drlk),
                                         0.1, params);
                logL = candidateScenarioLk - drlk->getValue();
            }
            else {
                logL = candidateScenarioLk - bestSequenceLogL;
            }
            //If the candidate tree has a DL + sequence Lk better than the current best
            if (logL - 0.01 > bestlogL) 
            {
                betterTree = true;
                bestlogL =logL;
                bestScenarioLk = candidateScenarioLk;
                if (computeSequenceLikelihoodForSPR) {
                    bestSequenceLogL = drlk->getValue();
                    
                    if (bestDrlk) {
                        delete bestDrlk;
                        bestDrlk = 0;
                    }
                    
                    bestDrlk = drlk->clone();
                }
                if (bestTree) {
                    delete bestTree;
                    bestTree = 0;
                }
                
                bestTree = treeForSPR->clone();
            }
            else {
                //   copyContentsFrom(*bestTreeLogLk);
                //  std::cout << "\t\t\tSPRs: No improvement : "<< logL << " compared to current best: "<< bestlogL << std::endl;
            }
        }
        
        if (betterTree) //If, among all the SPRs tried, a better tree has been found 
        {
            logL = bestlogL; 
            numIterationsWithoutImprovement = 0;
            if (treeForSPR) 
            {
                delete treeForSPR;
                treeForSPR = 0;
            }
            if (_rootedTree) 
            {
                delete _rootedTree;
                _rootedTree = 0;
            }
            _rootedTree = bestTree->clone();
            breadthFirstreNumber (*_rootedTree);
            
            if (bestTree) {
                delete bestTree;
                bestTree = 0;
            }
            if (computeSequenceLikelihoodForSPR) {                        
                if (nniLk_) {
                    delete nniLk_;
                    nniLk_ = 0;
                }
                nniLk_ = bestDrlk->clone();                        
            }
        }
        else 
        {
            numIterationsWithoutImprovement = _rootedTree->getNumberOfNodes() - 2;
        }
        
        if (treeForSPR) 
        {
            delete treeForSPR;
            treeForSPR = 0;
        }
        if (bestTree) {
            delete bestTree;
            bestTree = 0;
        }
        if (bestDrlk) {
            delete bestDrlk;
            bestDrlk = 0;
        }
        if (drlk) {
            delete drlk;
            drlk = 0;
        }
    }// End while (numIterationsWithoutImprovement < _rootedTree->getNumberOfNodes() - 2)
    
    _scenarioLikelihood = bestScenarioLk;
    //final branch lengths optimization
    
    _rootedTree->resetNodesId();
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    
    bestTree = _rootedTree->clone();
    
    drlk = new NNIHomogeneousTreeLikelihood (*bestTree, 
                                             *(nniLk_->getData()), 
                                             nniLk_->getSubstitutionModel(), 
                                             nniLk_->getRateDistribution(), 
                                             true, false);
    drlk->initialize();
    
    PhylogeneticsApplicationTools::optimizeParameters(drlk, drlk->getParameters(), params, "", true, false);
    
    if (nniLk_) {
        delete nniLk_;
        nniLk_ = 0;
    }
    
    nniLk_ = drlk->clone();
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    bestTree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
    
    vector<int> drlkNodes = bestTree->getNodesId();
    for (unsigned int j = 0 ; j < drlkNodes.size() ; j++) {
        if (_rootedTree->getNode(drlkNodes[j])->hasFather() && bestTree->getNode(drlkNodes[j])->hasFather())
            _rootedTree->getNode(drlkNodes[j])->setDistanceToFather(bestTree->getNode(drlkNodes[j])->getDistanceToFather());
    }
    
    //One more reconciliation, to update the "_num*Lineages" vectors.
    
    computeReconciliationLikelihood();
    
    Nhx *nhx = new Nhx();
    annotateGeneTreeWithDuplicationEvents (*_spTree, 
                                           *_rootedTree, 
                                           _rootedTree->getRootNode(), 
                                           _seqSp, _spId); 
    cout << "Reconciled tree: "<<endl;
    nhx->write(*_rootedTree, cout);
    if (bestTree) {
        delete bestTree;
        bestTree = 0;
    }
    if (drlk) {
        delete drlk;
        drlk = 0;
    }
    if (nhx) delete nhx;
    
}






/************************************************************************
 * Tries all NNIs, and accepts NNIs that improve the likelihood as soon as
 * they have been tried.
 ************************************************************************/
void DLGeneTreeLikelihood::refineGeneTreeNNIs(map<string, string> params, unsigned int verbose ) {
    if (ApplicationTools::getBooleanParameter("optimization.topology", params, true, "", false, false) == false ) {
        //We don't do NNIs
        computeReconciliationLikelihood();
        return;
    }
    bool test = true;
    do
    { 
        TreeTemplate<Node> * tree = dynamic_cast<const TreeTemplate<Node> *> (&(nniLk_->getTree()))->clone();
        vector<Node *> nodes = tree->getNodes();
                
        vector<Node *> nodesSub = nodes;
        for(unsigned int i = nodesSub.size(); i>0; i--)
        {
            // !!! must not reach i==0 because of unsigned int
            if(!(nodesSub[i-1]->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove root node.      
            else if(!(nodesSub[i-1]->getFather()->hasFather())) nodesSub.erase(nodesSub.begin()+i-1);//Remove son of root node. 
        }
        
        // Test all NNIs:
        test = false;
        for(unsigned int i = 0; !test && i < nodesSub.size(); i++)
        {
            Node* node = nodesSub[i];
            double diff = testNNI(node->getId());
            if(verbose >= 3)
            {
                ApplicationTools::displayResult("   Testing node " + TextTools::toString(node->getId())
                                                + " at " + TextTools::toString(node->getFather()->getId()),
                                                TextTools::toString(diff));
            }
            
            if(diff < 0.)
            { //Good NNI found...
                if(verbose >= 2)
                {
                    ApplicationTools::displayResult("   Swapping node " + TextTools::toString(node->getId())
                                                    + " at " + TextTools::toString(node->getFather()->getId()),
                                                    TextTools::toString(diff));
                }
                doNNI(node->getId());
                test = true;
                nniLk_->topologyChangeTested( 	TopologyChangeEvent ());
                nniLk_->topologyChangeSuccessful( 	TopologyChangeEvent ());
                if(verbose >= 1)
                    ApplicationTools::displayResult("   Current value", TextTools::toString(getValue(), 10));
            }
        }
    }
    while(test);
}


/************************************************************************
 * Tells if the gene family is single copy (1 gene per sp)
 ************************************************************************/
bool DLGeneTreeLikelihood::isSingleCopy() {
    vector<string> names = TreeTemplateTools::getLeavesNames(*(_geneTreeWithSpNames->getRootNode() ) );
    vector<string> uniqueNames = VectorTools::unique(names);
    if (uniqueNames.size() == names.size() ) {
        return true;
    }
    return false;
}





