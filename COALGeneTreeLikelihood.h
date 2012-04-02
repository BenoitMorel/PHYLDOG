//
// File: COALGeneTreeLikelihood.h
// Created by: Bastien Boussau 
// Created on: Tue October 04 14:16 2011
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


#ifndef COALGeneTreeLikelihood_h
#define COALGeneTreeLikelihood_h


#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>

// From NumCalc:
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Parametrizable.h>

#include "ReconciliationTools.h"
#include "COALTools.h"

#include "GeneTreeAlgorithms.h"
#include "mpi.h" 

namespace bpp 
{
    
    /**
     * @brief This class adds support for coalescence-based reconciliation to a species tree to the NNIHomogeneousTreeLikelihood class.
     */
    class COALGeneTreeLikelihood
    {
        NNIHomogeneousTreeLikelihood * nniLk_;
        //  TreeTemplate<Node> * _tree;
        TreeTemplate<Node> * _spTree;
        TreeTemplate<Node> * _rootedTree;
        TreeTemplate<Node> * _geneTreeWithSpNames;
        const std::map <std::string, std::string> _seqSp;
        std::map <std::string, int> _spId;
 
        //coalCounts: vector of genetreenbnodes vectors of 3 (3 directions) vectors of sptreenbnodes vectors of 2 ints
        std::vector < std::vector<std::vector<unsigned int> > > coalCounts;
        //coalBl: length of a branch of the species tree, in coalescent units (1 coalescent unit = N generations)
        std::vector < double > coalBl;
        
        /*
        std::vector <int> _duplicationNumbers;
        std::vector <int> _lossNumbers;
        std::vector <int>  _branchNumbers;
        std::vector <double> _duplicationProbabilities;
        std::vector <double> _lossProbabilities; 
        std::vector <int> _num0Lineages;
        std::vector <int> _num1Lineages;
        std::vector <int> _num2Lineages;
        mutable std::vector <int> _tentativeDuplicationNumbers;
        mutable std::vector <int> _tentativeLossNumbers; 
        mutable std::vector <int> _tentativeBranchNumbers; 
        mutable std::vector <int> _tentativeNum0Lineages;
        mutable std::vector <int> _tentativeNum1Lineages; 
        mutable std::vector <int> _tentativeNum2Lineages;
        mutable bool _DLStartingGeneTree;
         */
        
        std::set <int> _nodesToTryInNNISearch;
        double _scenarioLikelihood;
        //  mutable double _sequenceLikelihood;
        int _MLindex;
        bool _rootOptimization;
        mutable std::set <int> _tentativeNodesToTryInNNISearch;
        mutable int _tentativeMLindex;
        mutable double _tentativeScenarioLikelihood;
        mutable int _totalIterations;
        mutable int _counter;
        mutable std::vector <int> _listOfPreviousRoots;
        int _speciesIdLimitForRootPosition;
        int _heuristicsLevel;
        mutable bool _optimizeSequenceLikelihood;
        mutable bool _optimizeReconciliationLikelihood;
        mutable bool _considerSequenceLikelihood;
        unsigned int sprLimit_;
        
    public:
        /**
         * @brief Build a new ReconciliationTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param spTree The species tree
         * @param rootedTree rooted version of the gene tree
         * @param seqSp link between sequence and species names
         * @param spId link between species name and species ID
         * @param coalCounts vector to store coalescent numbers per branch
         * @param coalBl vector to give number of coalescent units per branch of the species tree
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        COALGeneTreeLikelihood(
                               const Tree & tree,
                               SubstitutionModel * model,
                               DiscreteDistribution * rDist,
                               TreeTemplate<Node> & spTree,  
                               TreeTemplate<Node> & rootedTree, 
                               TreeTemplate<Node> & geneTreeWithSpNames,
                               const std::map <std::string, std::string> seqSp,
                               std::map <std::string,int> spId,
                               std::vector < std::vector<std::vector<unsigned int> > > coalCounts,
                               std::vector < double > coalBl,
                               int speciesIdLimitForRootPosition,
                               int heuristicsLevel,
                               int & MLindex, 
                               bool checkRooted = true,
                               bool verbose = false,
                               bool rootOptimization = false, 
                               bool considerSequenceLikelihood = true, 
                               unsigned int sprLimit = 2)
        throw (Exception);
        
        /**
         * @brief Build a new ReconciliationTreeLikelihood object.
         *
         * @param tree The tree to use.
         * @param data Sequences to use.
         * @param model The substitution model to use.
         * @param rDist The rate across sites distribution to use.
         * @param spTree The species tree
         * @param rootedTree rooted version of the gene tree
         * @param seqSp link between sequence and species names
         * @param spId link between species name and species ID
         * @param coalCounts vector to store coalescent numbers per branch
         * @param coalBl vector to give number of coalescent units per branch of the species tree
         * @param speciesIdLimitForRootPosition limit for gene tree rooting heuristics
         * @param heuristicsLevel type of heuristics used
         * @param MLindex ML rooting position     
         * @param checkRooted Tell if we have to check for the tree to be unrooted.
         * If true, any rooted tree will be unrooted before likelihood computation.
         * @param verbose Should I display some info?
         * @throw Exception in an error occured.
         */
        COALGeneTreeLikelihood(
                               const Tree & tree,
                               const SiteContainer & data,
                               SubstitutionModel * model,
                               DiscreteDistribution * rDist,
                               TreeTemplate<Node> & spTree,  
                               TreeTemplate<Node> & rootedTree,  
                               TreeTemplate<Node> & geneTreeWithSpNames,
                               const std::map <std::string, std::string> seqSp,
                               std::map <std::string,int> spId,
                               std::vector < std::vector<std::vector<unsigned int> > > coalCounts,
                               std::vector < double > coalBl,
                               int speciesIdLimitForRootPosition,  
                               int heuristicsLevel,
                               int & MLindex, 
                               bool checkRooted = true,
                               bool verbose = false, 
                               bool rootOptimization = false, 
                               bool considerSequenceLikelihood = true, 
                               unsigned int sprLimit = 2)
        throw (Exception);
        
        /**
         * @brief Copy constructor.
         */ 
        COALGeneTreeLikelihood(const COALGeneTreeLikelihood & lik);
        
        COALGeneTreeLikelihood & operator=(const COALGeneTreeLikelihood & lik);
        
        virtual ~COALGeneTreeLikelihood();
        
        
        
#ifndef NO_VIRTUAL_COV
        COALGeneTreeLikelihood*
#else
        Clonable*
#endif
        clone() const { return new COALGeneTreeLikelihood(*this); }
        
        void initParameters();
        void resetMLindex() ;
        /**
         * @name The NNISearchable interface.
         *
         * Current implementation:
         * When testing a particular NNI, only the branch length of the parent node is optimized (and roughly).
         * All other parameters (substitution model, rate distribution and other branch length are kept at there current value.
         * When performing a NNI, only the topology change is performed.
         * This is up to the user to re-initialize the underlying likelihood data to match the new topology.
         * Usually, this is achieved by calling the topologyChangePerformed() method, which call the reInit() method of the LikelihoodData object.
         * @{
         */
        
        //double getLikelihood() const;
        
        double getLogLikelihood() const;
        
        void computeSequenceLikelihood();
        
        void computeReconciliationLikelihood();
        
        void computeTreeLikelihood();
        
        double getValue() const throw (Exception);
        
        void fireParameterChanged(const ParameterList & params);
        
        double getTopologyValue() const throw (Exception) { return getValue(); } 
        
        double getScenarioLikelihood() const throw (Exception) { return _scenarioLikelihood; }
        
        void setSpTree(TreeTemplate<Node> & spTree) { if (_spTree) delete _spTree; _spTree = spTree.clone(); }
        
        void setSpId(std::map <std::string, int> & spId) {_spId = spId;}
        
        double testNNI(int nodeId) const throw (NodeException);
        
        void doNNI(int nodeId) throw (NodeException);
        
        std::vector < std::vector<std::vector<unsigned int> > > getCoalNumbers() const;
        
        ParameterList getParameters() {return nniLk_->getParameters();}
        
        TreeTemplate<Node> & getSpTree() const {return *_spTree;}
        
        TreeTemplate<Node> & getRootedTree() const {return *_rootedTree;}
        
        TreeTemplate<Node> & getGeneTreeWithSpNames() const {return *_geneTreeWithSpNames;}
        
        std::map <std::string, std::string> getSeqSp() {return _seqSp;}
        
        void setCoalBl (std::vector < double > coalBl);
        
        int getRootNodeindex();
        
        //void resetSequenceLikelihood();
        
        double getSequenceLikelihood();
        
        void OptimizeSequenceLikelihood(bool yesOrNo) const  {
            _optimizeSequenceLikelihood = yesOrNo;
        }
        
        void OptimizeReconciliationLikelihood(bool yesOrNo) const {
            _optimizeReconciliationLikelihood = yesOrNo;
        }
        
        void initialize();
        
        void print() const;
        
        /************************************************************************
         * Tries all SPRs at a distance < dist for all possible subtrees of the subtree starting in node nodeForSPR, 
         * and executes the ones with the highest likelihood. 
         ************************************************************************/
        void refineGeneTreeSPRs(map<string, string> params);
        
        
        void refineGeneTreeSPRs2(map<string, string> params);
        
        
        /************************************************************************
         * Tries all NNIs, and accepts NNIs that improve the likelihood as soon as
         * they have been tried.
         ************************************************************************/
        void refineGeneTreeNNIs(map<string, string> params, unsigned int verbose = 0);
        
        /************************************************************************
         * Tells if the gene family is single copy (1 gene per sp)
         ************************************************************************/
        bool isSingleCopy();
        
        
    };
    
    
} //end of namespace bpp.



#endif
