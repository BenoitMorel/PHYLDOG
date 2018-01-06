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

#include <iostream>
#include <string>
#include <fstream>
#include "ReconciliationTools.h"
#include "FastReconciliationTools.h" 
#include <ctime>

void parseLinkFile(const std::string &filename, std::map<std::string, std::string > &genesToSpecies)
{
  std::ifstream input(filename);
  string line;
  while (std::getline(input, line)) {
    std::istringstream lineStream(line);
    std::string specie;
    std::string gene;
    std::getline(lineStream, specie, ':');
    while (std::getline(lineStream, gene, ';')) {
      genesToSpecies[gene] =  specie;
    }
  }
}

void initRates(std::vector<double> &lossRates, std::vector<double> &duplicationRates, unsigned int size) 
{
  lossRates.clear();
  lossRates.resize(size);
  duplicationRates.clear();
  duplicationRates.resize(size);
  for (unsigned int i = 0; i < size; ++i) {
    lossRates[i] = double(i) / double (2 * size) ;
    duplicationRates[i] = double(i) / double (size);
  }
}

void root(TreeTemplate<Node> *tree)
{
  if (tree->isRooted())
    return;
  std::vector<int> innerNodes = tree->getInnerNodesId();
  tree->newOutGroup(innerNodes[innerNodes.size()/2]);
}

int main(int args, char ** argv)
{
  std::cout << "Bench reconciliation..." << std::endl;
  
  if (args != 4) {
    std::cerr << "Error: syntax is " << std::endl;
    std::cerr << "./benchreconciliation speciesTree genesTree speciesToGenes" << std::endl;
  }

  std::string speciesTreeFilename = argv[1];
  std::string geneTreeFilename = argv[2];
  std::string linksFilename = argv[3];


  Newick newick;

  TreeTemplate<Node> *speciesTree = newick.read(speciesTreeFilename);
  TreeTemplate<Node> *geneTree = newick.read(geneTreeFilename);
  root(speciesTree);
  root(geneTree);
  breadthFirstreNumber(*speciesTree);
  std::map<std::string, std::string > genesToSpecies;
  std::map<std::string, int > speciesIDs = computeSpeciesNamesToIdsMap(*speciesTree);
  std::vector<double> lossRates;
  std::vector<double> duplicationRates;
  int MLindex = 0;
  std::vector<int> num0lineages(speciesTree->getNumberOfNodes(), 0);
  std::vector<int> num1lineages(speciesTree->getNumberOfNodes(), 0);
  std::vector<int> num2lineages(speciesTree->getNumberOfNodes(), 0);
  std::set<int> nodesToTryInNNISearch;

  parseLinkFile(linksFilename, genesToSpecies);
  //initSpeciesIDs(speciesTree, speciesIDs);
  initRates(lossRates, duplicationRates, speciesTree->getNumberOfNodes());
 
  double ll = 0;

  clock_t begin = clock();

  /////////////////////////////////////////////
  ////////////// BEEENCH //////////////////////
  /////////////////////////////////////////////
  
  std::vector<int> nodes = speciesTree->getNodesId();
  const unsigned int ITERATIONS = 30;
  for (unsigned int i = 0; i < nodes.size(); ++i) {
    speciesTree->newOutGroup(nodes[i]);
    for (unsigned int iteration = 0; iteration < ITERATIONS; ++iteration) {
      FastReconciliationTools rc(speciesTree, 
        geneTree, 
        genesToSpecies,
        speciesIDs, 
        lossRates,
        duplicationRates,
        num0lineages,
        num1lineages,
        num2lineages,
        nodesToTryInNNISearch,
        false);
      ll = rc.findMLReconciliationDR(MLindex);
      /*
      ll = findMLReconciliationDR(speciesTree, 
        geneTree, 
        genesToSpecies,
        speciesIDs, 
        lossRates,
        duplicationRates,
        MLindex,
        num0lineages,
        num1lineages,
        num2lineages,
        nodesToTryInNNISearch,
        true);
        */
      if (iteration == ITERATIONS - 1) 
        std::cout << " Reconciliation ll: " << ll << std::endl;
    }
  }
  
  /////////////////////////////////////////////
  ////////////// end BENCH /////////////////////
  /////////////////////////////////////////////
  
  clock_t end = clock();
  
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  
  std::cout << "ELAPSED " << elapsed_secs << std::endl;

  return 0;
}




