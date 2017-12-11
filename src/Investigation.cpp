
extern "C" {
#include <pll/pll.h>
#include <pllmodules/pll.h>
#include <pllmodules/pllmod_algorithm.h>
#include <pllmodules/pll_binary.h>
#include <pllmodules/pll_msa.h>
#include <pllmodules/pll_optimize.h>
#include <pllmodules/pll_tree.h>
#include <pllmodules/pllmod_util.h>
}

#include <iostream>

#define TOLERANCE 0.5

struct PLLData {
  pllInstanceAttr attributes;
  pllInstance *instance;
  pllAlignmentData *alignmentData;
  pllNewickTree *newick;
  partitionList *partitions;
  pllQueue *partitionInfo;
};

void execute_pll(const char *newick, const char *part, const char *msa)
{
  PLLData plldata;
  plldata.attributes.rateHetModel     = PLL_GAMMA;
  plldata.attributes.fastScaling      = PLL_TRUE;
  plldata.attributes.saveMemory       = PLL_FALSE;
  plldata.attributes.useRecom         = PLL_FALSE;
  plldata.attributes.randomNumberSeed = 0xDEADBEEF;
  plldata.attributes.numberOfThreads  = 8;            /* This only affects the pthreads version */
  plldata.instance = pllCreateInstance (&plldata.attributes);
  
  plldata.alignmentData = pllParseAlignmentFile (PLL_FORMAT_FASTA, msa);
  if (!plldata.alignmentData) {
    std::cerr << "Error in pllParseAlignmentFile " << msa << std::endl;
    return;
  }
  
  
  plldata.partitionInfo = pllPartitionParse (part);

  if (!pllPartitionsValidate (plldata.partitionInfo, plldata.alignmentData)) {
    std::cerr << "Error in pllPartitionParse " << part << std::endl;
    return;
  }

  plldata.partitions = pllPartitionsCommit (plldata.partitionInfo, plldata.alignmentData);
  plldata.partitions->partitionData[0]->frequencies = 0;
  pllQueuePartitionsDestroy (&plldata.partitionInfo);
  pllAlignmentRemoveDups (plldata.alignmentData, plldata.partitions);

  plldata.newick = pllNewickParseFile (newick);
  if (!plldata.newick || !pllValidateNewick (plldata.newick)) {
    std::cerr << "Error in pllNewickParseFile " << newick << std::endl;
    return;
  }

  
  pllTreeInitTopologyNewick (plldata.instance, plldata.newick, PLL_FALSE);
  if (!pllLoadAlignment (plldata.instance, plldata.alignmentData, plldata.partitions)) {
    std::cerr << "Error in pllLoadAlignment" << std::endl;
    return;
  }
  pllNewickParseDestroy(&plldata.newick);

  pllInitModel(plldata.instance, plldata.partitions);
  
  std::cout << "PLL initial likelihood: " << plldata.instance->likelihood << std::endl;
  pllOptimizeModelParameters(plldata.instance, plldata.partitions, TOLERANCE);
  std::cout << "PLL final likelihood: " << plldata.instance->likelihood << std::endl;
}


void execute_libpll2(const char *newick, const char *part, const char *msa)
{

}

int main(int args, char ** argv)
{
  const char *newick = "data_investigation/plltree.newick";
  const char *msa = "data_investigation/pllmsa.fasta";
  const char *part = "data_investigation/partition.txt";
  execute_pll(newick, part,  msa);
  std::cout << "end of the story" << std::endl;
  return 0;
}



