#ifndef BENOITPRINTER
#define BENOITPRINTER

#include<string>
#include<map>

#include<Bpp/Phyl/Node.h>
#include<Bpp/Phyl/TreeTemplate.h>
#include<Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>


extern "C" {
#include <pll/pll.h>
#include <pllmodules/pll.h>
#include <pllmodules/pll_tree.h>
#include <pllmodules/pllmod_util.h>
}

template<class T>
void print(const char* msg, T *data, unsigned int size)
{
  std::cout << msg << " ";
  for (unsigned int i = 0; i < size; ++i) 
    std::cout << data[i] << " ";
  std::cout << std::endl;
}

template<class T>
void print_v(vector<T> v)
{
  for (unsigned int i = 0; i < v.size(); ++i)
    std::cout << v[i] << " ";
  std::cout << std::endl;
}

enum LabelConversion {
  NO_CONVERSION,
  STRICT_TO_REAL,
  REAL_TO_STRICT
};

class BenoitPrinter {
public:
  BenoitPrinter() {}
  BenoitPrinter(std::map<std::string,std::string> &realToStrict,
      std::map<std::string,std::string> &strictToReal):
    realToStrict(realToStrict), strictToReal(strictToReal) {}


  std::string getTreeinfoString(pllmod_treeinfo_t *treeinfo, 
      LabelConversion conversion=NO_CONVERSION,
      bool printBL=false)
  {
    std::map<unsigned int, bool> printLeftFirst;
    std::string leftStr, rightStr;
    pll_unode_t *root = treeinfo->root;
    treeInfoSortRec(root, conversion, leftStr, printLeftFirst);
    treeInfoSortRec(root->back, conversion, rightStr, printLeftFirst);
  

    std::ostringstream os;
    os << "(";
    printTreeinfoRec((leftStr < rightStr) ? root : root->back, 
        conversion, printBL, printLeftFirst, os);
    os << ",";
    printTreeinfoRec(!(leftStr < rightStr) ? root : root->back, 
        conversion, printBL, printLeftFirst, os);
    os << ");";
    return os.str();
  }

  std::string getBPPNodeString(bpp::Node *node, 
      LabelConversion conversion=NO_CONVERSION,
      bool printBL=false)
  {
    std::map<unsigned int, bool> printLeftFirst;
    std::string leftStr, rightStr;
    
    bppNodeSortRec(node, conversion, leftStr, printLeftFirst);
  

    std::ostringstream os;
    printBPPNodeRec(node, conversion, printBL, printLeftFirst, os);
    os << ";";
    return os.str();
  }



private:
  std::map<std::string,std::string> realToStrict;
  std::map<std::string,std::string> strictToReal;

  void treeInfoSortRec(pll_unode_t *node, 
      LabelConversion conversion,
      std::string &firstLeaf,
      std::map<unsigned int, bool> &printLeftFirst)
  {
    if (!node->next) {
      firstLeaf = getLabel(node, conversion);
      return;
    }
    std::string leftStr;
    std::string rightStr;
    treeInfoSortRec(node->next->back, conversion, leftStr, printLeftFirst);
    treeInfoSortRec(node->next->next->back, conversion, rightStr, printLeftFirst);
    if (leftStr < rightStr) {
      firstLeaf = leftStr;
      printLeftFirst[node->node_index] = true;
    } else {
      firstLeaf = rightStr;
      printLeftFirst[node->node_index] = false;
    }
  }

  void bppNodeSortRec(bpp::Node *node, 
      LabelConversion conversion,
      std::string &firstLeaf,
      std::map<unsigned int, bool> &printLeftFirst)
  {
    if (node->isLeaf()) {
      firstLeaf = getBPPLabel(node, conversion);
      return;
    }
    std::string leftStr;
    std::string rightStr;
    bppNodeSortRec(node->getSon(0), conversion, leftStr, printLeftFirst);
    bppNodeSortRec(node->getSon(1), conversion, rightStr, printLeftFirst);
    if (leftStr < rightStr) {
      firstLeaf = leftStr;
      printLeftFirst[node->getId()] = true;
    } else {
      firstLeaf = rightStr;
      printLeftFirst[node->getId()] = false;
    }
  }

  void printBPPNodeRec(bpp::Node *node,
      LabelConversion conversion,
      bool printBL,
      std::map<unsigned int, bool> printLeftFirst,
      std::ostream &os)
  {
    if (node->isLeaf())
    {
      os << getBPPLabel(node, conversion);
    } else {
      os << "(";

      bool leftFirst = printLeftFirst[node->getId()];
      printBPPNodeRec(node->getSon(leftFirst ? 0 : 1), conversion, printBL, printLeftFirst, os);
      os << ",";
      printBPPNodeRec(node->getSon(!leftFirst ? 0 : 1), conversion, printBL, printLeftFirst, os);
      os << ")";
    }
    if (printBL) {
      os << ":" << node->getDistanceToFather();
    }
  }
  
  void printTreeinfoRec(pll_unode_t *node,
      LabelConversion conversion,
      bool printBL,
      std::map<unsigned int, bool> printLeftFirst,
      std::ostream &os)
  {
    if (!node->next)
    {
      os << getLabel(node, conversion);
    } else {
      os << "(";

      bool leftFirst = printLeftFirst[node->node_index];
      printTreeinfoRec(leftFirst ? node->next->back : node->next->next->back, 
          conversion, printBL, printLeftFirst, os);
      os << ",";
      printTreeinfoRec(!leftFirst ? node->next->back : node->next->next->back, 
          conversion, printBL, printLeftFirst, os);
      os << ")";
    }
    if (printBL) {
      os << ":" << node->length;
    }
  }

  std::string getLabel(pll_unode_t *node,
      LabelConversion conversion)
  {
    if (conversion == NO_CONVERSION)
      return node->label;
    else if (conversion == STRICT_TO_REAL)
      return strictToReal[node->label];
    else 
      return realToStrict[node->label];
  }
  
  std::string getBPPLabel(bpp::Node *node,
      LabelConversion conversion)
  {
    if (conversion == NO_CONVERSION)
      return node->getName();
    else if (conversion == STRICT_TO_REAL)
      return strictToReal[node->getName()];
    else 
      return realToStrict[node->getName()];
  }


};





#endif
