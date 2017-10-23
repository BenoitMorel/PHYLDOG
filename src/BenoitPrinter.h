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


   struct TreeinfoSortCell {
    pll_unode_t *node;
    std::string str;
    std::map<unsigned int, bool> printLeftFirst;
    friend bool operator< (const TreeinfoSortCell &c1, const TreeinfoSortCell &c2) {
      return c1.str < c2.str;
    }
  };

  std::string getTreeinfoString(pllmod_treeinfo_t *treeinfo, 
      bool printBL=false,
      bool printId=false,
      bool rooted=true,
      LabelConversion conversion=NO_CONVERSION)
  {
    return getUnodeString(treeinfo->root, printBL, printId, rooted, conversion);

  }

  std::string getUnodeString(pll_unode_t *node, 
      bool printBL=false,
      bool printId=false,
      bool rooted=true,
      LabelConversion conversion=NO_CONVERSION)
  {
    std::vector<TreeinfoSortCell> cells(rooted ? 2 : 3);
    if (rooted) {
      cells[0].node = node;
      cells[1].node = node->back;
    } else {
      cells[0].node = node->back;
      cells[1].node = node->next->back;
      cells[2].node = node->next->next->back;
    }
    for (unsigned int i = 0; i < cells.size(); ++i) {
      treeInfoSortRec(cells[i].node, conversion, cells[i].str, cells[i].printLeftFirst);
    }
    std::sort(cells.begin(), cells.end());
    std::ostringstream os;
    os << "(";
    for (unsigned int i = 0; i < cells.size(); ++i) {
      printTreeinfoRec(cells[i].node, printBL, printId, rooted, conversion, cells[i].printLeftFirst, os);
      if (i < cells.size() - 1) {
        os << ",";
      }
    }
    os << ");";
    return os.str();
  }


  std::string getBPPNodeString(const bpp::Node *node, 
      bool printBL=false,
      bool printId=false,
      LabelConversion conversion=NO_CONVERSION)
  {
    std::map<unsigned int, std::vector<unsigned int> > printLeftFirst;
    std::string leftStr;
    
    bppNodeSortRec(node, conversion, leftStr, printLeftFirst);
  

    std::ostringstream os;
    printBPPNodeRec(node, printBL, printId, conversion, printLeftFirst, os);
    os << ";";
    return os.str();
  }

void printUnode(pll_unode_t *node) {
  if (!node) {
    std::cout << "(null)" << std::endl;
  } else if (!node->next) {
    std::cout << "(";
    std::cout << node->node_index;
    std::cout << ")";
    std::cout << std::endl;
  } else {
    std::cout << "(";
    std::cout << node->node_index;
    std::cout << ",";
    std::cout << node->next->node_index;
    std::cout << ",";
    std::cout << node->next->next->node_index;
    std::cout << ")";
    std::cout << std::endl;
  }
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

  
  struct StringSortCell {
    unsigned int sonPos;
    std::string str;
    friend bool operator< (const StringSortCell &c1, const StringSortCell &c2) {
      return c1.str < c2.str;
    }
  };

  void bppNodeSortRec(const bpp::Node *node, 
      LabelConversion conversion,
      std::string &firstLeaf,
      std::map<unsigned int, std::vector<unsigned int> > &printLeftFirst)
  {
    if (node->isLeaf()) {
      firstLeaf = getBPPLabel(node, conversion);
      return;
    }
    std::vector<StringSortCell> cells;
    std::vector< bpp::Node * > & sons = ((bpp::Node*)node)->getSons(); //horrible const cast
    for (unsigned int i = 0; i < sons.size(); ++i) {
      StringSortCell cell; 
      bppNodeSortRec(sons[i], conversion, cell.str , printLeftFirst);
      cell.sonPos = i;
      cells.push_back(cell);
    }
    std::sort(cells.begin(), cells.end());
    firstLeaf = cells[0].str;
    std::vector<unsigned int> sonsPos;
    for (unsigned int i = 0; i < cells.size(); ++i)
      sonsPos.push_back(cells[i].sonPos);
    printLeftFirst[node->getId()] = sonsPos;
  }

  void printBPPNodeRec(const bpp::Node *node,
      bool printBL,
      bool printId,
      LabelConversion conversion,
      std::map<unsigned int, std::vector<unsigned int> > &printLeftFirst,
      std::ostream &os)
  {
    if (node->isLeaf())
    {
      os << getBPPLabel(node, conversion);
    } else {
      os << "(";
      std::vector<unsigned int> &sonsPos = printLeftFirst[node->getId()];
      for (unsigned int i = 0; i < sonsPos.size(); ++i) {
        printBPPNodeRec(node->getSon(sonsPos[i]), printBL, printId, conversion, printLeftFirst, os);
        if (i < sonsPos.size() - 1)
          os << ",";
      }
      os << ")";
    }
    if (printBL) {
      if (node->hasDistanceToFather()) {
        if (node->getFather()->hasFather()) 
          os << ":" << node->getDistanceToFather();
        else
          os << ":" << node->getDistanceToFather() * 2; // todobenoit this should be + distance uncle
      }
    }
    if (printId) {
      os << ":" << node->getId();
    }
  }
  
  void printTreeinfoRec(pll_unode_t *node,
      bool printBL,
      bool printId,
      bool root,
      LabelConversion conversion,
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
          printBL, printId, false, conversion, printLeftFirst, os);
      os << ",";
      printTreeinfoRec(!leftFirst ? node->next->back : node->next->next->back, 
          printBL, printId, false, conversion, printLeftFirst, os);
      os << ")";
    }
    if (printBL) {
      if (root) {
        os << ":" << node->length / 2;
      } else {
        os << ":" << node->length;
      }
    }
    if (printId) {
      os << ":" << node->node_index;
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
  
  std::string getBPPLabel(const bpp::Node *node,
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
