#ifndef __QCSTATS_H__
#define __QCSTATS_H__

#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#include <stdint.h>
#include "Sequence.h"
#include "SamFile.h"
//#include "SamFileHeader.h"
//#include "SamRecord.h"
#include <vector>
#include <map>
#include "GenomeSequence.h"
#include "GCContent.h"
#include "QSamFlag.h"
#include "DepthCounter.h"

using namespace std;

class QCStats
{
  void constructorClear();
public:
  GenomeSequence *referencegenome;
  std::vector<bool> * regionIndicator;
  std::vector<bool> * dbSNP;
  GCContent *GC;
  // ReadDepth *depthVec;
  DepthCounter* depthVec;
  uint32_t coveredGCFreq[101];
  uint64_t depthTotalVsGC[101];
  double   depthVsGC[101];
  double   depthVsGC_norm[101];
  
  std::map<int, uint64_t> matchCountByQual;
  std::map<int, uint64_t> misMatchCountByQual;
  std::map<int, double> misMatchRateByQual;
  std::map<int, uint64_t> *baseCountByQualByCycle;
  std::vector<uint64_t> baseQ20CountByCycle;
  std::vector<uint64_t> baseReportedQ20CountByCycle;
  std::map<int, uint64_t> qualCount;
  std::map<int, uint64_t> Q20QualScores;
  std::vector<int> qual;
  std::vector<bool> *genomePosCovered;
  std::map<int32_t, uint64_t> insertSize;
  std::map<int, uint64_t> cycles;
  std::map<int, uint64_t> depthDist;
  static const int depthThreshold[7];
  uint32_t depthDistribution[7]; // count >=1, >=5, ... sites
      
  int size;	 	 	 // Max length of reads 
  int size_reserved;	 	 // Size of memory allocated
  uint64_t *matchCountByCycle;
  uint64_t *misMatchCountByCycle;
  double  *misMatchRateByCycle;
  uint64_t totalMappedBases;	// calculated but never used
  uint64_t regionMappedBases;	// used at BamQC.cpp line 570
  double coverage;  
  uint64_t matchMatrix[6][6];
  uint64_t **baseCountByCycle;
  double **baseCompositionByCycle;
  double baseComposition[6];

  //general stats
  uint64_t nReads;
  uint64_t nUnMapped;
  uint64_t nUnMapped_Filter;
  uint64_t nZeroMapQual;
  uint64_t nLT10MapQual;
  uint64_t nReadsMapped2TargetRegions;
  uint64_t nQ20;
  double   pQ20; //proportion of Q20 bases
  uint64_t nDup;
  uint64_t nQCFail;
  uint64_t nPaired;
  uint64_t nProperPaired;
  uint64_t nSecondary;
  double genomeCoverage;
  double gcBiasStat;
  double insertSize_mean;
  double insertSize_var;
  int insertSize_mode;
  int insertSize_median;
  int MAX_ISIZE_ALLOWED;
  uint64_t nWarnings;

public:
  QCStats();
  QCStats(int);
  ~QCStats();
 
  void Init(int);  
  void ReAllocateMemory(int);
  void SetReferenceGenome(GenomeSequence *ref) { referencegenome = ref; };
  void SetdbSNPIndicator(vector<bool> *db){ dbSNP = db; };
  void SetGenomePositionCoveredIndicator(vector<bool> *gpc){ genomePosCovered = gpc;};
  void SetGCContent(GCContent *gc) { GC = gc; }
  void SetDepth(DepthCounter *depth) {depthVec = depth;}
  void SetRegionIndicator(std::vector<bool> *region){regionIndicator = region;};
  void CalcGenomeCoverage(uint64_t);
  void CalcMisMatchRateByCycle();
  double CalcMisMatchRateByCycle_MEAN();
  void CalcMisMatchRateByQual();
  double CalcMisMatchRateByQual_MSE();
  void CalcQ20Bases();
  void CalcQ20BasesByCycle();
  void CalcReportedQ20BasesByCycle();
  void CalcBaseComposition();
  void CalcDepthGC(GCContent &gc, std::vector<bool> &);
  void CalcDepthDist();
  void CalcGCBias(int s, int e);
  void CalcInsertSize_mean();
  void CalcInsertSize_var();
  void CalcInsertSize_mode();
  void CalcInsertSize_median();
  void PrintSamRecord(SamRecord &);
  void ReportWarningCount();
  void UpdateStats(SamRecord &, QSamFlag &filter, double, std::map<int, int> &, std::vector<std::string>& readGroups);
};

#endif
